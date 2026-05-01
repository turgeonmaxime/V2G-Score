import pandas as pd
import numpy as np

class V2GAnnotator:
    def __init__(self):
        self.gwas = None
        self.eqtl = None
        self.gene_map = None
        self.results = None

    def load_data(self, gwas_df, eqtl_df, gene_map_df):
        """
        Load input data into the annotator.
        """
        self.gwas = gwas_df.copy()
        self.eqtl = eqtl_df.copy() if eqtl_df is not None else None
        self.gene_map = gene_map_df.copy()
        
    def map_distance(self):
        """
        Calculate distance score for genes within 500kb of a SNP.
        Distance Score = exp(-d / 100,000)
        """
        if self.gwas is None or self.gene_map is None:
            raise ValueError("GWAS and Gene Map data must be loaded.")

        # Cross join GWAS and Gene Map on chromosome
        merged = pd.merge(self.gwas, self.gene_map, on='chromosome', suffixes=('_gwas', '_gene'))
        
        # Calculate absolute distance to TSS
        merged['distance'] = np.abs(merged['position'] - merged['tss_position'])
        
        # Filter to genes within 500kb
        within_500kb = merged[merged['distance'] <= 500000].copy()
        
        # Calculate distance score: exp(-d / 100,000)
        within_500kb['distance_score'] = np.exp(-within_500kb['distance'] / 100000)
        
        return within_500kb[['rsid', 'gene_id', 'distance', 'distance_score']]
        
    def map_eqtl(self):
        """
        Calculate eQTL score based on negative log10 p-value.
        eQTL Score = min(1, -log10(P) / 10)
        """
        if self.eqtl is None or self.eqtl.empty:
            return pd.DataFrame(columns=['rsid', 'gene_id', 'eqtl_score'])

        eqtl_scores = self.eqtl.copy()
        
        # Handle 0 p-values to avoid log(0) warnings
        min_p = eqtl_scores.loc[eqtl_scores['p_value'] > 0, 'p_value'].min()
        if pd.isna(min_p):
            min_p = 1e-300
        eqtl_scores['p_value'] = np.where(eqtl_scores['p_value'] == 0, min_p, eqtl_scores['p_value'])
        
        # Calculate eQTL score: min(1, -log10(P) / 10)
        eqtl_scores['eqtl_score'] = np.minimum(1, -np.log10(eqtl_scores['p_value']) / 10)
        
        # If multiple tissues exist for the same rsid-gene_id pair, take the max score
        max_eqtl = eqtl_scores.groupby(['rsid', 'gene_id'])['eqtl_score'].max().reset_index()
        
        return max_eqtl
        
    def generate_report(self, output_file=None):
        """
        Generate global priority score and output top 3 candidate genes for each GWAS rsid.
        Priority Score = 0.3 * Distance score + 0.7 * eQTL score
        """
        if self.gwas is None or self.gene_map is None:
            raise ValueError("Data not loaded. Run load_data() first.")
            
        dist_df = self.map_distance()
        eqtl_df = self.map_eqtl()
            
        # Outer join to combine distance and eQTL evidence
        combined = pd.merge(dist_df, eqtl_df, on=['rsid', 'gene_id'], how='outer')
        
        # Fill missing scores with 0
        combined['distance_score'] = combined['distance_score'].fillna(0)
        combined['eqtl_score'] = combined['eqtl_score'].fillna(0)
        
        # Calculate Global Priority Score
        combined['priority_score'] = 0.3 * combined['distance_score'] + 0.7 * combined['eqtl_score']
        
        # Filter to pairs where we have at least one piece of evidence
        combined = combined[combined['priority_score'] > 0]
        
        # Sort and rank (descending priority score)
        combined = combined.sort_values(['rsid', 'priority_score'], ascending=[True, False])
        
        # Get top 3 per rsid
        top_candidates = combined.groupby('rsid').head(3)
        
        self.results = top_candidates
        
        if output_file:
            top_candidates.to_csv(output_file, sep='\t', index=False)
            
        return top_candidates

def generate_mock_data():
    """
    Generate dummy data for 5 SNPs and 10 Genes.
    """
    np.random.seed(42)
    
    # 5 SNPs on chr1
    gwas_data = {
        'rsid': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
        'chromosome': ['chr1'] * 5,
        'position': [1000000, 1500000, 2000000, 2500000, 3000000],
        'p_value': [5e-8, 1e-7, 5e-9, 2e-6, 1e-10]
    }
    gwas_df = pd.DataFrame(gwas_data)
    
    # 10 Genes on chr1
    gene_data = {
        'gene_id': [f'GENE{i}' for i in range(1, 11)],
        'chromosome': ['chr1'] * 10,
        'tss_position': [
            1050000,  # 50kb from rs1 -> distance score: 0.606
            1400000,  # 100kb from rs2 -> distance score: 0.367
            1550000,  # 50kb from rs2 -> distance score: 0.606
            2010000,  # 10kb from rs3 -> distance score: 0.904
            2600000,  # 100kb from rs4 -> distance score: 0.367
            3000000,  # 0kb from rs5 -> distance score: 1.0
            4000000,  # Far
            4500000,  # Far
            5000000,  # Far
            5500000   # Far
        ]
    }
    gene_map_df = pd.DataFrame(gene_data)
    
    # eQTL evidence
    eqtl_data = {
        'rsid': ['rs1', 'rs1', 'rs2', 'rs3', 'rs4'],
        'gene_id': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE10'], # GENE10 is far away from rs4 but has eQTL
        'p_value': [1e-15, 1e-5, 1e-20, 1e-8, 1e-2],
        'tissue': ['Liver', 'Liver', 'Brain', 'Blood', 'Blood']
    }
    eqtl_df = pd.DataFrame(eqtl_data)
    
    return gwas_df, eqtl_df, gene_map_df

if __name__ == "__main__":
    print("Generating mock data...")
    gwas, eqtl, gene_map = generate_mock_data()
    
    print("Initializing V2GAnnotator...")
    annotator = V2GAnnotator()
    annotator.load_data(gwas, eqtl, gene_map)
    
    output_path = "v2g_results.tsv"
    print(f"Generating report and saving to {output_path}...")
    results = annotator.generate_report(output_path)
    
    print("\nTop 3 Candidate Genes per GWAS Hit:")
    print(results.to_string(index=False))
