import pandas as pd
import numpy as np
import logging
import argparse

logging.basicConfig(level=logging.INFO, format='%(asctime)s - %(levelname)s - %(message)s')

class V2GAnnotator:
    def __init__(self):
        self.gwas = None
        self.eqtl = None
        self.gene_map = None
        self.results = None
        self.logger = logging.getLogger(__name__)

    def load_data(self, gwas_df, eqtl_df, gene_map_df):
        """
        Load input data into the annotator.
        """
        self.gwas = gwas_df.copy()
        self.eqtl = eqtl_df.copy() if eqtl_df is not None else None
        self.gene_map = gene_map_df.copy()
        
        num_snps = self.gwas['rsid'].nunique()
        self.logger.info(f"Loaded data for {num_snps} unique GWAS SNPs.")
        
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
        
        num_mapped = within_500kb['rsid'].nunique()
        self.logger.info(f"Mapped {num_mapped} SNPs to candidate genes via distance (<=500kb).")
        
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
        
        num_mapped = max_eqtl['rsid'].nunique()
        self.logger.info(f"Mapped {num_mapped} SNPs to candidate genes via eQTL evidence.")
        
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
        
        num_prioritized = top_candidates['rsid'].nunique()
        self.logger.info(f"Prioritized candidate genes for {num_prioritized} unique SNPs.")
        
        self.results = top_candidates
        
        if output_file:
            top_candidates.to_csv(output_file, index=False)
            
        return top_candidates

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='V2G-Score: Variant-to-Gene prioritizing tool.')
    parser.add_argument('--gwas', required=True, help='Path to GWAS TSV file')
    parser.add_argument('--eqtl', help='Path to eQTL TSV file (optional)')
    parser.add_argument('--gene_map', required=True, help='Path to Gene Map TSV file')
    parser.add_argument('--out', required=True, help='Path to output CSV file')
    
    args = parser.parse_args()
    
    logging.info("Starting V2G-Score process...")
    
    gwas_df = pd.read_csv(args.gwas, sep='\t')
    eqtl_df = pd.read_csv(args.eqtl, sep='\t') if args.eqtl else None
    gene_map_df = pd.read_csv(args.gene_map, sep='\t')
    
    annotator = V2GAnnotator()
    annotator.load_data(gwas_df, eqtl_df, gene_map_df)
    
    annotator.generate_report(args.out)
    logging.info(f"Report saved to {args.out}")
