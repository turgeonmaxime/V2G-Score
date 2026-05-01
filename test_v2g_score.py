import unittest
import pandas as pd
import numpy as np
from v2g_score import V2GAnnotator

class TestV2GAnnotator(unittest.TestCase):
    def setUp(self):
        np.random.seed(42)
        
        # 5 SNPs on chr1
        gwas_data = {
            'rsid': ['rs1', 'rs2', 'rs3', 'rs4', 'rs5'],
            'chromosome': ['chr1'] * 5,
            'position': [1000000, 1500000, 2000000, 2500000, 3000000],
            'p_value': [5e-8, 1e-7, 5e-9, 2e-6, 1e-10]
        }
        self.gwas_df = pd.DataFrame(gwas_data)
        
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
        self.gene_map_df = pd.DataFrame(gene_data)
        
        # eQTL evidence
        eqtl_data = {
            'rsid': ['rs1', 'rs1', 'rs2', 'rs3', 'rs4'],
            'gene_id': ['GENE1', 'GENE2', 'GENE3', 'GENE4', 'GENE10'], # GENE10 is far away from rs4 but has eQTL
            'p_value': [1e-15, 1e-5, 1e-20, 1e-8, 1e-2],
            'tissue': ['Liver', 'Liver', 'Brain', 'Blood', 'Blood']
        }
        self.eqtl_df = pd.DataFrame(eqtl_data)

    def test_report_generation(self):
        annotator = V2GAnnotator()
        annotator.load_data(self.gwas_df, self.eqtl_df, self.gene_map_df)
        
        results = annotator.generate_report()
        
        # Test that it successfully generates some results
        self.assertIsNotNone(results)
        self.assertFalse(results.empty)
        
        # Test top candidates logic
        rs1_results = results[results['rsid'] == 'rs1']
        self.assertGreater(len(rs1_results), 0)
        self.assertEqual(rs1_results.iloc[0]['gene_id'], 'GENE1')
        
        rs5_results = results[results['rsid'] == 'rs5']
        self.assertGreater(len(rs5_results), 0)
        self.assertEqual(rs5_results.iloc[0]['gene_id'], 'GENE6')
        
        # Ensure no more than 3 candidates per rsid
        counts = results['rsid'].value_counts()
        self.assertTrue(all(count <= 3 for count in counts))
        
    def test_missing_eqtl(self):
        annotator = V2GAnnotator()
        # Pass None for eqtl
        annotator.load_data(self.gwas_df, None, self.gene_map_df)
        
        results = annotator.generate_report()
        
        self.assertIsNotNone(results)
        self.assertFalse(results.empty)
        self.assertIn('eqtl_score', results.columns)
        self.assertTrue((results['eqtl_score'] == 0).all())

if __name__ == '__main__':
    unittest.main()
