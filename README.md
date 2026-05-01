# V2G-Score (Variant-to-Gene)

V2G-Score is a Python prototype tool for prioritizing potential causal genes for GWAS hits. It integrates distance-based mapping with eQTL evidence to rank candidate genes near a trait-associated genetic variant (SNP).

## Methodology

The tool assigns candidate genes to a given GWAS SNP by calculating a **Global Priority Score**, combining two independent pieces of evidence:

1. **Distance Score (30% weight):** 
   - Evaluates any gene within a 500kb window of the SNP.
   - Calculates a score using exponential decay: `exp(-d / 100,000)`, where `d` is the absolute distance in base pairs from the SNP to the gene's Transcription Start Site (TSS).
2. **eQTL Score (70% weight):**
   - Incorporates functional evidence via expression quantitative trait loci (eQTLs).
   - Normalizes the eQTL p-value to a 0-1 scale: `min(1, -log10(P) / 10)`.

**Global Priority Score = 0.3 * Distance Score + 0.7 * eQTL Score**

For each GWAS SNP, the tool ranks the mapped genes based on this priority score and outputs the top 3 candidates.

## Requirements

- Python 3.x
- `pandas`
- `numpy`

## Input Data Format

The script requires standard TSV files as input.

- **GWAS Data** (`--gwas`): Must contain columns `rsid`, `chromosome`, `position`, `p_value`.
- **Gene Map Data** (`--gene_map`): Must contain columns `gene_id`, `chromosome`, `tss_position`.
- **eQTL Data** (`--eqtl`) (Optional): Must contain columns `rsid`, `gene_id`, `p_value`, `tissue`.

## Usage

You can run the script directly from the command line:

```bash
python v2g_score.py --gwas data/gwas.tsv \
                    --gene_map data/gene_map.tsv \
                    --eqtl data/eqtl.tsv \
                    --out results/v2g_output.csv
```

### Arguments

- `--gwas`: Path to the GWAS TSV file (Required).
- `--gene_map`: Path to the Gene Map TSV file (Required).
- `--eqtl`: Path to the eQTL TSV file (Optional). If omitted, only distance scores will be used.
- `--out`: Path to save the final ranked CSV output file (Required).

## Running Tests

Unit tests are included to ensure scoring and prioritization logic remains robust. To run the tests:

```bash
python -m unittest test_v2g_score.py
```
