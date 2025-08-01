# Variants' Extraction and Analysis Project

This repository contains code for extracting and analyzing variants data of different genes.

## Conda environment

```
conda create -n variants_info python=3.11
conda activate variants_info

pip install gql[all] pandas numpy
```

## Data Extraction from GnomAD graphQL service
### Variants' Overall Information Extraction

Given a list of gene symbols in an Excel file with header of "Gene_Symbol" and we extract the variants' information such as Expected Loss of function, Observed Loss of function, Observed/Expected Loss of function, Probability of being loss-of-function (pLI), Loss of function Z score, Observed Missense Variants, Expected Missense Variants, Observed synonymous variants, Expected synonymous variants, and their corresponding Z scores.

Simply replace the name of your Excel file to "GENE-LIST-FILE.xlsx" in line 19 and run the following code:

```
python gnomad_constraint_metrics_extractor.py
```

### Details of the variants of various genes Extraction

To extract the details of each variant in the genes, again, given a list of gene symbols in an Excel file with header of "Gene_Symbol" and we extract the detail of the variants such as Variant ID, Chromosome, Position, RS IDs, HGVS Consequence, Coding Change, Protein Change, VEP Consequence, Clinical Significance, ClinVar Variation ID, Allele Count/Frequency/Number for genome and exome, Homozygote/Hemizygote Count for genome and exome, CADD score, SpliceAI Ds Max score, Pangolin Largest Ds score, and PhyLop score. 
For each gene, we generate a new file with suffix of "_variants.csv" containing the details of the variants.

Simply replace the name of your Excel file to "GENE-LIST-FILE.xlsx" in line 26 and run the following code:

```
python gnomad_variants_details_extractor.py
```

### Compute the Ratio of Pathogenicity and Missense of the genes


If you are interested in to compute the ratio of number of pathogenic;likely pathogenic over the number of missense variants, the ratio of number of pathogenic;likely pathogenic over the length of the protein, and the ratio of the number of missense variants over the length of the protein, you can use the following command:

```
python compute_gene_missense_pathogenic_ratios.py
```

**Note:** Remember to replace the name of the directory in line No. 14, and also remember that the name of the files in the directory must have suffix of "variants.csv".
