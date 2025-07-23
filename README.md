# Variants' Extraction and Analysis Project

This repository contains code for extracting and analyzing variants data of different genes.

## Conda environment

```
conda create -n variants_info python=3.11
conda activate variants_info

pip install gql[all] pandas
```

## Data Extraction from GnomAD graphQL service
### Variants Information Extraction

Given a list of gene symbols in an Excel file with header of "Gene_Symbol" and we extract the variants' information such as Expected Loss of function, Observed Loss of function, Observed/Expected Loss of function, Probability of being loss-of-function (pLI), Loss of function Z score, Observed Missense Variants, Expected Missense Variants, Observed synonymous variants, Expected synonymous variants, and their corresponding Z scores.

Simply replace the name of your Excel file to "GENE-LIST-FILE.xlsx" in line 19 and run the following code:

```
python gnomad_variants_extractor.py
```
