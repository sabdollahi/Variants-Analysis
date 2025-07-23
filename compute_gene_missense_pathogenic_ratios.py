'''
Author: Sina Abdollahi
Release Date: 23.07.2025
Title: Extract Pathogenic Variants Ratio and Missense Variants Ratio

Description: Given a name of the directory containing the detail of the variants of various genes and compute the ratio of number of pathogenic;likely pathogenic over the number of missense variants, the ratio of number of pathogenic;likely pathogenic over the length of the protein, and the ratio of the number of missense variants over the length of the protein.
The results are saved into a new CSV file called "Protein_Prioritization_Values.csv".
'''


import pandas as pd
import numpy as np

gene_variants_dir = "VARINATS-DETAIL-FILES-DIRECTORY"
data = {}
data["Protein Name"] = []
data["num_misse/prot_len"] = []
data["num_patho/prot_len"] = []
data["num_patho/num_misse"] = []

ii = 1
print(ii, end="\t")
for gene_variants in os.listdir(gene_variants_dir):
    if(gene_variants.find('variants.csv') != -1):
        if((ii % 500) == 0):
            df = pd.DataFrame(data)
            df.to_csv('Protein_Prioritization_Values.csv', index=False)
            print(ii, end="\t")
        gene_name = gene_variants.split("_")[0]
        print(gene_name, end=' ')
        
        variants_df = pd.read_csv(f'{gene_variants_dir}/{gene_variants}')
        max_pos = variants_df['Position'].max()
        min_pos = variants_df['Position'].min()
        #Protein approximate length
        prot_len = (max_pos - min_pos)/3
        
        # We only keep the variants that are Missense Variants
        variants_df = variants_df[variants_df['VEP Consequence']=='missense_variant']
        #Number of Missense Variants
        num_misse = len(variants_df)
        variants_df = variants_df[['Variant ID', 'Position', 'Coding Change', 'Protein Change', 'Clinical Significance']]
        pathogenicity_group = ["Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic; other", "Pathogenic; Affects", 
                               "Pathogenic; association", "Pathogenic; drug response", "Pathogenic; risk factor", 
                               "Pathogenic; confers sensitivity", "Pathogenic/Likely pathogenic; other", 
                               "Pathogenic/Likely pathogenic; association", "Pathogenic/Likely pathogenic; risk factor", 
                               "Pathogenic/Likely pathogenic/Pathogenic, low penetrance", "Pathogenic/Likely pathogenic/Pathogenic, low penetrance; other", "Pathogenic/Pathogenic, low penetrance; other; risk factor", "Pathogenic/Pathogenic, low penetrance; other", 
                               "Pathogenic/Likely pathogenic/Likely risk allele", "Pathogenic/Likely risk allele", "Likely pathogenic", "Likely pathogenic; other", "Likely pathogenic; Affects", "Likely pathogenic; risk factor", 
                               "Likely pathogenic; drug response", "Likely pathogenic/Likely risk allele", "Likely pathogenic, low penetrance"]
        # Number of Pathogenic and likely pathogenic variants
        num_patho = len(variants_df[variants_df['Clinical_score'].isin(pathogenicity_group)])
        nm_pl = 0
        np_pl = 0
        np_nm = 0
        if(prot_len != 0):
            nm_pl = num_misse/prot_len
            np_pl = num_patho/prot_len
        if(num_misse != 0):
            np_nm = num_patho/num_misse
        
        data["Protein Name"].append(gene_name)
        data["num_misse/prot_len"].append(nm_pl)
        data["num_patho/prot_len"].append(np_pl)
        data["num_patho/num_misse"].append(np_nm)
        
df = pd.DataFrame(data)
df.to_csv('Protein_Prioritization_Values.csv', index=False)
