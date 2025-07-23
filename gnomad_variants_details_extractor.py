'''
Author: Sina Abdollahi
Release Date: 23.07.2025
Title: Extracting the variants' detailed information from gnomAD 

Description: Given a list of gene symbols in an Excel file with header of "Gene_Symbol" and we extract 
the detailed variants' information such as Variant ID, Chromosome, Position, RS IDs, HGVS Consequence, Coding Change, Protein Change, VEP Consequence, Clinical Significance, ClinVar Variation ID, Allele Count/Frequency/Number for genome and exome, Homozygote/Hemizygote Count for genome and exome, CADD score, SpliceAI Ds Max score, Pangolin Largest Ds score, and PhyLop score. 
For each gene, we generate a new file with suffix of "_variants.csv" containing the details of the variants. 
'''


import asyncio
from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import time
import pandas as pd
import numpy as np

# Define the gnomAD API endpoint
transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api", ssl=True)

# Create a GraphQL client
client = Client(transport=transport, fetch_schema_from_transport=True)

# List of genes to query
excel_file = 'GENE-LIST-FILE.xlsx'
rare_prots_df = pd.read_excel(excel_file, engine="openpyxl")
genes = rare_prots_df['Gene_Symbol'].unique()

# Function to generate a GraphQL query for a specific gene
def generate_query(gene_symbol):
    return gql(
        f"""
        query {gene_symbol}Variants {{
          gene(gene_symbol: "{gene_symbol}", reference_genome: GRCh38) {{
            clinvar_variants {{
              variant_id
              clinical_significance
              clinvar_variation_id
            }}

            variants(dataset: gnomad_r4) {{
              variant_id
              chrom
              pos
              rsids
              hgvs
              hgvsc
              hgvsp
              consequence
              
              in_silico_predictors {{
                id
                value
                flags
              }}

              exome {{
                ac
                an
                af
                homozygote_count
              }}

              genome {{
                ac
                an
                af
                homozygote_count
                hemizygote_count
              }}
            }}
          }}
        }}
        """
    )

# Async function to fetch variant data for multiple genes
async def fetch_gnomad_variants():
    for gene in genes:
        if((gene.find('-') == -1)):
            print(gene, end='\t')
            query = generate_query(gene)
            result = await client.execute_async(query)
            
            gene_data = result.get("gene")
            if not gene_data:
                print(f"\nNo data found for {gene}")
                continue
    
            data = {}
            data["Variant ID"] = []
            data["Chromosome"] = []
            data["Position"] = []
            data["RS IDs"] = []
            data["HGVS Consequence"] = []
            data["Coding Change"] = []
            data["Protein Change"] = []
            data["VEP Consequence"] = []
            data["Clinical Significance"] = []
            data["ClinVar Variation ID"] = []
            data["Allele Count (Genome)"] = []
            data["Allele Number (Genome)"] = []
            data["Allele Frequency (Genome)"] = []
            data["Homozygote Count (Genome)"] = []
            data["Hemizygote Count"] = []
            data["Allele Count (Exome)"] = []
            data["Allele Number (Exome)"] = []
            data["Allele Frequency (Exome)"] = []
            data["Homozygote Count (Exome)"] = []
            data["cadd"] = []
            data["spliceai_ds_max"] = []
            data["pangolin_largest_ds"] = []
            data["phylop"] = []
            variants = gene_data["variants"]
            clinvar_variants = gene_data["clinvar_variants"]
    
            # Map ClinVar data by variant_id for easy lookup
            clinvar_dict = {cv["variant_id"]: cv for cv in clinvar_variants}
    
            for variant in variants:
                # Variant ID
                if(variant["variant_id"] != None):
                    data["Variant ID"].append(variant["variant_id"])
                else:
                    data["Variant ID"].append("-")
                # Chromosome
                if(variant["chrom"] != None):
                    data["Chromosome"].append(variant["chrom"])
                else:
                    data["Chromosome"].append("-")
                # Position
                if(variant["pos"] != None):
                    data["Position"].append(variant["pos"])
                else:
                    data["Position"].append("-")
                # RS IDs: An RS ID (Reference SNP cluster ID) is a unique identifier assigned by dbSNP (Database of Single Nucleotide Polymorphisms) to a specific genetic variant
                if(variant["rsids"] != None):
                    if(len(variant["rsids"]) > 0):
                        rsids = ""
                        ii = 0
                        for rsid in variant["rsids"]:
                            ii += 1
                            if(ii == len(variant["rsids"])):
                                rsids = rsids + str(rsid)
                            else:
                                rsids = rsids + str(rsid) + ","
                        data["RS IDs"].append(rsids)
                    else:
                        data["RS IDs"].append("-")
                else:
                    data["RS IDs"].append("-")
                # HGVS Consequence
                if(variant["hgvs"] != None):
                    data["HGVS Consequence"].append(variant["hgvs"])
                else:
                    data["HGVS Consequence"].append("-")    
                # Coding Change
                if(variant["hgvsc"] != None):
                    data["Coding Change"].append(variant["hgvsc"])
                else:
                    data["Coding Change"].append("-")
                # Protein Change
                if(variant["hgvsp"] != None):
                    data["Protein Change"].append(variant["hgvsp"])
                else:
                    data["Protein Change"].append("-")
                # VEP Consequence
                if(variant["consequence"] != None):
                    data["VEP Consequence"].append(variant["consequence"])
                else:
                    data["VEP Consequence"].append("-")
                
                #  Match variant with ClinVar data
                clinvar_data = clinvar_dict.get(variant["variant_id"])
                if(clinvar_data):
                    data["Clinical Significance"].append(clinvar_data.get("clinical_significance", "-"))
                    data["ClinVar Variation ID"].append(clinvar_data.get("clinvar_variation_id", "-"))
                else:
                    data["Clinical Significance"].append("-")
                    data["ClinVar Variation ID"].append("-")
    
                # Population Data (Genome) 
                if variant["genome"]:
                    data["Allele Count (Genome)"].append(variant["genome"].get("ac", "-"))
                    data["Allele Number (Genome)"].append(variant["genome"].get("an", "-"))
                    data["Allele Frequency (Genome)"].append(variant["genome"].get("af", "-"))
                    data["Homozygote Count (Genome)"].append(variant["genome"].get("homozygote_count", "-"))
                    data["Hemizygote Count"].append(variant["genome"].get("hemizygote_count", "-"))
                else:
                    data["Allele Count (Genome)"].append("-")
                    data["Allele Number (Genome)"].append("-")
                    data["Allele Frequency (Genome)"].append("-")
                    data["Homozygote Count (Genome)"].append("-")
                    data["Hemizygote Count"].append("-")
    
                # Population Data (Exome)
                if variant["exome"]:
                    data["Allele Count (Exome)"].append(variant["exome"].get("ac", "-"))
                    data["Allele Number (Exome)"].append(variant["exome"].get("an", "-"))
                    data["Allele Frequency (Exome)"].append(variant["exome"].get("af", "-"))
                    data["Homozygote Count (Exome)"].append(variant["exome"].get("homozygote_count", "-"))
                else:
                    data["Allele Count (Exome)"].append("-")
                    data["Allele Number (Exome)"].append("-")
                    data["Allele Frequency (Exome)"].append("-")
                    data["Homozygote Count (Exome)"].append("-")
    
                is_cadd = False
                is_spliceai_ds_max = False
                is_pangolin_largest_ds = False
                is_phylop = False
                for predictor in variant.get("in_silico_predictors", []):
                    if(('id' in predictor) and ('value' in predictor)):
                        if((predictor['id'] == 'cadd') and (not is_cadd)):
                            data["cadd"].append(predictor['value'])
                            is_cadd = True
                        elif((predictor['id'] == 'spliceai_ds_max') and (not is_spliceai_ds_max)):
                            data["spliceai_ds_max"].append(predictor['value'])
                            is_spliceai_ds_max = True
                        elif((predictor['id'] == 'pangolin_largest_ds') and (not is_pangolin_largest_ds)):
                            data["pangolin_largest_ds"].append(predictor['value'])
                            is_pangolin_largest_ds = True
                        elif((predictor['id'] == 'phylop') and (not is_phylop)):
                            data["phylop"].append(predictor['value'])
                            is_phylop = True
    
                if(not is_cadd):
                    data["cadd"].append("-")
                if(not is_spliceai_ds_max):
                    data["spliceai_ds_max"].append("-")
                if(not is_pangolin_largest_ds):
                    data["pangolin_largest_ds"].append("-")
                if(not is_phylop):
                    data["phylop"].append("-")
            gene_variants_df = pd.DataFrame(data)
            gene_variants_df.to_csv(f'Rare_Diseases_Genes_Variants/{gene}_variants.csv', index=False)
            await asyncio.sleep(60)

# Run the async function properly
asyncio.run(fetch_gnomad_variants())
