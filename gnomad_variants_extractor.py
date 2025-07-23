'''
Author: Sina Abdollahi
Release Date: 23.07.2025
Title: Extracting the variants information from gnomAD 

Description: Given a list of gene symbols in an Excel file with header of "Gene_Symbol" and we extract 
the variants' information such as Expected Loss of function, Observed Loss of function, 
Observed/Expected Loss of function, Probability of being loss-of-function (pLI), 
Loss of function Z score, Observed and Expected missense and synonymous variants and their Z scores.
'''

import asyncio
from gql import gql, Client
from gql.transport.aiohttp import AIOHTTPTransport
import pandas as pd
import time

#Put your main Excel file containing Gene Symbols next to this Python file
excel_file = 'GENE-LIST-FILE.xlsx'
all_genes_df = pd.read_excel(excel_file, engine="openpyxl")

# Define the gnomAD API endpoint
transport = AIOHTTPTransport(url="https://gnomad.broadinstitute.org/api", ssl=True)

# Create a GraphQL client
client = Client(transport=transport, fetch_schema_from_transport=True)

mutinfo = {}
mutinfo['Gene_Symbol'] = []
mutinfo['Expected_LoF_SNVs'] = []
mutinfo['Observed_LoF_SNVs'] = []
mutinfo['O/E_LoF'] = []
mutinfo['pLI'] = []
mutinfo['LoF_Z_Score'] = []
mutinfo['Expected_Missense_SNVs'] = []
mutinfo['Observed_Missense_SNVs'] = []
mutinfo['O/E_Missense'] = []
mutinfo['Missense_Z_Score'] = []
mutinfo['Expected_Synonymous_SNVs'] = []
mutinfo['Observed_Synonymous_SNVs'] = []
mutinfo['O/E_Synonymous'] = []
mutinfo['Synonymous_Z_Score'] = []

gene_names = all_genes_df['Gene_Symbol'].unique()
i = 0
for gene_name in gene_names:
    query = gql("""
    query GeneConstraint($gene: String!) {
      gene(gene_symbol: $gene, reference_genome: GRCh38) {
        gnomad_constraint {
          lof_z
          mis_z  # Corrected field name for missense
          syn_z  # Corrected field name for synonymous
          pli
          exp_lof    # Expected LoF SNVs
          obs_lof    # Observed LoF SNVs
          exp_mis    # Expected Missense SNVs
          obs_mis    # Observed Missense SNVs
          exp_syn    # Expected Synonymous SNVs
          obs_syn    # Observed Synonymous SNVs
        }
      }
    }
    """)

    # Async function to fetch data
    async def fetch_gnomad_data(gene_name):
        return await client.execute_async(query, variable_values={"gene": gene_name})

    try:
        result = await fetch_gnomad_data(gene_name)
        exp_lof = result["gene"]["gnomad_constraint"]["exp_lof"]  
        obs_lof = result["gene"]["gnomad_constraint"]["obs_lof"]  
        pli = result["gene"]["gnomad_constraint"]["pli"]
        lof_z = result["gene"]["gnomad_constraint"]["lof_z"]

        mis_z = result["gene"]["gnomad_constraint"]["mis_z"]
        exp_mis = result["gene"]["gnomad_constraint"]["exp_mis"]  
        obs_mis = result["gene"]["gnomad_constraint"]["obs_mis"]

        syn_z = result["gene"]["gnomad_constraint"]["syn_z"]
        exp_syn = result["gene"]["gnomad_constraint"]["exp_syn"]  
        obs_syn = result["gene"]["gnomad_constraint"]["obs_syn"]

        mutinfo['Gene_Symbol'].append(gene_name)

        if(exp_lof is None):
            mutinfo['Expected_LoF_SNVs'].append('-')
        else:
            mutinfo['Expected_LoF_SNVs'].append(exp_lof)

        if(obs_lof is None):
            mutinfo['Observed_LoF_SNVs'].append('-')
        else:
            mutinfo['Observed_LoF_SNVs'].append(obs_lof)

        if((exp_lof == 0) or (obs_lof is None) or (exp_lof is None)):
            mutinfo['O/E_LoF'].append('-')
        else:
            mutinfo['O/E_LoF'].append(obs_lof/exp_lof)

        if(pli is None):
            mutinfo['pLI'].append('-')
        else:
            mutinfo['pLI'].append(pli)

        if(lof_z is None):
            mutinfo['LoF_Z_Score'].append('-')
        else:
            mutinfo['LoF_Z_Score'].append(lof_z)

        if(exp_mis is None):
            mutinfo['Expected_Missense_SNVs'].append('-')
        else:
            mutinfo['Expected_Missense_SNVs'].append(exp_mis)

        if(obs_mis is None):
            mutinfo['Observed_Missense_SNVs'].append('-')
        else:
            mutinfo['Observed_Missense_SNVs'].append(obs_mis)

        if((exp_mis == 0) or (exp_mis is None) or (obs_mis is None)):
            mutinfo['O/E_Missense'].append('-')
        else:
            mutinfo['O/E_Missense'].append(obs_mis/exp_mis)

        if(mis_z is None):
            mutinfo['Missense_Z_Score'].append('-')
        else:
            mutinfo['Missense_Z_Score'].append(mis_z)

        if(exp_syn is None):
            mutinfo['Expected_Synonymous_SNVs'].append('-')
        else:
            mutinfo['Expected_Synonymous_SNVs'].append(exp_syn)

        if(obs_syn is None):
            mutinfo['Observed_Synonymous_SNVs'].append('-')
        else:
            mutinfo['Observed_Synonymous_SNVs'].append(obs_syn)

        if((exp_syn == 0) or (exp_syn is None) or (obs_syn is None)):
            mutinfo['O/E_Synonymous'].append('-')
        else:
            mutinfo['O/E_Synonymous'].append(obs_syn/exp_syn)

        if(syn_z is None):
            mutinfo['Synonymous_Z_Score'].append('-')
        else:
            mutinfo['Synonymous_Z_Score'].append(syn_z)

        i += 1
        # Due to reduce the load over the GnomAD server, the program needs to be suspended for 1 minutes after extracting 20 genes' information
        if((i % 20) == 0):
            mut_df = pd.DataFrame(mutinfo)
            mut_df = pd.concat([already_saved_df, mut_df], axis=0, ignore_index=True)
            mut_df.to_csv('Rare_diseases_prots_gnomAD_scores.csv',index=False)
            print(i, end=', ')
            time.sleep(60)
    except Exception as e:
        print(f"FAILURE: '{gene_name}' does not exist!", end=' ')
