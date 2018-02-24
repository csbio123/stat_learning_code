
# coding: utf-8

# In[32]:


import pandas as pd 
import codecs
import numpy as np
from Bio import SeqIO
from os import listdir
from os.path import isfile, join
import sys


# In[33]:


def join_string(strings):
    return ",".join(strings)

def generate_results(output, sequence_dict, output_file):

    peta_peptides = pd.DataFrame.from_dict(sequence_dict, orient="index")
    peta_peptides.columns = ["peta-peptide"]

    merged_output = pd.merge(output, peta_peptides, left_on = ["protein_pos"], right_index=True)

    ids = [v.split("_")[0] for v in merged_output["protein_pos"].values]
    aa_pos = [v.split("_")[2] for v in merged_output["protein_pos"].values]
    protein = [v.split("|")[1] for v in ids]
    gene = [v.split("|")[2] for v in ids]

    merged_output["protein"] = protein
    merged_output["gene"] = gene
    merged_output["number_x_peptides"] = aa_pos

    protein_AA_counts = merged_output.groupby("protein").count().reset_index()
    total = pd.Series((protein_AA_counts.sum()))
    total.protein = "total"
    protein_AA_counts = protein_AA_counts.append(total, ignore_index=True)

    merged_output = merged_output.dropna(subset=["from"], axis=0)

    peta_peptide_matches = pd.pivot_table(merged_output, values='peta-peptide', index=['protein'], columns=['match'], aggfunc=join_string).fillna(0)
    peta_peptide_matches_number = pd.pivot_table(merged_output, values='peta-peptide', index=['protein'], columns=['match'], aggfunc=len).fillna(0)
    total_matches = peta_peptide_matches_number.sum(0)
    total_matches.name = "total"
    peta_peptide_matches = peta_peptide_matches.append(total_matches)
    peta_peptide_matches = peta_peptide_matches.reset_index()

    peta_peptide_matches = pd.merge(peta_peptide_matches, protein_AA_counts[["number_x_peptides", "protein"]], on="protein")

    final_output = pd.merge(peta_peptide_matches, merged_output[["protein", "gene"]].drop_duplicates(), on="protein", how="left")
    final_output = final_output.set_index(["number_x_peptides", "protein", "gene"])
    print("output_file:" + output_file + "_result.csv")
    final_output.to_csv(output_file+"_result.csv")
    


# In[34]:


mypath = sys.argv[1]
onlyfiles = [f for f in listdir(mypath) if isfile(join(mypath, f))]


# In[35]:


for input_fasta in onlyfiles:
    if input_fasta.endswith(".fasta"):
        input_fasta = mypath +"/" + input_fasta
        print("Processing:" + input_fasta)
        input_r = pd.read_csv(input_fasta+".out", sep='\t')
        input_r.columns = ["protein_pos", "match", "length", "from", "to", "ignore"]
        sequence_dict = {}
        for seq_record in SeqIO.parse(input_fasta, "fasta"):
            seq = str(seq_record.seq)
            sequence_dict[seq_record.id] = str(seq)
        generate_results(input_r, sequence_dict, input_fasta)

