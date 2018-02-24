import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
import pathlib
import os 
import subprocess
import sys

# In[224]:

#This is the FASTA file which contains all the non-brain-protein-sequences from UNIPROT
input_fasta = sys.argv[1]#"/Users/ti1/Downloads/uniprot-non-brain.fasta"
#This is the 'reference' sequence (e.g. influenza A) that we will use to align our sequences to in FASTA format
reference = sys.argv[2]#'/Users/ti1/Downloads/influenza.fasta'
#Any output directory
output_files = sys.argv[3]#"/Users/ti1/Downloads/query_files/"

#MINIMUM total length of amino acids
minimal_length = int(sys.argv[4])
#MAXIMUM total length of amino acids
maximal_length = int(sys.argv[5])
#Number of subsets to generate
number_subsets = int(sys.argv[6])


# In[ ]:
pathlib.Path(output_files).mkdir(parents=True, exist_ok=True)
index_name = os.path.basename(reference)


# In[225]:
#Put FASTA sequence into a better datatype
sequence_dict = {}
for seq_record in SeqIO.parse(input_fasta, "fasta"):
    seq = str(seq_record.seq)
    gene = 
    sw_sequence= [seq[i:i+5] for i in range(len(seq)-4)]
    sequence_dict[seq_record.id] = [len(seq), sw_sequence]


# In[226]:


#Generate subset, aloowing one protein to be present across subsets, but only once within the set
all_subsets = []
all_ids = np.array(list(sequence_dict.keys()))
for i in range(0, number_subsets-1):
    found = False
    while(not found):
        sample_ids = all_ids[np.random.choice(len(all_ids), 150, replace=False)]
        total_length = sum([ sequence_dict[x][0] for x in sample_ids])
        found = minimal_length < total_length < maximal_length
       	sw_sequences = ([ sequence_dict[x][1] for x in sample_ids])
        ids_p = ([ np.repeat(x, len(sequence_dict[x][1])) for x in sample_ids])

    all_subsets.append([ids_p, sw_sequences])
    
#Save each query as one output file
paths = [ output_files+"query_%s.fasta"%(str(i)) for i,e in enumerate(all_subsets)]
for num, subset in enumerate(all_subsets):
   length = len(subset[0])
   string = ""
   for i in range(0, length):
       ids = subset[0][i]
       seqs = subset[1][i]
       for j in range(0, len(ids)):
          seq = seqs[j]
          id = ids[j]
          string += ">"+id+"\n"+seq+"\n"
   with open(paths[num], "w") as text_file:
     text_file.write(string)
# In[230]:


#indexing of reference sequence
index_cmd = 'java -jar /Users/ti1/Downloads/PeptideMatchCMD_1.0.jar -a index -d %s  -i %s -f' % (reference, index_name)
p = subprocess.Popen(index_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
for line in p.stdout.readlines():
    print (line)
retval = p.wait()


# In[244]:


#Alignments of each subset to the reference sequence. outputing the result to a file
for i, query_path in enumerate(paths):
    print(query_path)
    output = str(query_path+".out")
    index_cmd = 'java -jar /Users/ti1/Downloads/PeptideMatchCMD_1.0.jar -a query -i %s -Q %s -o %s' % (index_name, query_path ,output)
    p = subprocess.Popen(index_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
     print (line)
    retval = p.wait()

