#!/opt/apps/general/python/3.5.1/bin/python3


import pandas as pd
import numpy as np
import subprocess
from Bio import SeqIO
import pathlib
import os 
import subprocess
import sys
import random
random.seed(10)





# In[224]:
#Running Manually
#This is the FASTA file which contains all the non-brain-protein-sequences from UNIPROT
#input_fasta = "/users/spjtcoi/git/stat_learning_code/gwas/query_proteins.fasta"
#This is the 'reference' sequence (e.g. influenza A) that we will use to align our sequences to in FASTA format
#reference = '/users/spjtcoi/git/stat_learning_code/gwas/InfluenzaA_1918_Brevig.fasta  '
#Any output directory
#output_files = "/users/spjtcoi/git/stat_learning_code/gwas/output_Non-CNS/"

#MINIMUM total length of amino acids
#minimal_length = 90000
#MAXIMUM total length of amino acids
#maximal_length = 98000
#Number of subsets to generate
#number_subsets = 10


def create_sequence_dict(input_fasta, length_peptide):
    #Put FASTA sequence into a better datatype
    sequence_dict = {}
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(seq_record.seq)

        sw_sequence= [seq[i:i+length_peptide] for i in range(len(seq)-length_peptide-1)]
        #print(seq_record.id)
        sequence_dict[seq_record.id] = [len(seq), sw_sequence]
    return sequence_dict


def create_sequence_dict_target(input_fasta, length_peptide):
    #Put FASTA sequence into a better datatype
    sequence_dict = {}
    for seq_record in SeqIO.parse(input_fasta, "fasta"):
        seq = str(seq_record.seq)
        sw_sequence= [seq[i:i+length_peptide] for i in range(len(seq)-length_peptide-1)]
        sequence_dict[seq_record.id] = [len(seq), sw_sequence]
    return sequence_dict


def write_control(number_subsets, minimal_length, maximal_length, sequence_dict, output_files):
    #Generate subset, allowing one protein to be present across subsets, but only once within the set
    all_subsets = []
    paths = []
    all_ids = np.array(list(sequence_dict.keys()))

    for i in range(0, number_subsets):
        iterations = 0
        found = False
        while((not found) | (iterations < 10000)):
            iterations += 1
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
                string += ">"+str(id)+ '_'+str(j)+"\n"+seq+"\n"
        with open(paths[num], "w") as text_file:
            text_file.write(string)
    return paths

def write_target(input_fasta, sequence_dict, output_files):
    string = ""
    paths = []
    out_f = output_files + "/" + os.path.basename(input_fasta)
    #Save each query as one output file
    for id in sequence_dict.keys():
        seq = sequence_dict[id]
        length = seq[0]
        for i in range(0, len(seq[1])):
            string += ">"+str(id)+ '_'+str(i)+"\n"+seq[1][i]+"\n"
    with open(out_f, "w") as text_file:
        text_file.write(string)
    paths.append(out_f)
    return paths

def  run_index_building(software, reference, index_name):
    #indexing of reference sequence
    index_cmd = 'java -jar %s -a index -d %s  -i %s -f' % (software, reference, index_name)
    p = subprocess.Popen(index_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
    for line in p.stdout.readlines():
        continue
        #`print (line)
    retval = p.wait()
    return retval


def run_query(software, index_name, paths):
    #Alignments of each subset to the reference sequence. outputing the result to a file
    for i, query_path in enumerate(paths):
        print(query_path)
        output = str(query_path+".out")
        index_cmd = 'java -jar %s -a query -i %s -Q %s -o %s' % (software, index_name, query_path, output)
        print(index_cmd)
        p = subprocess.Popen(index_cmd, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
        for line in p.stdout.readlines():
            continue
            #print (line)
        retval = p.wait()


def run(input_fasta, reference, minimal_length, maximal_length, number_subsets, length_peptide, software, output_files):

    pathlib.Path(output_files).mkdir(parents=True, exist_ok=True)
    index_name = os.path.basename(reference)
    index_name = os.path.splitext(index_name)[0]

    print("Processing input-queries \"{}\" ...".format(input_fasta))
    seq_dict = create_sequence_dict(input_fasta, length_peptide)

    if (not len(seq_dict)):
        print("Query seuence couldn't be read. Please check input format and file name.")
        exit()

    print("Splitting protein into {} peptides...".format(length_peptide))

    if (number_subsets > 0):
        print("Generating randomly {} subsets...".format(number_subsets))
        paths = write_control(number_subsets, minimal_length, maximal_length, seq_dict, output_files)

    else:
        paths = write_target(input_fasta, seq_dict, output_files)

    print(paths)
    if (len(paths) == 0):
        print("No input files have been generated")
        exit()

    print("Building index from reference \"{}\" sequence...".format(reference))
    val = run_index_building(software, reference, index_name)

    print("Running queries...")
    run_query(software, index_name, paths)

    print("All done. Enjoy your day!")

if __name__ == "__main__":

    # Running on cluster
    # This is the FASTA file which contains all the non-brain-protein-sequences from UNIPROT
    input_fasta = "/Users/ti1/Downloads/config/gwas/target_test.fasta" #sys.argv[1]  #
    # This is the 'reference' sequence (e.g. influenza A) that we will use to align our sequences to in FASTA format
    reference = '/Users/ti1/Downloads/config/gwas/influenza.fasta' # sys.argv[2]  #
    # Any output directory
    output_files = "/Users/ti1/Downloads/config/gwas/output_test/" # sys.argv[3]  #

    # MINIMUM total length of amino acids
    minimal_length = 90000 #int(sys.argv[4])
    # MAXIMUM total length of amino acids
    maximal_length = 100000 #int(sys.argv[5])
    # Number of subsets to generate
    number_subsets = 0 # int(sys.argv[6])
    length_peptide = 5# int(sys.argv[7])

    software = "/Users/ti1/Downloads/config/gwas/PeptideMatchCMD_1.0.jar"
    # In[ ]:
    run(input_fasta, reference, minimal_length, maximal_length, number_subsets, length_peptide, software, output_files)



# python3 multi-peptide-match.py query_proteins.fasta InfluenzaA_1918_Brevig.fasta /users/spjtcoi/git/stat_learning_code/gwas/gails/ 90000 100000 11 5
#				query_sequences	influenza_sequence output_directory   amino_acid_length_range   subset_number query_sequence_length