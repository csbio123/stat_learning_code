import multi_peptide_match
import process
import combine_results
import os



# In[34]:
if __name__ == "__main__":

    # This is the FASTA file which contains all the non-brain-protein-sequences from UNIPROT
    target_fasta = "/users/spjtcoi/brc_scratch/gwas_influenza/target_test.fasta"  # sys.argv[1]  #
    control_fasta = "/users/spjtcoi/brc_scratch/gwas_influenza/non_brain.fasta"  # sys.argv[2]  #
    # This is the 'reference' sequence (e.g. influenza A) that we will use to align our sequences to in FASTA format
    reference = '/users/spjtcoi/brc_scratch/gwas_influenza/influenza.fasta'  # sys.argv[3]  #
    # Any output directory
    output_files = "/users/spjtcoi/brc_scratch/gwas_influenza/pipeline_output/"  # sys.argv[4]  #
    summary_output_files = "/users/spjtcoi/brc_scratch/gwas_influenza/pipeline_output_summary/"  # sys.argv[5]  #can accumulate summary files for all analysis groups in ths directpry
    software = "/users/spjtcoi/brc_scratch/gwas_influenza/PeptideMatchCMD_1.0.jar"#sys.argv[6]

    # MINIMUM total length of amino acids
    minimal_length = 90000  # int(sys.argv[7])
    # MAXIMUM total length of amino acids
    maximal_length = 100000  # int(sys.argv[8])
    # Number of subsets to generate
    number_subsets = 1  # int(sys.argv[9])
    length_peptide = 5  # int(sys.argv[10])

    output_files = os.path.normpath(output_files)
    summary_output_files = os.path.normpath(summary_output_files)

    print("Running multi peptide match on control...")
    multi_peptide_match.run(control_fasta, reference, minimal_length, maximal_length, number_subsets, length_peptide, software, output_files+"/")

    print("Running multi peptide match on target...")
    multi_peptide_match.run(target_fasta, reference, minimal_length, maximal_length, 0, length_peptide, software, output_files+"/")


    print("Processing outputs...")
    process.run(output_files)

    print("Combining results...")
    target_csv = os.path.basename(target_fasta) + "_result.csv"
    combine_results.run(output_files, output_files+"/"+target_csv, summary_output_files)


