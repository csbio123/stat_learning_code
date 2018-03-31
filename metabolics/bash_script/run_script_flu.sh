#!/bin/sh
#$ -N pep_match
#$ -q LowMemShortterm.q
#$ -t 1
#$ -pe smp 10
#$ -l h_rt=00:10:00
#$ -o /dev/null
#$ -e /dev/null

#LowMemLongterm.q
#HighMemLongterm.q
#HighMemShortterm.q
#LowMemShortterm.q


echo $SGE_TASK_ID

input_dir=/users/spjtcoi/brc_scratch/gwas_influenza
target_fasta=${input_dir}/target_test.fasta 
control_fasta=${input_dir}/non_brain.fasta
reference=${input_dir}/influenza.fasta
output_files=${input_dir}/pipeline_output/
summary_output_files=${input_dir}/pipeline_output_summary/ #can accumulate summary files for all analysis groups in ths directpry
software=${input_dir}/PeptideMatchCMD_1.0.jar
# MINIMUM total length of amino acids
minimal_length=92000 
# MAXIMUM total length of amino acids
maximal_length=97000
# Number of subsets to generate
number_subsets=10
length_peptide=5


echo "input_dir = $input_dir"
echo "target_fasta = $target_fasta"
echo "control_fasta = $control_fasta"
echo "reference= $reference"
echo "output_files = $output_files"
echo "summary_output_files = $summary_output_files"
echo "software = $software"
echo "minimal_length = $minimal_length"
echo "maximal_length = $maximal_length"
echo "number_subsets = $number_subsets"
echo "length_peptide = $length_peptide"



module load general/python/3.5.1


mkdir -p $output_files
mkdir -p $summary_output_files

python3 /users/spjtcoi/brc_scratch/gwas_influenza/all_in_one.py $target_fasta $control_fasta $reference $output_files $summary_output_files $software $minimal_length $maximal_length $number_subsetsn $length_peptide $output_dir> $output_dir/test.${SGE_TASK_ID}.out 2> $output_dir/test.${SGE_TASK_ID}.err $genes_named


