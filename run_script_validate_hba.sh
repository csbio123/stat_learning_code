#!/bin/sh
#$ -N hbapred
#$ -q LowMemShortterm.q
#$ -t 1-500
#$ -pe smp 10
#$ -l h_rt=00:10:00
#$ -o /dev/null
#$ -e /dev/null

#LowMemLongterm.q
#HighMemLongterm.q
#HighMemShortterm.q
#LowMemShortterm.q


echo $SGE_TASK_ID

training_input="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/IMPUTATION_SETS/input_files_HBA_2catg/dataset_${SGE_TASK_ID}.csv"
test_input="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/PREDICTION/validation_input/hba/dataset_${SGE_TASK_ID}.csv"
output_dir="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/PREDICTION/output/HBA/"


# input parameters via commandline (trouble-shooting mode only)
#output_dir=$1 #comment this out if running as batch
#removing_features=$2 # comment this out if running as batch


mkdir -p $output_dir #-p create multiple subdirectories on the fly (ie they dont already exist)

module load bioinformatics/R/3.4.1
Rscript /users/spjtcoi/git/stat_learning_code/classification_prediction_hba.R  $training_input $test_input $output_dir> $output_dir/test.${SGE_TASK_ID}.out 2> $output_dir/test.${SGE_TASK_ID}.err


#Commandline demo
#SGE_TASK_ID=50 ./run_script_amended.sh "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/genes_test" "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/remove.csv"

#Batch submission example
#qsub run_script_amended.sh
