#!/bin/sh
#$ -N CVvaldhba
#$ -q LowMemShortterm.q
#$ -t 1-500
#$ -pe smp 10
#$ -l h_rt=00:20:00
#$ -o /dev/null
#$ -e /dev/null

#LowMemLongterm.q
#HighMemLongterm.q
#HighMemShortterm.q
#LowMemShortterm.q

echo $SGE_TASK_ID

data="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/PREDICTION/validation_input/models_using_chipID_well/hba/dataset_${SGE_TASK_ID}.csv"
feature_list="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/FEATURE_LIST/featurelist_hba_2categ_validation_chip_well.csv"
features_to_remove="EMPTY" 
#if nothing to remove then just add EMPTY (or any word) without brackets eg.features_to_remove="EMPTY" 
output_dir="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/PREDICTION/output/CV_in_validationset_chip_well/hba"


#demographic_only.csv
#clinical_only.csv
#chip_only.csv
#technical_only.csv
#transcripts_only.csv
#unattributed_only.csv

# input parameters via commandline (trouble-shooting mode only)
#output_dir=$1 #comment this out if running as batch
#removing_features=$2 # comment this out if running as batch


mkdir -p $output_dir #-p create multiple subdirectories on the fly (ie they dont already exist)

module load bioinformatics/R/3.4.1
Rscript /users/spjtcoi/git/stat_learning_code/run_dataset_generic_KCL.R $data $feature_list $features_to_remove $output_dir> $output_dir/test.${SGE_TASK_ID}.out 2> $output_dir/test.${SGE_TASK_ID}.err


#Commandline demo
#SGE_TASK_ID=50 ./run_script_amended.sh "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/genes_test" "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/remove.csv"

#Batch submission example
#qsub run_script_amended.sh
