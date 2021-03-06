#!/bin/sh
#$ -N purple7_only
#$ -q LowMemShortterm.q
#$ -t 1-500
#$ -pe smp 10
#$ -l h_rt=00:15:00
#$ -o /dev/null
#$ -e /dev/null

#LowMemLongterm.q
#HighMemLongterm.q
#HighMemShortterm.q
#LowMemShortterm.q

echo $SGE_TASK_ID

data="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/IMPUTATION_SETS/input_files_BMI_2catg/dataset_${SGE_TASK_ID}.csv"
feature_list="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/WGCNA/all_features_bmi_2categ_modules.csv"
features_to_remove="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/WGCNA/remove_files/purple7_only.csv" 
#if nothing to remove then just add EMPTY (or any word) without brackets eg.features_to_remove="EMPTY" 
output_dir="/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/WGCNA/categorical_BMI/purple7_only"


#chip
#technical
#clinical
#demographic

#tan
#black
#magenta
#turquoise
#blue
#pink
#purple
#midnightblue
#red
#yellow
#lightcyan
#grey
#salmon
#greenyellow
#green
#cyan
#brown


# input parameters via commandline (trouble-shooting mode only)
#output_dir=$1 #comment this out if running as batch
#removing_features=$2 # comment this out if running as batch


mkdir -p $output_dir #-p create multiple subdirectories on the fly (ie they dont already exist)

module load bioinformatics/R/3.4.1
Rscript /users/spjtcoi/git/stat_learning_code/metabolics/predictions/run_dataset_generic_KCL.R $data $feature_list $features_to_remove $output_dir> $output_dir/test.${SGE_TASK_ID}.out 2> $output_dir/test.${SGE_TASK_ID}.err


#Commandline demo
#SGE_TASK_ID=50 ./run_script_amended.sh "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/genes_test" "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/remove.csv"

#Batch submission example
#qsub run_script_amended.sh
