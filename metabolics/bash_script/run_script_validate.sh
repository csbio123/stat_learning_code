#!/bin/sh
#$ -N clin_only
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

input_dir=/users/spjtcoi/git/stat_learning_code/metabolics/validation_prediction_new
training_input=${input_dir}/batch_normalized_data_regression/training_set_${SGE_TASK_ID}.csv
test_input=${input_dir}/batch_normalized_data_regression/validation_set_${SGE_TASK_ID}.csv
training_GE=${input_dir}/batch_normalized_data/training_set_gene_expression.csv
test_GE=${input_dir}/batch_normalized_data/validation_set_gene_expression.csv
feature_list=${input_dir}/validation_features_symbols.csv
remove_feature_list=${input_dir}/clinical_only.txt #If nothing to remove then add random (null) filename which will be ignored
output_dir=/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/predict-validate_regression/LEAVE_ONE_IN/clinical_only
genes_named=${input_dir}/gene_names.csv

echo "input_dir = $input_dir"
echo "training_input  = $training_input"
echo "test_input = $test_input"
echo "training_GE = $training_GE"
echo "test_GE = $test_GE"
echo "feature_list = $feature_list"
echo "remove_feature_list = $remove_feature_list"
echo "output_dir = $output_dir"
echo "gene_names = $genes_named"


module load bioinformatics/R/3.4.1
#export SGE_TASK_ID=50

mkdir -p $output_dir
Rscript /users/spjtcoi/git/stat_learning_code/metabolics/predictions/validation_classification_prediction_tomi_regression.R $training_input $test_input $training_GE $test_GE $feature_list $remove_feature_list  $output_dir> $output_dir/test.${SGE_TASK_ID}.out 2> $output_dir/test.${SGE_TASK_ID}.err $genes_named

#Commandline demo
#SGE_TASK_ID=50 ./run_script_amended.sh "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/genes_test" "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/remove.csv"


