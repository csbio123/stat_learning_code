#!/bin/sh
#$ -N bmipred
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

input_dir=ASDAWDA
training_input=${input_dir}/training_set_${SGE_TASK_ID}.csv
test_input=${input_dir}/validation_set_${SGE_TASK_ID}.csv
training_GE=${input_dir}/training_set_gene_expression.csv
test_GE=${input_dir}/validation_set_gene_expression.csv
feature_list=${input_dir}/validation_features.csv
remove_feature_list=${input_dir}/demo.txt
output_dir=my_output_dir

echo "input_dir = $input_dir"
echo "training_input  = $training_input"
echo "test_input = $test_input"
echo "training_GE = $training_GE"
echo "test_GE = $test_GE"
echo "feature_list = $feature_list"
echo "remove_feature_list = $remove_feature_list"
echo "output_dir = $output_dir"

module load bioinformatics/R/3.4.1
Rscript /users/spjtcoi/git/stat_learning_code/classification_prediction.R $training_input $test_input $training_GE $test_GE $feature_list $remove_feature_list  $output_dir> $output_dir/test.${SGE_TASK_ID}.out 2> $output_dir/test.${SGE_TASK_ID}.err


#Commandline demo
#SGE_TASK_ID=50 ./run_script_amended.sh "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/genes_test" "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/remove.csv"

#Batch submission example
#qsub run_script_amended.sh
