
#####USE AS FOLLOWS#####:
##Rscript p_value_analysis.R /users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/PREDICTION/output/matched_GXdata/leave_one_in/bmi####

args <-commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc
print(paste0('Input directory for p value analysis:', args[1]))#This line will tell you the input directory

path = "/Users/ti1/Google\ Drive/non_normalised/"


library(caret)
filenames <-
  list.files(path,
             pattern = "*predictions*",
             full.names = TRUE,
             recursive = TRUE)

files = sapply(
  filenames, function(x){ 
  d = read.csv(x, sep = " ",  header = TRUE)
  name = basename(dirname(x))
  d = (cbind(name, d))
  return (d)
  }, simplify=FALSE)

files_df = do.call(rbind, files)

hist(files_df[["alpha"]])

accuracies = sapply( files, function(x) {
    confusionMatrix(as.factor(x[,4]),as.factor(x[,5]))$byClass
  }, simplify = FALSE)
accuracies = do.call(rbind, accuracies)


meta = sapply(
  files, function(x){ 
    name = as.character((x[,1][1]))
    id = x[,2][1]
    return (c(name, id))
  }, simplify=FALSE)

meta = do.call(rbind, meta)
rownames(meta) = c()

run_type = as.data.frame(meta[,1])
run_id = as.data.frame(as.numeric(meta[,2]))

rownames(accuracies) = c()
accuracies = cbind(run_type, run_id, accuracies)
cols = colnames(accuracies)
cols[1] = "run_type"
cols[2] = "run_id"
colnames(accuracies) = gsub( " ", "_", cols)


my_plot = ggplot(data = accuracies, aes(run_type, Sensitivity), outlier.shape = NA) + geom_boxplot(aes(color =
                                                                                                         run_type)) 


my_plot = ggplot(data = accuracies, aes(run_type, Specificity), outlier.shape = NA) + geom_boxplot(aes(color =
                                                                                                         run_type)) 
