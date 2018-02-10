
library(glmnet)
library(caret)

before <- Sys.time() 

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc

#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 

#args[1] = "/Users/ti1/Google\ Drive/raw\ data/training_data/dataset_1.csv"
#args[2] = "/Users/ti1/Google\ Drive/raw\ data/validation_data_classification/dataset_1.csv"
#args[3] = "/Users/ti1/Google\ Drive/raw\ data/output/"

#args[1] = "/Users/ti1/Google\ Drive/raw data/batch_normalized_data/training_set_1.csv"
#args[2] = "/Users/ti1/Google\ Drive/raw data/batch_normalized_data/validation_set_1.csv"
#args[3] = "/Users/ti1/Google\ Drive/raw data/batch_normalized_data/training_set_gene_expression.csv"
#args[4] = "/Users/ti1/Google\ Drive/raw data/batch_normalized_data/validation_set_gene_expression.csv"
#args[5] = "~/Downloads/test.csv"
#args[6] = "/Users/ti1/Google\ Drive/raw\ data/remove_demo.txt"
#args[7] = "/Users/ti1/Google\ Drive/raw\ data/test_output-dir"

print(paste0('GLMNET analaysis'))
print(paste0('Dataset train:', args[1]))
print(paste0('Dataset test: ', args[2] ))
print(paste0('Genes train: ', args[3] ))
print(paste0('Genes test: ', args[4] ))
print(paste0('All features: ', args[5] ))
print(paste0('To remove features: ', args[6] ))
print(paste0('Output-dir: ', args[7] ))

number = Sys.getenv(c("SGE_TASK_ID"))


data_train = read.csv(args[1], stringsAsFactors = FALSE)
data_test = read.csv(args[2], stringsAsFactors = FALSE)
gene_expression_train = read.csv(args[3], stringsAsFactors = FALSE)
gene_expression_test = read.csv(args[4], stringsAsFactors = FALSE)

#This is the fike which contains all features and its corresponding class
#Column 1: Feature name
#Column 2: Feature class
if(!file.exists(args[5])) {
  stop("Feature reference file not found. Please check.")
} else {
  feature_list = read.csv(args[5], stringsAsFactors = FALSE)
}

#This file lists all classes and/or features that should be removed. It is a 1-dim list.
#If this file doesn't exist, nothing will be removed
if(file.exists(args[6])) {
  print(paste0('Features to remove: ', args[6] ))
  features_to_remove = unlist(read.csv(args[6], stringsAsFactors = FALSE, header=FALSE))
} else {
  features_to_remove = c()
}

data_train = cbind(data_train, gene_expression_train)
data_test = cbind(data_test, gene_expression_test)

same_columns = intersect(colnames(data_train), colnames(data_test))

#write.csv(same_columns, "~/Downloads/test.csv")
data_train = data_train[, which(colnames(data_train) %in% same_columns)]
data_test = data_test[, which(colnames(data_test) %in% same_columns)]

target_columns = "BMI.catg" #Choose BMI.catg or hba1c.catg
print(paste0('Target equals: ', target_columns))#

Y_train = factor(data_train[, which(colnames(data_train) %in% target_columns)])
Y_test = factor(data_test[, which(colnames(data_test) %in% target_columns)])

X_train = data_train[, -which(colnames(data_train) %in% target_columns)]
X_test = data_test[, -which(colnames(data_test) %in% target_columns)]

#Remove features from the class
if(length(features_to_remove) >0) {
  
  features_to_remove=unlist(features_to_remove)
  print("Following features will be removed:")
  
  top_level_to_remove = feature_list[feature_list[,2] %in% features_to_remove, ]
  specific_to_remove = feature_list[feature_list[,1] %in% features_to_remove,]
  
  print("Classess:")
  print(unique(top_level_to_remove[,2]))
  
  print("Specific features:")
  print(unique(specific_to_remove[,1]))
  
  remove = unique(c(top_level_to_remove[,1], specific_to_remove[,1]))

  X_train = X_train[,-which(colnames(X_train) %in% remove)]
  X_test = X_test[,-which(colnames(X_test) %in% remove)]
}

previous_na_action <- options('na.action')
options(na.action='na.pass')

design_matrix = model.matrix(Y_train ~ ., data=X_train)
null_rows =  rowSums(is.na(design_matrix)) > 0

X_train = design_matrix[!null_rows, ]
Y_train  = Y_train[!null_rows]

# Do your stuff...
design_matrix_test = model.matrix(~.-Y_test, data=X_test)
null_rows =  rowSums(is.na(design_matrix_test)) > 0
X_test = design_matrix_test[!null_rows, ]
Y_test  = Y_test[!null_rows]

same_columns = intersect(colnames(X_train), colnames(X_test))

X_train = X_train[, which(colnames(X_train) %in% same_columns)]
X_test = X_test[, which(colnames(X_test) %in% same_columns)]

# Train the model
print("Training GLMNet model")#


alphas <- seq(from = 0.1, to = 1, by = 0.001)

alpha_errors = unlist(foreach(i = 1:length(alphas), .inorder = FALSE,.packages = "glmnet", .multicombine = TRUE,.export= c("target","alphas","features")) %dopar% {
  #Run GLMNET witt our features and our target variable
  output = cv.glmnet(
    x = X_train,
    y = Y_train,
    family = 'binomial',
    type.measure = "auc", nfolds = 4)
  
  error_mean =  mean(output$cvm)
  return(error_mean)
})

after<-Sys.time() #stop timing
print(after-before)

best_alpha = alphas[which.max(alpha_errors)]
glmnet.fit <- cv.glmnet(x=X_train, y=Y_train, family='binomial', alpha = best_alpha ) #0.009

# Generate predictions
print("generating predictions (test sample) based on training data")#
preds <- as.numeric(unlist(predict(glmnet.fit, newx = X_test, type='class', s='lambda.min')))

# Put results into dataframe for plotting.
print("Compiling results table")
results <- data.frame(run_id=number,prediction=as.factor(preds), actual=as.factor(Y_test), alpha=best_alpha)

print("Writing file") #
output_dir=args[7]
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)
write.table(results, file=paste0(output_dir, "/predictions_", number))
