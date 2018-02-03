
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc

#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 

#args[1] = "/Users/ti1/Google\ Drive/raw\ data/training_data/dataset_1.csv"
#args[2] = "/Users/ti1/Google\ Drive/raw\ data/validation_data_classification/dataset_1.csv"
#args[3] = "/Users/ti1/Google\ Drive/raw\ data/output/"

args[1] = "/Users/ti1/Google\ Drive/raw data/batch_normalized_data/training_set_1.csv"
args[2] = "/Users/ti1/Google\ Drive/raw data/batch_normalized_data/validation_set_1.csv"
args[3] = "/Users/ti1/Google\ Drive/raw\ data/output/"

print(paste0('GLMNET analaysis'))
print(paste0('Dataset train:', args[1]))
print(paste0('Dataset test: ', args[2] ))
print(paste0('Output dir: ', args[3] ))

data_train = read.csv(args[1], stringsAsFactors = FALSE)
data_test = read.csv(args[2], stringsAsFactors = FALSE)

gene_expression_train = read.csv("/Users/ti1/Google\ Drive/raw data/batch_normalized_data/training_set_gene_expression.csv", stringsAsFactors = FALSE)
gene_expression_test = read.csv("/Users/ti1/Google\ Drive/raw data/batch_normalized_data/validation_set_gene_expression.csv", stringsAsFactors = FALSE)


data_train = cbind(data_train, gene_expression_train)
data_test = cbind(data_test, gene_expression_test)

same_columns = intersect(colnames(data_train), colnames(data_test))

data_train = data_train[, which(colnames(data_train) %in% same_columns)]
data_test = data_test[, which(colnames(data_test) %in% same_columns)]

target_columns = "BMI.catg" #Choose BMI.catg or hba1c.catg
print(paste0('Target equals: ', target_columns))#

Y_train = factor(data_train[, which(colnames(data_train) %in% target_columns)])
Y_test = factor(data_test[, which(colnames(data_test) %in% target_columns)])

X_train = data_train[, -which(colnames(data_train) %in% target_columns)]
X_test = data_test[, -which(colnames(data_test) %in% target_columns)]

X_train=X_train[,1:7]
X_test=X_test[,1:7]

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


alphas <- c(seq(from = 0, to = 0.01, by = 0.0001), seq(from = 0.1, to = 1, by = 0.1))

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

best_alpha = alphas[which.max(alpha_errors)]
print(best_alpha)
glmnet.fit <- cv.glmnet(x=X_train, y=Y_train, family='binomial', alpha = best_alpha ) #0.009

# Generate predictions
print("generating predictions (test sample) based on training data")#
preds <- as.numeric(unlist(predict(glmnet.fit, newx = X_test, type='class', s='lambda.min')))

# Put results into dataframe for plotting.
print("Compiling results table")
results <- data.frame(prediction=preds, actual=Y_test)
colnames(results) = c("pred", "actual")

accuracy = Accuracy(results[["pred"]], results[["actual"]])
recall =  Recall(results[["pred"]], results[["actual"]])
precision =  Precision(results[["pred"]], results[["actual"]])
print(accuracy)
print(recall)
print(precision)

output_dir=arsgs[3]
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

number = Sys.getenv(c("SGE_TASK_ID"))

print("Writing file") #
write.table(results, , file=paste0(output_dir, "/predictions_", number))
