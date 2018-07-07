
library(glmnet)
library(caret)
library(randomForest)

before <- Sys.time() 

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc
#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 
#folder = "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/downloaded/"
#folder="/Users/ti1/Downloads/config/validation_prediction"
#args[1] = paste0(folder, "/train_sva_9_03_18/training_set_13.csv")
#args[2] = paste0(folder, "/validation_sva_9_03_18/validation_set_13.csv")
#args[3] = paste0(folder, "/train_sva_9_03_18/training_set_gene_expression.csv")
#args[4] = paste0(folder, "/validation_sva_9_03_18/validation_set_gene_expression.csv")
#args[5] = paste0(folder, "/feature_list_10_03_18.csv")
#args[6] = paste0(folder, "/remove_class.txt")
#args[7] = paste0(folder, "/test_output-dir")
#args[8] = 10#label swapping permutation to establish null prediction rate

#args[1] = "/users/spjtcoi/brc_scratch/wgcna/Input/SVA/train_sva_9_03_18/training_set_59.csv" 
#args[2] = "/users/spjtcoi/brc_scratch/wgcna/Input/SVA/validation_sva_9_03_18/validation_set_59.csv"
#args[3] = "/users/spjtcoi/brc_scratch/wgcna/Input/SVA/train_sva_9_03_18/training_set_gene_expression.csv"
#args[4] = "/users/spjtcoi/brc_scratch/wgcna/Input/SVA/validation_sva_9_03_18/validation_set_gene_expression.csv"
#args[5] = "/users/spjtcoi/brc_scratch/wgcna/Input/FEATURES/feature_list_bmi_sva.csv"
#args[6] = "/users/spjtcoi/brc_scratch/wgcna/Input/REMOVE/SVA_rm_features/3plus_1color/3plus_brown.txt"
#args[7] = "/users/spjtcoi/brc_scratch/wgcna/Output/SVA/July_Analysis/3plus_brown/"

iterations = 1 #as.numeric(args[8])
print(paste0('GLMNET analaysis'))
print(paste0('Dataset train:', args[1]))
print(paste0('Dataset test: ', args[2] ))
print(paste0('Genes train: ', args[3] ))
print(paste0('Genes test: ', args[4] ))
print(paste0('All features: ', args[5] ))
print(paste0('To remove features: ', args[6] ))
print(paste0('Output-dir: ', args[7] ))
#print(paste0('Permutations: ', iterations ))

number = Sys.getenv(c("SGE_TASK_ID"))


data_train = read.csv(args[1], stringsAsFactors = FALSE)
data_test = read.csv(args[2], stringsAsFactors = FALSE)
gene_expression_train = read.csv(args[3], stringsAsFactors = FALSE, check.names = FALSE)
gene_expression_test = read.csv(args[4], stringsAsFactors = FALSE, check.names = FALSE)

#This is the fike which contains all features and its corresponding class
#Column 1: Feature name
#Column 2: Feature class
if(!file.exists(args[5])) {
  stop("Feature reference file not found. Please check.")
} else {
  feature_list = read.csv(args[5], stringsAsFactors = FALSE,  check.names = FALSE)
}

#This file lists all classes and/or features that should be removed. It is a 1-dim list.
#If this file doesn't exist, nothing will be removed
if(file.exists(args[6])) {
  print(paste0('Features to remove: ', args[6] ))
  features_to_remove = unlist(read.csv(args[6], stringsAsFactors = FALSE, header=FALSE))
} else {
  features_to_remove = c()
}

target_columns = feature_list[feature_list[, 2] == "target",][["feature"]]

data_train = cbind(data_train, gene_expression_train)
data_test = cbind(data_test, gene_expression_test)

data_train = data_train[, which(colnames(data_train) %in% feature_list[,1])]
data_test = data_test[, which(colnames(data_test) %in% feature_list[,1])]

same_columns = intersect(colnames(data_train), colnames(data_test))
print("Total number of columns:")
print(length(same_columns))

#write.csv(same_columns, "~/Downloads/test.csv")
data_train = data_train[, which(colnames(data_train) %in% same_columns)]
data_test = data_test[, which(colnames(data_test) %in% same_columns)]


#write.table(colnames(data_train), file="/Users/ti1/Downloads/feature_list_10_03_18.csv",  quote = FALSE, row.names = FALSE)
print(paste0('Target equals: ', target_columns))#

Y_train = (data_train[, which(colnames(data_train) %in% target_columns)])
Y_test = (data_test[, which(colnames(data_test) %in% target_columns)])


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

X_train = design_matrix[!null_rows, -1]
Y_train  = Y_train[!null_rows]

# Do your stuff...
design_matrix_test = model.matrix(Y_test~., data=X_test)
null_rows =  rowSums(is.na(design_matrix_test)) > 0
X_test = design_matrix_test[!null_rows, -1]
Y_test  = Y_test[!null_rows]


same_columns = intersect(colnames(X_train), colnames(X_test))

X_train = X_train[, which(colnames(X_train) %in% same_columns)]
X_test = X_test[, which(colnames(X_test) %in% same_columns)]

#X_train[["gender"]] = as.factor( X_train[["gender"]] )
#X_test[["gender"]] = as.factor( X_test[["gender"]] )



# Train the model
print("Training Random Forest and GLMNet prediction models")#
d = cbind(Y_train, X_train)

test_rf = randomForest(Y_train~., data=d, ntree=100, proximity=T)
preds_rf = predict(test_rf, newdata=X_test)

sample_results = sapply(1:iterations, function(x) {
  
  #random_features = sample( colnames(data_train), ncol(X_train))
  X_train_r = X_train
  X_test_r = X_test
  
  #X_train_r = data_train[, which(colnames(data_train) %in% random_features)]
  #X_test_r = data_test[, which(colnames(data_test) %in% random_features)]
  Y_rand = data.frame(sample(Y_train))
  colnames(Y_rand) = "Y_rand"
  d_r = cbind(Y_rand, X_train_r)
  test_rf = randomForest(Y_rand~., data=d_r, ntree=100, proximity=T)
  preds_rf = predict(test_rf, newdata=X_test_r)
  return(mean(100*(abs((preds_rf - Y_test))/Y_test)))
})

varImpPlot(test_rf,  
        sort = T,
          n.var=20,
         main="Top 10 - Variable Importance")

var.imp = data.frame(importance(test_rf, type=2))
var.imp$Variables = row.names(var.imp)  
top = var.imp[order(var.imp$IncNodePurity,decreasing = T),]


alphas <- seq(from = 0.0, to = 1, by = 0.1)

alpha_errors = unlist(foreach(i = 1:length(alphas), .inorder = FALSE,.packages = "glmnet", .multicombine = TRUE) %dopar% {
  #Run GLMNET witt our features and our target variable
  output = cv.glmnet(
    x = X_train,
    y = Y_train,
    family = 'gaussian',
    type.measure = "mse", nfolds = 10, alpha=alphas[i])
  
  error_mean =  mean(output$cvm)
  return(error_mean)
})

after<-Sys.time() #stop timing
print(after-before)

best_alpha = alphas[which.max(alpha_errors)]
glmnet.fit <- cv.glmnet(x=X_train, y=Y_train, family='gaussian', alpha = best_alpha ) #0.009

# Generate predictions
print("generating predictions (test sample) based on training data")#
preds_glmnet <- as.numeric(unlist(predict(glmnet.fit, newx = X_test, s='lambda.min')))

# Put results into dataframe for plotting.
print("Compiling results table")
results <- data.frame(run_id=number,prediction_glmnet=(preds_glmnet), prediction_rf=preds_rf, actual=Y_test, alpha=best_alpha)
results["prediction_glmnet_ape"] = 100*(abs((results["prediction_glmnet"] - results["actual"]))/results["actual"])
results["prediction_rf_ape"] = 100*(abs((results["prediction_rf"] - results["actual"]))/results["actual"])
mean(results[["prediction_glmnet_ape"]])
mean(results[["prediction_rf_ape"]])

print("Writing file") #
output_dir=args[7]
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)
write.table(results, file=paste0(output_dir, "/predictions_", number))
rownames(top) = NULL
write.table(top, file=paste0(output_dir, "/variable_importance_", number), row.names = FALSE)
write.table(sample_results, file=paste0(output_dir, "/permut_results_", number), row.names = FALSE)
