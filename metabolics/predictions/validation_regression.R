
library(glmnet)
library(caret)
library(randomForest)

before <- Sys.time() 

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc
#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 
folder = "/users/spjtcoi/brc_scratch/wgcna/Input/COMBAT"
#folder="/Users/ti1/Downloads/config/validation_prediction"
args[1] = paste0(folder, "/train_combat_9_03_18/training_set_13.csv")
args[2] = paste0(folder, "/validation_combat_9_03_18/validation_set_13.csv")
args[3] = paste0(folder, "/train_combat_9_03_18/training_set_gene_expression.csv")
args[4] = paste0(folder, "/validation_combat_9_03_18/validation_set_gene_expression.csv")
args[5] = "/users/spjtcoi/brc_scratch/wgcna/Input/FEATURES/feature_list_bmi_combat.csv"
args[6] = paste0(folder, "/NULL.txt")
args[7] = paste0(folder, "/test_output-dir_28-04-18")
args[8] = 10#label swapping permutation to establish null prediction rate

iterations = 10 #as.numeric(args[8])
print(paste0('GLMNET analaysis'))
print(paste0('Dataset train:', args[1]))
print(paste0('Dataset test: ', args[2] ))
print(paste0('Genes train: ', args[3] ))
print(paste0('Genes test: ', args[4] ))
print(paste0('All features: ', args[5] ))
print(paste0('To remove features: ', args[6] ))
print(paste0('Output-dir: ', args[7] ))
print(paste0('Permutations: ', iterations ))

number = Sys.getenv(c("SGE_TASK_ID"))


data_train = read.csv(args[1], stringsAsFactors = FALSE)
data_test = read.csv(args[2], stringsAsFactors = FALSE)
gene_expression_train = read.csv(args[3], stringsAsFactors = FALSE, check.names=FALSE)
gene_expression_test = read.csv(args[4], stringsAsFactors = FALSE, check.names=FALSE)
g_test_col = (colnames(gene_expression_test))
g_train_col = (colnames(gene_expression_train))

colnames(gene_expression_test) = unlist(lapply(g_test_col, function(x) { gsub('\`','', x)}))
colnames(gene_expression_train) = unlist(lapply(g_train_col, function(x) { gsub('\`','', x)}))


#This is the fike which contains all features and its corresponding class
#Column 1: Feature name
#Column 2: Feature class
if(!file.exists(args[5])) {
  stop("Feature reference file not found. Please check.")
} else {
  feature_list = read.csv(args[5], stringsAsFactors = FALSE)
}
feature_list[,1] = (unlist(lapply(feature_list[,1], function(x) { gsub('\`','', x)})))

	#This file lists all classes and/or features that should be removed. It is a 1-dim list.
#If this file doesn't exist, nothing will be removed
if(file.exists(args[6])) {
  print(paste0('Features to remove: ', args[6] ))
  features_to_remove = unlist(read.csv(args[6], stringsAsFactors = FALSE, header=FALSE))
} else {
  features_to_remove = c()
}

print(paste0('All number of features in feature list: ', nrow(feature_list)))
print(paste0('All number of fetures in gene expression train: ', ncol(gene_expression_train)))
print(paste0('All number of features in gene expression test: ', ncol(gene_expression_test)))

target_columns = feature_list[feature_list[, 2] == "target",][["feature"]]

data_train = cbind(data_train, gene_expression_train)
data_test = cbind(data_test, gene_expression_test)

data_train = data_train[, which(colnames(data_train) %in% feature_list[,1])]
data_test = data_test[, which(colnames(data_test) %in% feature_list[,1])]



print(paste0('Number of features in train set also in feature list: ', ncol(data_train)))
print(paste0('Number of features in test set also in feature list: ', ncol(data_test)))
same_columns = intersect(colnames(data_train), colnames(data_test))
print("Total number of features that are same between test and train:")
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



print(paste0('All number of features after feature removal', ncol(X_train)))
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


colnames(X_test) = unlist(lapply(colnames(X_test), function(x) { gsub('\`','feature_', x)}))
colnames(X_train) = unlist(lapply(colnames(X_train), function(x) { gsub('\`','feature_', x)}))

same_columns = intersect(colnames(X_train), colnames(X_test))
print(same_columns)
X_train = X_train[, which(colnames(X_train) %in% same_columns)]
X_test = X_test[, which(colnames(X_test) %in% same_columns)]

# Train the model
print("Training Random Forest and GLMNet prediction models")#
d = cbind(Y_train, X_train)
print('Run RF')
test_rf = randomForest(Y_train~., data=d, ntree=100, proximity=T)

print('fit')
preds_rf = predict(test_rf, newdata=X_test)
print('predict')
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
