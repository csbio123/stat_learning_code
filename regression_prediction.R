
library(glmnet)

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc

#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 

args[1] = "TRAIN IMPUTATION DATASET"
args[2] = "VALIDATION IMPUTATED DATASET"
args[3] = "/Users/ti1/Google\ Drive/raw\ data/output/"

print(paste0('GLMNET analaysis'))
print(paste0('Dataset train:', args[1]))
print(paste0('Dataset test: ', args[2] ))
print(paste0('Output dir: ', args[3] ))

data_train = read.csv(args[1], stringsAsFactors = FALSE)
data_test = read.csv(args[2], stringsAsFactors = FALSE)

same_columns = intersect(colnames(data_train), colnames(data_test))

data_train = data_train[, which(colnames(data_train) %in% same_columns)]
data_test = data_test[, which(colnames(data_test) %in% same_columns)]

target_columns = "BMI.norm"
Y_train = data_train[, which(colnames(data_train) %in% target_columns)]
Y_test = data_test[, which(colnames(data_test) %in% target_columns)]

X_train = data_train[, -which(colnames(data_train) %in% target_columns)]
X_test = data_test[, -which(colnames(data_test) %in% target_columns)]
previous_na_action <- options('na.action')
options(na.action='na.pass')

design_matrix = model.matrix(Y_train ~ ., data=X_train)
null_rows =  rowSums(is.na(design_matrix)) > 0

X_train = design_matrix[!null_rows, ]
Y_train  = Y_train[!null_rows]
# Do your stuff...

design_matrix_test = model.matrix(~.-Y_test, data=X_test)

# Train the model
glmnet.fit <- cv.glmnet(x=X_train, y=Y_train, family='gaussian', alpha=0)

# Generate predictions
preds <- predict(glmnet.fit, newx=design_matrix_test, type='response', s='lambda.min')

# Put results into dataframe for plotting.
results <- data.frame(pred=preds, actual=Y_test)

output_dir=args[3]
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

number = Sys.getenv(c("SGE_TASK_ID"))
write.table(results, , file=paste0(output_dir, "/predictions_", number))
