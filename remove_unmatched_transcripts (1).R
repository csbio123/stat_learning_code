
#library(glmnet)

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc

#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 

args[1] = "C:/Users/spjtcoi/Google Drive/bmi_train_dataset_nochip.csv"
args[2] = "C:/Users/spjtcoi/Google Drive/bmi_test_dataset_nochip.csv"
args[3] = "/Users/ti1/Google\ Drive/raw\ data/output/"
args[4] = 14
args[5] = 21

print(paste0('GLMNET analaysis'))
print(paste0('Dataset train:', args[1]))
print(paste0('Dataset test: ', args[2] ))
print(paste0('Output dir: ', args[3] ))
print(paste0('gene ex starts (training):', args[4] ))
print(paste0('gene ex starts (test):', args[5] ))

data_train = read.csv(args[1], stringsAsFactors = FALSE)
data_test = read.csv(args[2], stringsAsFactors = FALSE)

same_columns = intersect(colnames(data_train[args[4]:length(data_train)]), colnames(data_test[args[5]:length(data_test)]))



#X_train = data_train[, -which(colnames(data_train) %in% target_columns)]
data_test2 = data_test[, which(colnames(data_test) %in% same_columns)] #isolated genes in separate file
datasets_merged<-cbind(data_test2[,(2:args[4])], data_test2)

previous_na_action <- options('na.action')
options(na.action='na.pass')


output_dir=args[3]
dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)

number = Sys.getenv(c("SGE_TASK_ID"))
write.csv(datasets_merged, file=paste0(output_dir, "/dataset_", number, ".csv"))
