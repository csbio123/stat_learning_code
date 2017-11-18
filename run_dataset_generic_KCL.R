#source("http://bioconductor.org/biocLite.R")
#biocLite(c("parallel", "doParallel", "glmnet"), lib="/nfs/users/nfs_t/ti1/r_libs")
.libPaths( c( .libPaths(), "/users/spjtcoi/brc_scratch/Rlibs") )
library(parallel)#package 'parallel' is a base package, and should not be updated
library(doParallel)
library(glmnet)
library(purrr) #for map function for lists. See tutorial by Hadley Wickham: http://r4ds.had.co.nz/lists.html
library(ggplot2)


#Parallelise foreach onto multiple cores on a SINGLE node
register_parallel<-function(){
  #register the workers
  cores=detectCores()

  cl <- makeCluster(cores-1)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("/users/spjtcoi/brc_scratch/Rlibs"))#Register libraries onto each core
  return (cl)
}

#RUN GLMNET FOR n iterations and all parameters alpha on our dataset. 
run_glmnet<-function(target, features, alphas, iterations, random=FALSE){
  
  #RANDOMISE TARGET IF NECCESSARY
  if(random){ 
    target = sample(target, length(target), replace = FALSE)
  }
  
  #Target variable - what we want to predict
  #These are our features that we train the algorithm with
  cl = register_parallel()
  
  folds<-sample(rep(seq(10), length.out=nrow(features)))#Create f folds
  folds <- unique(t(replicate(10000, sample(folds)))) #Create 10FoldCV a 10000 times
  folds_matrix<-folds[sample(nrow(folds), iterations), ] #100 iterations of 10 folds
  
  #For each alpha in the list of alphas, .export the new fold dataset we want to use
  #For each iteration, copy the variables via .export to each CPU
  glmnet_results = foreach(j = 1:iterations)  %:% foreach(i = 1:length(alphas), .inorder = FALSE,.packages = "glmnet", .multicombine = TRUE,.export= c("target","alphas","features", "folds_matrix")) %dopar% {
    folds = folds_matrix[j, ]
    #Run GLMNET witt our features and our target variable
    output = cv.glmnet(
      x = features,
      y = target,
      family = "gaussian",
      type.measure = "mse",
      foldid = folds)#use system lambdas
    #Return GLM-Net results
    return(output)
  }
  stopCluster(cl)
  #Return GLMNET results for each iteration, each alpha
  return(glmnet_results)
}




extractGlmnetInfo <-function(object) {
  #find lambdas
  lambdaMin <- object$lambda.min
  lambda1se <- object$lambda.1se
  
  #figure out where thise lambdas fall in the path
  whichMin <- which(object$lambda ==lambdaMin)
  which1se <- which(object$lambda ==lambda1se)
  
  coefficents = as.data.frame(as.matrix(coef(object, s = "lambda.min")))
  coef_features = rownames(coefficents)
  negative_class = which(coefficents < 0)
  postitive_class = which(coefficents > 0)
  negative_features = coef_features[negative_class]
  positive_features = coef_features[postitive_class]
  
  #build a one line data.frame with each of the selected lambdas and
  #its corresponding error figures
  df=data.frame(lambda.min=lambdaMin, error.min=object$cvm[whichMin],
             lambda.1se=lambda1se, error.1se=object$cvm[which1se])
  
  return (df)
}


extractFeatureImportance <-function(object) {

  coefficents = as.data.frame(as.matrix(coef(object, s = "lambda.min")))
  coef_features = rownames(coefficents)
  negative_class = which(coefficents < 0)
  postitive_class = which(coefficents > 0)
  negative_features = coef_features[negative_class]
  positive_features = coef_features[postitive_class]
  
  feature_list = c(positive_features, negative_features)

  #build a one line data.frame with each of the selected lambdas and
  #its corresponding error figures
  df=data.frame(feature_name = feature_list, class_importance = c(rep("positive", length(positive_features)),rep("negative", length(negative_features))))
  return (df)
  
}


feature_importance_stats <- function(model_results) {
  
  exp<-function(x) Reduce(rbind, lapply(x, extractFeatureImportance))
  feature_importance<-map(model_results, exp) #Extract imputation sets one at a time
  feature_importance = do.call("rbind", feature_importance)
  feature_importance_freq = table(feature_importance)
  
  total_runs =  length(model_results)*length(model_results[[1]])
  feature_importance_freq_relative = feature_importance_freq/total_runs
  return(list(feature_importance_freq, feature_importance_freq_relative))
}


#GETTING THE RESULTS
summariseGLMNETResults <- function (full_results) {
  exp<-function(x) Reduce(rbind, lapply(x, extractGlmnetInfo))
  ExtractedSet.tmp<-map(full_results, exp) #Extract imputation sets one at a time
  full_results = do.call("rbind", ExtractedSet.tmp)
  alphas <- seq(from = 0.5, to = 1, by = 0.05)
  full_results = cbind(full_results,rep.int(alphas, 100))
  full_results = cbind(full_results,sort(rep(1:100, 11)))
  
  full_results = cbind(full_results, rep(1:11, 100))
  colnames(full_results) = c("lambda_min", "error_min", "lambda_dev", "error_dev", "alpha", "iteration", "parameter_combi")
  return(full_results)
}

#Generate counts matrix for feature importance
generate_counts<-function(feat_import_orig, feat_import_perm){
  original_fi_counts = (feat_import_orig[[1]]) #contains both positive ad negative
  permutation_fi_counts = (feat_import_perm[[1]]) #contains both positive ad negative
  original_fi_counts_df <-convert_counts_to_dataframe(original_fi_counts)
  permutation_fi_counts_df <-convert_counts_to_dataframe(permutation_fi_counts)
  all_fi_counts<-merge(x = original_fi_counts_df, y = permutation_fi_counts_df, by = "Genes", all.x = TRUE)
	colnames(all_fi_counts) = c("Genes", "Positive_real", "Negative_real", "Positive_Permutation", "Negative_Permutation")
  all_fi_counts[is.na(all_fi_counts)] <- 0
  all_fi_counts=all_fi_counts[c("Genes", "Positive_real", "Positive_Permutation", "Negative_real", "Negative_Permutation")]
  return(all_fi_counts)
}

#Generate counts df for feature importance
convert_counts_to_dataframe <- function(data) {


  col_names=(colnames(data))
  if (!"negative" %in% col_names  || !"positive" %in% col_names){
	positive = data[,1]
	data = cbind(data, rep(0, length(positive)))
 } 

  df = data.frame(matrix(data, ncol = 2))
	
  df["Genes"] = rownames(data)

  colnames(df) = c("Positive", "Negative", "Genes")
  return(df)
}

#Calculate p value for a single gene
calc_pval_per_gene <- function(t){ 
  
  pr = t[["Positive_real"]]
  pp = t[["Positive_Permutation"]]
  
  nr = t[["Negative_real"]]
  np = t[["Negative_Permutation"]]
  
  positive_case <-
    matrix(c(pr, 1100-pr, pp, 1100-pp),
           nrow = 2,
           dimnames =
             list(c("Effect", "No Effect"),
                  c("Real Data", "Permutations")))
  
  negative_case <-
    matrix(c(nr, 1100-nr, np, 1100-np),
           nrow = 2,
           dimnames =
             list(c("Effect", "No Effect"),
                  c("Real Data", "Permutations")))
  
    positive_f = fisher.test(positive_case, alternative = "greater")$p.value
  negative_f = fisher.test(negative_case, alternative = "greater")$p.value
  return(c(positive_f, negative_f))
}

#Calculate p value for all gene
calc_pval_for_genelist<-function(x){

  remove_i=which(x["Genes"] == c("(Intercept)"))
  p_values = t(apply(x[, -1], 1, calc_pval_per_gene))
  
  positive_binary=p_values[,1]<0.05
  negative_binary=p_values[,2]<0.05
  
  x["positive_binary"] =  positive_binary
  x["negative_binary"] =  negative_binary
  
  x["positive_p_value"] =  p_values[,1]
  x["negative__value"] =  p_values[,2]
  return(x)
}

generate_predictions <- function(data, feature_list, features_to_remove, output_dir, number) {
  
  before <- Sys.time()  #keep track of timing
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)
  
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
    data = data[, -which(colnames(data) %in% remove)]
  }
  
  f = feature_list[feature_list[, 2] == "target",]
  target_column =  as.character(f$feature_name)
  print(paste0('This is my target_column: ',target_column))
  target = as.matrix(data[, which(colnames(data) %in% target_column)])
  features = as.matrix(data[, -which(colnames(data) %in% target_column)])

print("Running predictions")
  
  alphas <- seq(from = 0.5, to = 1, by = 0.05)
  results = run_glmnet(target,
                       features,
                       alphas,
                       100)
  
  full_results = summariseGLMNETResults(results)
  write.table(full_results, paste0(output_dir, "/accuracy", number, ".txt"),  col.names =FALSE, row.names=FALSE)
  
  mean_per_parameter_combi = data.frame(do.call("rbind", sapply(c(
    unique(full_results["parameter_combi"])
  )[[1]], function (x) {
    apply(full_results[subset(full_results["parameter_combi"] == x), ], 2, mean)
  }, simplify = FALSE)))


  #Feature importance calculations
  original_fi = feature_importance_stats(results)#extraction of transcript-level counts
print("Running permutated predictions")
#We want to increase permutations to probably a 1000
  p_result = run_glmnet(target,features,alphas,
                                         iterations = 100, random = TRUE)
  permutation_fi = feature_importance_stats(p_result) #extaction of transcript-level
  count_results <- generate_counts(original_fi, permutation_fi) #consolidation of

  #print(count_results)
  pvals_all_genes <- calc_pval_for_genelist(count_results)
  write.table(pvals_all_genes, paste0(output_dir, "/pvalues.", number), sep="\t")

  after<-Sys.time() #stop timing
  print(after-before) #time difference(depends on speed, mem and number of cores available on the machine)
}

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc

#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 
#args[1] = "/Users/ti1/Documents/data_test/dataset_3.csv"
#args[2] = "/Users/ti1/Downloads/all_features2.csv"
#args[3] = "/Users/ti1/Downloads/demo_features_to_remove.csv"
#args[4] = "~/Downloads/my_dir_test"
print(paste0('GLMNET analaysis'))
print(paste0('Dataset:', args[1]))
print(paste0('Feature reference: ', args[2] ))
print(paste0('Features to remove: ', args[3] ))
print(paste0('Output dir: ', args[4] ))

#data = read.csv(args[1], stringsAsFactors = FALSE, sep=' ' )# for space separated file 
data = read.csv(args[1], stringsAsFactors = FALSE)#if comma delimited the remove

#This is the fike which contains all features and its corresponding class
#Column 1: Feature name
#Column 2: Feature class
if(!file.exists(args[2])) {
  stop("Feature reference file not found. Please check.")
} else {
  feature_list = read.csv(args[2], stringsAsFactors = FALSE)
}

#This file lists all classes and/or features that should be removed. It is a 1-dim list.
#If this file doesn't exist, nothing will be removed
if(file.exists(args[3])) {
  features_to_remove = read.csv(args[3], stringsAsFactors = TRUE, header=FALSE)
} else {
  features_to_remove = c()
}
number = Sys.getenv(c("SGE_TASK_ID"))

generate_predictions(data=data, feature_list, features_to_remove, output_dir = args[4], number)
