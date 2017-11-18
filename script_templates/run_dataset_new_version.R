#source("http://bioconductor.org/biocLite.R")
#biocLite(c("parallel", "doParallel", "glmnet"), lib="/nfs/users/nfs_t/ti1/r_libs")
.libPaths( c( .libPaths(), "/nfs/users/nfs_t/ti1/r_libs/") )
library(parallel)#package 'parallel' is a base package, and should not be updated
library(doParallel)
library(glmnet)
library(purrr) #for map function for lists. See tutorial by Hadley Wickham: http://r4ds.had.co.nz/lists.html
library(ggplot2)


prepare_feature_sets<-function(data_list, remove_features, imputed_set_number) {
  
  gene_expression = subset(data_list$bmi_enet, select=-c(gender, ChipID_2B,  ChipID_2C,  ChipID_2D,  ChipID_2E,  ChipID_2F, ChipID_2G,  ChipID_2H,  ChipID_2I,  ChipID_2J,  ChipID_2K,  ChipID_2L,  Timepoint, 
                                                              sv4, sv5, hlth, ovwgt, obese, Age.norm, PC1, PC2, med.days, GLUCOSE_GRS, OBSESITY_GRS))
  
  technical_cofounders = subset(data_list$bmi_enet, select=c(ChipID_2B,  ChipID_2C,  ChipID_2D,  ChipID_2E,  ChipID_2F, ChipID_2G,  ChipID_2H,  ChipID_2I,  ChipID_2J,  ChipID_2K,  ChipID_2L,  Timepoint, 
                                                           sv4, sv5))
  
  imputed_data = data_list$imputed_sets[[as.numeric(imputed_set_number)]]
  target = imputed_data[, ('BMI.norm')]
  imputed_features = imputed_data[ , -which(colnames(imputed_data) %in% c("BMI.norm"))]
  
  personal_cofounders = subset(data_list$bmi_enet, select=c(gender))
  personal_cofounders = cbind(personal_cofounders,imputed_features)
  
  if (length(remove_features) > 0) {
    
    print("Removing following features from the dataset:")
    #print(paste0(remove_features, sep=""))
    
    print(paste0("personal_cofounders before: ", ncol(personal_cofounders)))
    print(paste0("technical_cofounders before: ", ncol(technical_cofounders)))
    print(paste0("gene_expression before: ", ncol(gene_expression)))
    
    personal_cofounders = personal_cofounders[,!colnames(personal_cofounders) %in% remove_features]
    technical_cofounders = technical_cofounders[,!colnames(technical_cofounders)  %in% remove_features]
    gene_expression = gene_expression[,!colnames(gene_expression)  %in% remove_features]
    
    print(paste0("personal_cofounders after: ", ncol(personal_cofounders)))
    print(paste0("technical_cofounders after: ", ncol(technical_cofounders)))
    print(paste0("gene_expression after: ", ncol(gene_expression)))
    
  }
  return(list(target=target, personal_cofounders=personal_cofounders, technical_cofounders=technical_cofounders, gene_expression=gene_expression))
}

#Parallelise foreach onto multiple cores on a SINGLE node
register_parallel<-function(){
  #register the workers
  cores=detectCores()

  print("Number of cores:")
  print(cores)
  cl <- makeCluster(cores-1)
  registerDoParallel(cl)
  clusterEvalQ(cl, .libPaths("/nfs/users/nfs_t/ti1/r_libs"))#Register libraries onto each core
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
      x =features,
      y = target,
      family = "gaussian",
      type.measure = "mse")#use system lambdas
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

process_single_imputeset <- function(input_data, imputed_set_number, remove_features, output_dir) {
  
  before <- Sys.time()  #keep track of timing
  dir.create(output_dir, showWarnings = TRUE, recursive = TRUE)
  
  data_list = prep_data(
    input_data$imp_list_Tall,
    input_data$GX_sva,
    input_data$svmod_bmi_catg_sv,
    input_data$svmod_hba1c_catg_sv,
    input_data$pheno_mini,
    input_datatpheno
  )#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file
  #Set the folds for each BMI
   
  prepared_data = prepare_feature_sets(data_list, remove_features, 1)
  alphas <- seq(from = 0.5, to = 1, by = 0.05)
  
  features = cbind(prepared_data$personal_cofounders, prepared_data$gene_expression, prepared_data$technical_cofounders) 
  #features = cbind(sample(c(2323:100000),size = 97), sample(c(2323:5000),size = 97))
  results = run_glmnet(prepared_data$target,
                       features,
                       alphas,
                       100)
  
  full_results = summariseGLMNETResults(results)
  write.table(full_results, paste0(output_dir, "/accuracy", imputed_set_number, ".txt"),  col.names =FALSE, row.names=FALSE)
 
  
  mean_per_parameter_combi = data.frame(do.call("rbind", sapply(c(
    unique(full_results["parameter_combi"])
  )[[1]], function (x) {
    apply(full_results[subset(full_results["parameter_combi"] == x), ], 2, mean)
  }, simplify = FALSE)))

  #Feature importance calculations
  original_fi = feature_importance_stats(results)#extraction of transcript-level counts
  #We want to increase permutations to probably a 1000
  p_result = run_glmnet(prepared_data$target,features,alphas,
                                         iterations = 100, random = TRUE)
  permutation_fi = feature_importance_stats(p_result) #extaction of transcript-level
  count_results <- generate_counts(original_fi, permutation_fi) #consolidation of

  #print(count_results)
  pvals_all_genes <- calc_pval_for_genelist(count_results)
  write.table(pvals_all_genes, paste0(output_dir, "/pvalues.", imputed_set_number), sep="\t")

  after<-Sys.time() #stop timing
  print(after-before) #time difference(depends on speed, mem and number of cores available on the machine)
}


#Read data from server
read_data <-function(my_dir, input_filename){ 
  setwd(my_dir)
  imputation_data_generated = load(input_filename)
  data = list(hba1c_descrip=hba1c.descrip, GX_sva=GX.sva, gene_symbols=gene_symbols, random_seed=.Random.seed, Tri_revert=TriG.revert,
              svmod_bmi_catg_sv=svmod.bmi.catg.sv, svmod_hba1c_catg=svmod.hba1c.catg, t_expr_fep=t.expr.fep, fglu_revert=fglu.revert, 
              Age_revert=Age.revert, Height_revert=Height.revert,
              svmod_hba1c_catg_sv=svmod.hba1c.catg.sv,  pump_imp=pump.imp, gx_prejoin_1=gx.prejoin.1, pheno_sva=pheno.sva, hba1c_revert=hba1c.revert,
              imp_out=imp.out, imp_out_all_list=imp.out.all.list,
              BMI_revert=BMI.revert, svobj_hba1c_catg=svobj.hba1c.catg, svmod_bmi_catg=svmod.bmi.catg, svobj_bmi_catg=svobj.bmi.catg, pheno_mini=pheno.mini, 
              test_obj=test.obj, imp_list_Tall=imp.list.Tall)
  
  return(data)
}

#input vars
#Inplicit way to pass parameters over positions
#data_list=prep_data(imp.list.Tall, GX.sva, svmod.bmi.catg.sv, svmod.hba1c.catg.sv, pheno.mini, )#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file
#NOT OK TO DO THIS: THE FUNCTION WILL ASSUME THAT covars IS GX.sva, because it looks at the position
#data_list=prep_data( GX.sva, imp.list.Tall, svmod.bmi.catg.sv, svmod.hba1c.catg.sv, pheno.mini, )#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file

#Explicit way to pass parameters
#data_list=prep_data(covars=imp.list.Tall, sva_file=GX.sva, bmi_sv_file=svmod.bmi.catg.sv, hba_sv_file=svmod.hba1c.catg.sv, pheno=pheno.mini, tpheno=)#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file\
#This is the same as the line before, as we explictily assign the variables to the parameters so we 
#are allowed to change the order 
#data_list=prep_data( sva_file=GX.sva, covars=imp.list.Tall,hba_sv_file=svmod.hba1c.catg.sv, bmi_sv_file=svmod.bmi.catg.sv, tpheno=, pheno=pheno.mini)#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file

prep_data<-function(covars, sva_file, bmi_sv_file, hba_sv_file, pheno, tpheno ){
  test.obj2<-lapply(covars,`[`, c("ChipID", "PC1", "PC2", "Age.norm", "BMI.norm", "GLUCOSE_GRS", "OBSESITY_GRS", "med.days"))#DELETE FOLLOWING VARS AS APPROPRIATE (ACCORDING TO TRAIT)
  #2=CHipID,9=PC1,10=PC2,15=Age.norm,16=BMI.norm,20=hba.norm,43=BMI.obese,52=hba.dbx, 24=GLUCOSE_GRS, 25=OBSESITY_GRS, 14=med.days  #output_check<-test.obj2[[1]]
  myvec <- sapply(test.obj2, NROW)#check rownums make sense and of equal length. have dup sample ids been removed?
  tpheno<-t(pheno)#create list of ChipIDs with same sequence as 'pheno'. To be used for post-mice filtering
  tpheno.list <- split(tpheno, seq(nrow(tpheno))) #http://stackoverflow.com/questions/3492379/data-frame-rows-to-a-list
  list.97<-tpheno.list[[2]]
  rm(tpheno, tpheno.list)
  #final househeeping on test.obj
  abc<-lapply(test.obj2, function(x) x[x$ChipID %in% list.97,]) #subset list of dataframes based on list of 90 final IDs
  abc<-lapply(abc, function(x) x[!duplicated(x["ChipID"]), ]) #remove duplicate ChipIDs
  abc<-lapply(abc, function(x) x[!(names(x) %in% "ChipID")]) #remove column in listed dataframe:http://stackoverflow.com/questions/12664430/delete-a-column-in-a-data-frame-within-a-list
  abc<-lapply(abc, function(x) sapply(x, function(x) unclass(x))) # r base function #http://stackoverflow.com/questions/19836464/apply-subset-function-to-a-list-of-dataframes , http://stackoverflow.com/questions/22407494/remove-duplicates-from-list-elements
  #Creation of Elastic Net data object
  Gx.t<-t(sva_file) #Check obs order is preserved via: Gx.t<-cbind(rows=rownames(Gx.t), Gx.t)
  bmi.enet = cbind(bmi_sv_file, Gx.t)
  hba.enet = cbind(hba_sv_file, Gx.t)
  
  return(list(imputed_sets=abc, bmi_enet=bmi.enet, hba_enet=hba.enet)) #returning multiple things requires wrapping as a list
}


#Read input data containing everything 
i = Sys.getenv(c("LSB_JOBINDEX"))

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc

#This is just for testing purposes. If you don't make it as comment, it will overwrite any values that you have given over the command line 
args[1] = "~/Downloads/my_dir_test"
#args[2] = "/Users/ti1/Google\ Drive/results/default_accuracy_analysis_gene_only_removal/remove_top__random_100.csv"
#args[2] = ""

print(paste0('GLMNET analaysis on dataset:', i))
print(paste0('Output directory: ', args[1] ))
print(paste0('Removing: ', args[2] ))#

#If file doesnt exist or nothing has been specified, assume no features should be removed
if(!file.exists(args[2])) {
  remove=c()
} else {
  remove=read.csv(args[2], stringsAsFactors = FALSE)[,2]
}

input_data=read_data("/Users/ti1/Google\ Drive/raw\ data", "imputation_data_generated.RData")
process_single_imputeset(input_data, i, remove, output_dir = args[1])
