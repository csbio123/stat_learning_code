library(RUVSeq)
library(sva)


#Read data from server
read_data <-function(my_dir, input_filename){ 
  setwd(my_dir)
  imputation_data_generated = load(input_filename)
  data = list(hba1c_descrip=hba1c.descrip, GX_sva=GX.sva, gene_symbols=gene_symbols, random_seed=.Random.seed, 
              svmod_bmi_catg_sv=svmod.bmi.catg.sv, svmod_hba1c_catg=svmod.hba1c.catg,
          
              svmod_hba1c_catg_sv=svmod.hba1c.catg.sv, gx_prejoin_1=gx.prejoin.1, pheno_sva=pheno.sva,
              imp_out=imp.out, imp_out_all_list=imp.out.all.list,
              svobj_hba1c_catg=svobj.hba1c.catg, svmod_bmi_catg=svmod.bmi.catg, svobj_bmi_catg=svobj.bmi.catg, pheno_mini=pheno.mini, 
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

prep_data<-function(covars, sva_file, bmi_sv_file, hba_sv_file, pheno ){
  test.obj2<-lapply(covars,`[`, c("ChipID", "PC1", "PC2", "Age.norm", "BMI.catg", "BMI.norm"))#DELETE FOLLOWING VARS AS APPROPRIATE (ACCORDING TO TRAIT)
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





test_data=read_data("/Users/ti1/Google\ Drive/conrad/dataset2_chronic_scz/", "generate_metabolic_features_ALL_SVs_2catg_validationset_09-01-18.RData")

data_list_test = prep_data(
  test_data$imp_list_Tall,
  test_data$GX_sva,
  test_data$svmod_bmi_catg_sv,
  test_data$svmod_hba1c_catg_sv,
  test_data$pheno_mini)#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file
#Set the folds for each BMI


train_data=read_data("/Users/ti1/Google\ Drive/raw\ data/", "ALL_SV6_drugnaive_hba_chip_2factor_01-09-17.RData")

data_list_train = prep_data(
  train_data$imp_list_Tall,
  train_data$GX_sva,
  train_data$svmod_bmi_catg_sv,
  train_data$svmod_hba1c_catg_sv,
  train_data$pheno_mini
)

gene_expression_train = (train_data$GX_sva)
colnames(gene_expression_train) = paste0("train_", colnames(gene_expression_train))
category_train = as.factor(data_list_train$imputed_sets[[1]][,4]) #category ==condition eg bmi_catg

gene_expresion_test = (test_data$GX_sva)
colnames(gene_expresion_test) = paste0("test_", colnames(gene_expresion_test))
category_test = as.factor(data_list_test$imputed_sets[[1]][,4]) #category ==condition eg bmi_catg

##BATCH CORRECTION##
removing_unwanted_variation<-function(category, gene_expression_train, gene_expresion_test) {
  
  same_genes = intersect( row.names(gene_expression_train), row.names(gene_expresion_test))#intersection of genes in both training and test
  
  train_common = gene_expression_train[rownames(gene_expression_train) %in% same_genes,] # select intersecting genes
  test_common  = gene_expresion_test[rownames(gene_expresion_test) %in% same_genes,]
  
  sample_num = 1:length(category)
  sample_name = colnames(gene_expression_train)
  train_pheno = data.frame(sample_num, condition=category, row.names = sample_name)
  
  trainMod = model.matrix(~as.factor(condition), data=train_pheno) #Within model matrix the main phenotye is the condition
  trainMod0 = model.matrix(~1, data=train_pheno)
  trainSv = sva(train_common, trainMod, trainMod0)# apply sva to training data
  
  fsvaobj = fsva(train_common, trainMod, trainSv, test_common)#provide parameters for training sample and then matching test geneset
  mydataSv = list(x=fsvaobj$db, y=train_pheno$cat)# col1==residualised gene expression for test geneset col2 = phenotypic status
  
  train_common_norm = mydataSv$x # residualised gene expression for test geneset 
  test_common_norm = fsvaobj$new # corrected (residulaised) version of test geneset
  
  return(list(train_common_norm, test_common_norm))
}


prepare_feature_sets<-function(data_list, gene_expression, outptut_dir, prefix="dataset") {
  
  col_names = colnames(data_list$bmi_enet)
  colnames(data_list$bmi_enet) = gsub("MaleGender", "gender", col_names) # Data-specific part!
  clinical_cofounders = subset(data_list$bmi_enet, select=c(hba1c.norm))
  imputed_data = data_list$imputed_sets
  personal_cofounders = subset(data_list$bmi_enet, select=c(gender)) #decided to use gender alongside SVs4&5 as these SVs are more technical than demographic
  dir.create(outptut_dir, showWarnings = TRUE, recursive = TRUE)
  outptut_file_pre = paste0(outptut_dir, "/" , prefix)
  print(length(imputed_data))
  
  outptut_file = paste0(outptut_file_pre, "_gene_expression.csv")#gene expression outputted separately
  write.csv(t(gene_expression), file=outptut_file)
  
  sapply(seq(1, length(imputed_data)), function(x) {
    imp_d = imputed_data[[x]]
    features = cbind(imp_d, personal_cofounders, clinical_cofounders)
    print(dim(features))
    outptut_file = paste0(outptut_file_pre, "_", x, ".csv")
    write.csv(features, file=outptut_file)
  })
  
}



##BATCH CORRECTION##
removing_unwanted_variation_sva<-function(category, gene_expression_train, gene_expresion_test) {
  
  same_genes = intersect( row.names(gene_expression_train), row.names(gene_expresion_test))#intersection of genes in both training and test
  
  train_common = gene_expression_train[rownames(gene_expression_train) %in% same_genes,] # select intersecting genes
  test_common  = gene_expresion_test[rownames(gene_expresion_test) %in% same_genes,]
  
  sample_num = 1:length(category)
  sample_name = colnames(gene_expression_train)
  train_pheno = data.frame(sample_num, condition=category, row.names = sample_name)
  
  trainMod = model.matrix(~as.factor(condition), data=train_pheno) #Within model matrix the main phenotye is the condition
  trainMod0 = model.matrix(~1, data=train_pheno)
  trainSv = sva(train_common, trainMod, trainMod0)# apply sva to training data
  
  fsvaobj = fsva(train_common, trainMod, trainSv, test_common)#provide parameters for training sample and then matching test geneset
  mydataSv = list(x=fsvaobj$db, y=train_pheno$cat)# col1==residualised gene expression for test geneset col2 = phenotypic status
  
  train_common_norm = mydataSv$x # residualised gene expression for test geneset 
  test_common_norm = fsvaobj$new # corrected (residulaised) version of test geneset
  
  return(list(train_common_norm, test_common_norm))
}


##BATCH CORRECTION##
removing_unwanted_variation_combat<-function(category_train, category_test, gene_expression_train, gene_expresion_test) {
  
  
  same_genes = intersect( row.names(gene_expression_train), row.names(gene_expresion_test))#intersection of genes in both training and test
  train_common = gene_expression_train[rownames(gene_expression_train) %in% same_genes,] # select intersecting genes
  test_common  = gene_expresion_test[rownames(gene_expresion_test) %in% same_genes,]
  
  sample_num = 1:length(colnames(gene_expression_train))
  sample_name = colnames(gene_expression_train)
  batches_train = gsub("_.*", "",gsub("train_", "",sample_name))
  train_pheno = data.frame(sample_num,  row.names = sample_name, condition=category_train)
  #print(train_pheno)
  
  sample_num = 1:length(colnames(gene_expresion_test))
  sample_name = colnames(gene_expresion_test)
  batches_test = gsub("_.*", "",gsub("test_", "",sample_name))
  test_pheno = data.frame(sample_num, row.names = sample_name, condition=category_test)
  
  modcombat = model.matrix(~as.factor(condition), data=train_pheno)
  combat_normalized_train = ComBat(dat=gene_expression_train, batch=batches_train, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  modcombat = model.matrix(~as.factor(condition), data=test_pheno)
  combat_normalized_test = ComBat(dat=gene_expresion_test, batch=batches_test, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  return(list(combat_normalized_train, combat_normalized_test))
}


##BATCH CORRECTION##
removing_unwanted_variation_combat_multiple_batch<-function(category_train, category_test, gene_expression_train, gene_expresion_test, batch_train, batch_test) {
  
  
  same_genes = intersect( row.names(gene_expression_train), row.names(gene_expresion_test))#intersection of genes in both training and test
  train_common = gene_expression_train[rownames(gene_expression_train) %in% same_genes,] # select intersecting genes
  test_common  = gene_expresion_test[rownames(gene_expresion_test) %in% same_genes,]
  
  sample_num = 1:length(colnames(gene_expression_train))
  sample_name = colnames(gene_expression_train)
  batches_train = gsub("_.*", "",gsub("train_", "",sample_name))
  train_pheno = data.frame(sample_num,  row.names = sample_name, condition=category_train, date=batch_train)
  #print(train_pheno)
  
  sample_num = 1:length(colnames(gene_expresion_test))
  sample_name = colnames(gene_expresion_test)
  batches_test = gsub("_.*", "",gsub("test_", "",sample_name))
  test_pheno = data.frame(sample_num, row.names = sample_name, condition=category_test)
  
  modcombat = model.matrix(~as.factor(condition), data=train_pheno)
  combat_normalized_train = ComBat(dat=gene_expression_train, batch=batches_train, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  modcombat = model.matrix(~as.factor(condition), data=test_pheno)
  combat_normalized_test = ComBat(dat=gene_expresion_test, batch=batches_test, mod=modcombat, par.prior=TRUE, prior.plots=FALSE)
  
  return(list(combat_normalized_train, combat_normalized_test))
}




batch_test = "/Users/ti1/Google\ Drive/validation_prediction/batch_test.csv"
batch_train = "/Users/ti1/Google\ Drive/validation_prediction/batch_train.csv"

batch_test = read.csv(batch_test)
batch_train = read.csv(batch_train)


datasets_sva = removing_unwanted_variation_sva(category_train, gene_expression_train, gene_expresion_test)
datasets_combat = removing_unwanted_variation_combat(category_train, category_test, gene_expression_train, gene_expresion_test)
#datasets_combat = removing_unwanted_variation_combat_multiple_batch(category_train, category_test, gene_expression_train, gene_expresion_test, batch_train, batch_test)


prepare_feature_sets(data_list_test, datasets_sva[[2]], "/Users/ti1/Downloads/validation_prediction/validation_sva_9_03_18", "validation_set")
prepare_feature_sets(data_list_train, datasets_sva[[1]], "/Users/ti1/Downloads/validation_prediction/train_sva_9_03_18", "training_set")

prepare_feature_sets(data_list_test, datasets_combat[[2]], "/Users/ti1/Downloads/validation_prediction/validation_combat_9_03_18", "validation_set")
prepare_feature_sets(data_list_train, datasets_combat[[1]], "/Users/ti1/Downloads/validation_prediction/train_combat_9_03_18", "training_set")



