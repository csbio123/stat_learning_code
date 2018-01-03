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
  test.obj2<-lapply(covars,`[`, c("ChipID", "PC1", "PC2", "Age.norm", "BMI.norm"))#DELETE FOLLOWING VARS AS APPROPRIATE (ACCORDING TO TRAIT)
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



prepare_feature_sets<-function(data_list, outptut_dir, prefix="dataset") {
  
  
  
  gene_expression =data_list$bmi_enet[,60:ncol(data_list$bmi_enet)]
  
  technical_cofounders = subset(data_list$bmi_enet, select=c(sv1, sv2, sv3, sv4, sv5, sv6))
  imputed_data = data_list$imputed_sets
  personal_cofounders = subset(data_list$bmi_enet, select=c(MaleGender)) #decided to use gender alongside SVs4&5 as these SVs are more technical than demographic
  dir.create(outptut_dir, showWarnings = TRUE, recursive = TRUE)
  outptut_file_pre = paste0(outptut_dir, "/" , prefix)
  print(length(imputed_data))
  sapply(seq(1, length(imputed_data)), function(x) {
    
    imp_d = imputed_data[[x]]
    features = cbind(imp_d, personal_cofounders, technical_cofounders, gene_expression)
    print(dim(features))
    outptut_file = paste0(outptut_file_pre, "_", x, ".csv")
    write.csv(features, file=outptut_file)
  })
  
}

input_data=read_data("/Users/ti1/Google Drive/raw data/", "generate_metabolic_features_ALL_SVs_3catg_validationset.RData")

data_list = prep_data(
  input_data$imp_list_Tall,
  input_data$GX_sva,
  input_data$svmod_bmi_catg_sv,
  input_data$svmod_hba1c_catg_sv,
  input_data$pheno_mini)#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file
#Set the folds for each BMI

prepared_data = prepare_feature_sets(data_list, outptut_dir = "/Users/ti1/Google Drive/raw data/validation_data_regression")
names(prepared_data)
