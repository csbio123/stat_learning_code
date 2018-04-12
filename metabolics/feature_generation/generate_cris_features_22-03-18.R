#Read data from server
read_data <-function(my_dir, input_filename){ 
  setwd(my_dir)
  imputation_data_generated = load(input_filename)
  data = list(GX_sva=GX.sva, gene_symbols=gene_symbols, svmod_warddays_sv=svmod.warddays.sv, svmod_referraldays_sv=svmod.referraldays.sv, t_expr_fep=t.expr.fep, cris_data_imp=cris_data.imp, gx_prejoin_1=gx.prejoin.1, pheno_sva=pheno.sva, imp_out=imp.out, imp_out_all_list=imp.out.all.list, svobj_referraldays=svobj.referraldays, svobj_warddays=svobj.warddays, pheno_mini=pheno.mini, test_obj=test.obj, imp_list_Tall=imp.list.Tall)

  
  return(data)
}


#input vars
#Inplicit way to pass parameters over positions
#data_list=prep_data(imp.list.Tall, GX.sva, svmod.bmi.catg.sv, svmod.hba1c.catg.sv, pheno.mini, )#this works because order in function definition is maintained ie.covars, sva_file, warddays_file, referral_file
#NOT OK TO DO THIS: THE FUNCTION WILL ASSUME THAT covars IS GX.sva, because it looks at the position
#data_list=prep_data( GX.sva, imp.list.Tall, svmod.bmi.catg.sv, svmod.hba1c.catg.sv, pheno.mini, )#this works because order in function definition is maintained ie.covars, sva_file, warddays_file, referral_file

#Explicit way to pass parameters
#data_list=prep_data(covars=imp.list.Tall, sva_file=GX.sva, warddays_file=svmod.bmi.catg.sv, referral_file=svmod.hba1c.catg.sv, pheno=pheno.mini, tpheno=)#this works because order in function definition is maintained ie.covars, sva_file, warddays_file, referral_file\
#This is the same as the line before, as we explictily assign the variables to the parameters so we 
#are allowed to change the order 
#data_list=prep_data( sva_file=GX.sva, covars=imp.list.Tall,referral_file=svmod.hba1c.catg.sv, warddays_file=svmod.bmi.catg.sv, tpheno=, pheno=pheno.mini)#this works because order in function definition is maintained ie.covars, sva_file, warddays_file, referral_file

prep_data<-function(covars, sva_file, warddays_file, referral_file, pheno, tpheno ){
  test.obj2<-lapply(covars,`[`, c("ChipID", "african", "Indian", "Age.norm", "referraldays", "med.days", "Med_complian"))#DELETE FOLLOWING VARS AS APPROPRIATE (ACCORDING TO TRAIT)
  #2=CHipID,9=african,10=Indian,15=Age.norm,16=referraldays,20=hba.norm,43=BMI.obese,52=hba.dbx, 24=GLUCOSE_GRS, 25=OBSESITY_GRS, 14=med.days  #output_check<-test.obj2[[1]]
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
  warddays.enet = cbind(warddays_file, Gx.t)
  referral.enet = cbind(referral_file, Gx.t)
  
  return(list(imputed_sets=abc, warddays_enet=warddays.enet, referral_enet=referral.enet)) #returning multiple things requires wrapping as a list
}



prepare_feature_sets<-function(data_list, outptut_dir, prefix="dataset") {
  
  gene_expression = subset(data_list$warddays_enet, select=-c(warddays, Age.norm, gender, african, Indian, ChipID_19020374058, ChipID_19020374069, ChipID_19020374072,ChipID_19020374078, ChipID_19020374079, ChipID_19031356054, ChipID_19031356056, ChipID_19031356062, ChipID_19031356068, ChipID_19031356070, ChipID_1903135607, ChipID_19216457008, ChipID_19216457009, ChipID_19216457012, ChipID_19216457014, ChipID_19216457025, ChipID_19216457029, ChipID_19216457032, ChipID_19216457033,ChipID_19234921047, ChipID_19234921059, ChipID_19234921061, ChipID_19234921063, ChipID_19234921065, ChipID_19234921066, ChipID_19234921070, ChipID_19234921074,ChipID_19234921077, ChipID_19234921082, ChipID_19234921083, ChipID_19234921090, ChipID_19234921094, ChipID_19234921100, ChipID_19234921101, ChipID_19235792061,ChipID_19235792082, ChipID_19235792088, ChipID_19235792089, ChipID_19235792091, ChipID_19235792095, ChipID_19249896045, ChipID_19249896067, ChipID_19249896073,ChipID_19249896091, ChipID_19249907011, ChipID_19249907021, ChipID_19249907031, ChipID_19249907044, ChipID_19249907045, ChipID_19249907052, med.days, Med_complian1, Med_complian2, Med_complian3, SV1, SV2, SV4, SV5, SV6, SV7))
  
  technical_cofounders = subset(data_list$warddays_enet, select=c(ChipID_19020374058, ChipID_19020374069, ChipID_19020374072,ChipID_19020374078, ChipID_19020374079, ChipID_19031356054, ChipID_19031356056, ChipID_19031356062, ChipID_19031356068, ChipID_19031356070, ChipID_1903135607, ChipID_19216457008, ChipID_19216457009, ChipID_19216457012, ChipID_19216457014, ChipID_19216457025, ChipID_19216457029, ChipID_19216457032, ChipID_19216457033,ChipID_19234921047, ChipID_19234921059, ChipID_19234921061, ChipID_19234921063, ChipID_19234921065, ChipID_19234921066, ChipID_19234921070, ChipID_19234921074,ChipID_19234921077, ChipID_19234921082, ChipID_19234921083, ChipID_19234921090, ChipID_19234921094, ChipID_19234921100, ChipID_19234921101, ChipID_19235792061,ChipID_19235792082, ChipID_19235792088, ChipID_19235792089, ChipID_19235792091, ChipID_19235792095, ChipID_19249896045, ChipID_19249896067, ChipID_19249896073,ChipID_19249896091, ChipID_19249907011, ChipID_19249907021, ChipID_19249907031, ChipID_19249907044, ChipID_19249907045, ChipID_19249907052, SV1, SV2, SV4, SV5, SV6, SV7))

  imputed_data = data_list$imputed_sets
  personal_cofounders = subset(data_list$warddays_enet, select=c(gender))
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

input_data=read_data("/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct", "ALL_SV6_drugnaive_hba_chip_01-09-16.RData")

data_list = prep_data(
  input_data$imp_list_Tall,
  input_data$GX_sva,
  input_data$svmod_warddays_sv,
  input_data$svmod_referraldays_sv,
  input_data$pheno_mini,
  input_datatpheno
)#this works because order in function definition is maintained ie.covars, sva_file, warddays_file, referral_file
#Set the folds for each BMI

prepared_data = prepare_feature_sets(data_list, outptut_dir = "/users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/new_protocol_25thOct/IMPUTATION_SETS")
names(prepared_data)
