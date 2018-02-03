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

fun_bmi<-function(imputed_single, XL) {
  combination<-cbind(imputed_single, XL)
  bmi.enet.Y<<-combination[,'BMI.norm']  # << double assignment creates resulting object in the global envirnmental. Content of object will change with each 'm' in CV loop
  remaining = combination[,-4] #4:BMI.norm, 5:hba1c.norm, 6:BMI.obese, 7:hba1c.diabetes
  
  print("Removing following features from the dataset:")

  print(paste0("After removal:", ncol(remaining)))
  return(list(target=bmi.enet.Y, features=remaining))
}

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


set_folds_bmi<-function(bmi.enet, f, n){
  set.seed(1)
  theFolds<-sample(rep(seq(f), length.out=nrow(bmi.enet)))#Create f folds
  fold.10000 <- t(replicate(10000, sample(theFolds))) #Create 10FoldCV a 10000 times
  # print(sum(duplicated(fold.10000)))
  #out<-fold.10000[!(duplicated(fold.10000)), ][1:10000] #impt if sum(duplicated(fold.10000)) >0
  fold.10000<-unique(fold.10000);
  folds<-fold.10000[sample(nrow(fold.10000), n), ] #100 iterations of 10 folds
  AllImputedSets<-vector('list', 0)
  #bmi.XL<-subset(bmi.enet, select=-c(hlth, ovwgt, obese, Age.norm, PC1, PC2, hba1c.catgprediabetes, hba1c.catgdiabetes, med.days, GLUCOSE_GRS, OBSESITY_GRS))#removal of variables to be filled by loop
  bmi.XL<-subset(bmi.enet, select=-c(hlth, ovwgt, obese, Age.norm, PC1, PC2, med.days, GLUCOSE_GRS, OBSESITY_GRS)) #removal of variables to be filled by loop
  
  return(list(folds=folds, bmi_XL=bmi.XL))
}

input_data=read_data("/Users/ti1/Google\ Drive/raw\ data/", "imputation_data_generated.RData")

data_list = prep_data(
  input_data$imp_list_Tall,
  input_data$GX_sva,
  input_data$svmod_bmi_catg_sv,
  input_data$svmod_hba1c_catg_sv,
  input_data$pheno_mini,
  input_datat$pheno
)#this works because order in function definition is maintained ie.covars, sva_file, bmi_sv_file, hba_sv_file
#Set the folds for each BMI
fold_and_data <- set_folds_bmi(data_list$bmi_enet, 10, 100)


full_data = fun_bmi(data_list$imputed_sets[[as.numeric(1)]], fold_and_data$bmi_XL)



# apply PCA - scale. = TRUE is highly 
# advisable, but default is FALSE. 
features= full_data$feature[, seq(22, ncol(full_data$feature))]
#features= full_data$feature


ir.pca <- prcomp(t(fold_and_data$bmi_XL),
                 center = TRUE,
                 scale. = TRUE) 

df=data.frame(x=ir.pca$x[,1], y=ir.pca$x[,2], target = full_data$target, type = full_data$feature[,1], chip=input_data$svmod_hba1c_catg_sv[,7])
ggplot(df, aes(x,y, color=type)) + geom_point()


ir.pca <- prcomp(features,
                 center = TRUE,
                 scale. = TRUE) 

df=data.frame(x=ir.pca$x[,1], y=ir.pca$x[,2])
ggplot(df, aes(x,y)) + geom_point()

head(sort(ir.pca$x), 100)
head(sort(ir.pca$x, decreasing = TRUE), 100)


