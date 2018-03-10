
library(limma)

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

input_data=read_data("/Users/ti1/Google\ Drive/raw\ data", "imputation_data_generated.RData")


fat = (rowSums(input_data$svmod_bmi_catg_sv[,2:3]))
new_cat = input_data$svmod_bmi_catg_sv[,-3]
new_cat[,2]=fat
new_cat = new_cat[, 1:2]

gene_names = read.table("/Users/ti1/Google Drive/validation_prediction/gene_names.csv" , stringsAsFactors = FALSE)
nuID = gene_names[,1]
symbols = gene_names[,2]
gene_expression_train = t(input_data$GX_sva)
gene_expression_train = gene_expression_train[,which(colnames(gene_expression_train) %in% nuID)]
genes_train = symbols[which(nuID  %in%  colnames(gene_expression_train))]
colnames(gene_expression_train) = genes_train


fit1.bmi=lmFit(t(gene_expression_train), new_cat)
#fit1.hba1c=lmFit(input_data$GX_sva, input_data$svmod_hba1c_catg)



#Create contrast matix
#cols = colnames(input_data$svmod_bmi_catg_s) 
#cols[0:3]= c("norm", "pre", "dbx")
#colnames(input_data$svmod_bmi_catg_s) = cols
contrast.matx.bmi <- makeContrasts (ovwgt-hlth, levels =  new_cat) 

fit2.bmi=contrasts.fit(fit1.bmi, contrast.matx.bmi)
fit2.hba=contrasts.fit(fit1.hba1c, contrast.matx.hba)

fit2.bmi<-eBayes(fit2.bmi) ###use eBayes() to moderate the estimated error variances
#fit2.hba<-eBayes(fit2.hba)

#bmi: "hlth", "hlth.fu", "obese", "ovwgt", "ovwtplus.fu"
#BMI.obese = topTable(fit2.bmi, coef = "ovwgt - hlth", sort.by="p") #obese-healthy
best = topTable(fit2.bmi,adjust.method="fdr", n=1000, sort.by="P")

best = (best[best[,5] < 0.01,])
genes = best[,1]
write.csv(symbols[!symbols %in% genes], "/Users/ti1/Google\ Drive/validation_prediction/remove_non_diff_genes.txt", row.names = FALSE, quote = FALSE, col.names = NULL)
#on merg