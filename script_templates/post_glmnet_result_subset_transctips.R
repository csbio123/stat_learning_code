#P value analysis conrad
pval_anal<-function(path, top=0, type, all_genes_names){

		filenames <- list.files(path, pattern = "pvalues*", full.names = TRUE)
		files = sapply(filenames, read.csv, sep="\t", simplify = FALSE)

		p_value_results  = do.call(rbind, files)
		positive_freq = table(p_value_results[c("Genes","positive_binary")])
		positive_freq = cbind(positive_freq, positive_freq[,2]/length(files))
		genes=levels(droplevels(p_value_results[,1]))
		probs = sapply(unique(genes), function(x) { 

				sub = p_value_results[which(genes == x), ]

				if(nrow(sub) > 1) {
				p = mean(sub["positive_p_value"])
				} else {
				p = sub["positive_p_value"]
				}
				return (p)
				}, simplify=FALSE)

	probs = do.call(rbind, probs)
		probs = cbind(probs, rownames(probs))
		positive_freq[order(positive_freq[,1], positive_freq[,3], decreasing = TRUE),]
		top_genes = rownames(positive_freq)
		like = sapply(0:length(top_genes),function(x){
				p = probs[top_genes[x] == probs[,2],]
				p = as.numeric(p[1])
				k = as.numeric(positive_freq[x,2])
				likelihood = dbinom(k, as.numeric(length(files)), prob = p)
				return(likelihood)
				}, simplify = FALSE)
	like = do.call(rbind, like)

		positive_freq = cbind(positive_freq, like)
		positive_freq = cbind(abs(positive_freq[,2] -  positive_freq[,1]), positive_freq)
		positive_freq = positive_freq[order(positive_freq[,1], decreasing = TRUE),]

		if (top >0) {
		  

		  positive_genes = intersect(as.character(all_genes_names), as.character(rownames(positive_freq)))
		  
		  subset_genes = head(positive_genes, as.numeric(top))
		  print(subset_genes)
		  out_name = '/remove_top_'
		  
		  if (type == "random") {
		    print("random stuff")
		    subset_genes = sample(positive_genes, as.numeric(top))
		    out_name = paste0(out_name, "_random_")
		    print(subset_genes)
		  }
		  
			write.table(file=paste(path, out_name, top , '.csv', sep=""), subset_genes, quote = FALSE, row.names = FALSE, col.names = FALSE)
		}
}

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


#arg function below is in-built. It allows you to interact with the script from the command line.

#####USE AS FOLLOWS#####: 
#Rscript p_value_analysis.R /users/spjtcoi/brc_scratch/project_tomi/conrad/reanalyse/drug_naive/all_transcripts 476####

args <- commandArgs(trailingOnly = TRUE)#trailing only stops the argument function from requiring specification of too much information eg R version, etc
print(paste0('Input directory for p value analysis:', args[1]))#This line will tell you the input directory
args[1] = "/Users/ti1/Documents/conrad/results/default_data/all_features/"
args[2] = 10
args[3] = "ranasddom"

input_data=read_data("/Users/ti1/Google\ Drive/raw data/", "imputation_data_generated.RData")
my_genes = rownames(input_data$GX_sva)

if (length(args) == 3) {
	pval_anal(path = args[1], top = args[2], type = args[3], all_genes_names=my_genes)#if 2 arguments are entered then the first will be the path and the second will be the number of top genes to remove
}

