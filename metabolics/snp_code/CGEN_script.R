
#Recoding SNP genotype data with minor allele as the reference
##############################################################
#1) CGEN package
install.packages(CGEN)#1st time only
library(CGEN)#each session

#2) Copy SNP IDs
snps = features[31:82]
SNPIDs <- names(snps)

allells = read.table("/Users/ti1/Google\ Drive/study2_snps_as\ predictors/SNPNAMES_AND_ALLELES_FOR_RECODE.txt", header=1)


#3)Transformation script
recode_raw <- lapply(SNPIDs, function(x) {
  recode(snps[,x], allells[[x]], values=c(0,1,2))
})#Returns list of SNPs each with genotypes and alleles as sublist

recode_slim<-data.frame(lapply(recode_raw, function(x) x[1]))#remove allele lists
colnames(recode_slim)<-SNPIDs

recoded_Final<-recode_slim
data_merge.recodes<-cbind(data_merge, recoded_Final) #combines #original dataset with the newly-recoded data. ID ordering should be #consistent. This can be demonstrated in the next step..

#End of section housekeeping..
rm(recode_raw, recode_slim, data_merge, recoded_Final, SNPIDs)


CHECK VALIDITY OF RECODED DATA
#####################################################

#1)VISUAL CHECK
CompareSNPs<-data.frame(cbind(data_merge.recodes[3:54], data_merge.recodes[4838:4889]))
CompareSNPs.bycolumn<-CompareSNPs[order(colnames(CompareSNPs))] 

#2)GENOTYPIC CHECK
#Compare the (uncoded and recoded) genotype ratios for given SNP 
sort(table(CompareSNPs.bycolumn$rs823165), decreasing = TRUE)
sort(table(CompareSNPs.bycolumn$rs823165.1), decreasing = TRUE)

#3)MANUAL CHECK
#Carry out manual recodes and check for consistency with automated recodes
library(stringr)
CompareSNPs.bycolumn$rs823165<-as.factor(str_replace_all(CompareSNPs.bycolumn$rs823165, c("AA" = 0, "TA" = 1, "TT" = 2)))
#NB. alleles and #genotypes will vary for each SNP

identical(CompareSNPs.bycolumn$rs823165, CompareSNPs.bycolumn$rs823165.1)
