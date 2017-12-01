library(topGO)
library(biomaRt)


genes =read.table('/Users/ti1/Google\ Drive/results/default_accuracy_analysis/all/remove_top__random_100.csv', sep=",", stringsAsFactors = FALSE)
genes[,2]

length(gene_symbols$REFSEQ_ID)

genes_in = gene_symbols$REFSEQ_ID[which(gene_symbols$nuID  %in%  genes[,2])]
all_genes <- gene_symbols$REFSEQ_ID

genes_in=gsub("\\..*", "", genes_in)
all_genes=gsub("\\..*", "", all_genes)


ensembl = useEnsembl(biomart="ensembl", dataset="hsapiens_gene_ensembl")
genes_in = getBM(attributes=c("ensembl_gene_id", "refseq_mrna"), filters="refseq_mrna",
                 values= genes_in,
                 mart=ensembl)

all_genes = getBM(attributes=c("ensembl_gene_id", "refseq_mrna"), filters="refseq_mrna",
                 values= all_genes,
                 mart=ensembl)

geneList = all_genes[all_genes[,1] %in% genes_in[,1], ]

geneList <- factor(as.integer (all_genes[,1] %in% geneList[,1]))
names (geneList) <- all_genes[,1]

# then make a factor that is 1 if the probeset is "interesting" and 0 otherwise
# name the factor with the probeset names
# form the GOdata object

GOdata <-new ("topGOdata", 
              ontology = "BP", 
              allGenes = geneList, 
              annot = annFUN.org, mapping = "org.Hs.eg.db", 
              ID = "Ensembl"
)

resultFisher <- runTest(GOdata, algorithm = "classic", statistic = "fisher")
resultKS <- runTest(GOdata, algorithm = "classic", statistic = "ks")
resultKS.elim <- runTest(GOdata, algorithm = "elim", statistic = "ks")

allRes <- GenTable(GOdata, classicFisher = resultFisher,classicKS = resultKS, elimKS = resultKS.elim, orderBy = "classicFisher", ranksOf = "classicFisher", topNodes = 10)
showSigOfNodes(GOdata, score(resultFisher), firstSigNodes = 5, useInfo = 'all')
