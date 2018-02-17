#  Load/Prepare data
##=====================================================================================
rm(list=ls())
library(WGCNA);

setwd("C:/Users/spjtcoi/Google Drive/WGCNA")
load("C:/Users/spjtcoi/Google Drive/raw data/ALL_SV6_drugnaive_hba_chip_2factor_01-09-17.RData")

t.GX.sva<-t(GX.sva)



#check for genes and samples with too many missing values:
#=====================================================================================
good_samples_genes = goodSamplesGenes(t.GX.sva, verbose = 3);
good_samples_genes$allOK
#If the last statement returns TRUE, all genes have passed the cuts. (If not, remove the offending genes and samples
#from the data)


# Optionally, print the gene and sample names that were removed:
{
if (sum(!good_samples_genes$goodGenes)>0) 
  printFlush(paste("Removing genes:", paste(names(t.GX.sva)[!good_samples_genes$goodGenes], collapse = ", ")));
if (sum(!good_samples_genes$goodSamples)>0) 
  printFlush(paste("Removing samples:", paste(rownames(t.GX.sva)[!good_samples_genes$goodSamples], collapse = ", ")));
# Remove the offending genes and samples from the data:
t.GX.sva = t.GX.sva[good_samples_genes$goodSamples, good_samples_genes$goodGenes]
}


#cluster the samples (in contrast to clustering genes that will come later) to see if there are any obvious outliers
#===================================================================================================================

sampleTree = hclust(dist(t.GX.sva), method = "average");
# Plot the sample tree: Open a graphic output window of size 12 by 9 inches
# The user should change the dimensions if the window is too large or too small.
sizeGrWindow(12,9)# Omit
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect FEP outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 45, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 45, minSize = 90) #minSize minimum number of object on a branch to be considered a cluste
table(clust) #Read horizontally. Top row identifies clusters 0 and 1. Second row tells you their sizes
# clust 0 contains the samples we want to keep.
keepSamples = (clust==0)
t.GX.sva2 = t.GX.sva[keepSamples, ]
nGenes = ncol(t.GX.sva2)
nSamples = nrow(t.GX.sva2)
#The variable t.GX.sva2 now contains the expression data ready for network analysis.


#Read in the trait data and match the samples for which they were measured to the expression samples.
#=====================================================================================
#pheno.sva contains the complement phenotypes for GX.sva. Should already exist within the R workspace
#load("pheno_sva-gx_postjoin.RData") if this is not the case

head(data.frame(colnames(pheno.sva)), 54)#GX starts at col43
pheno_sva2<-pheno.sva[-c(1:2, 5:8, 11:13, 18:19, 22, 24, 25, 26:42)]
pheno_sva3 = model.matrix(~ 0 + Timepoint + gender + PC1 + PC2 + med.days  + Age.norm + BMI.norm  + BMI.catg + hba1c.norm + hba1c.catg + TriG.norm, data=pheno_sva2)
pheno_sva3<-data.frame(pheno_sva3[,-8])#remove dummy variable for non-obese subjects
rownames(pheno_sva3)=rownames(t.GX.sva2)
identical(rownames(pheno_sva3), rownames(t.GX.sva2))


#visualize how clinical traits #relate to the sample dendrogram.
#=====================================================================================

# Re-cluster samples
sampleTree = hclust(dist(t.GX.sva2), method = "average")
# Convert traits to a color representation: white means low, red means high, grey means missing entry
traitColors = numbers2colors(pheno_sva3, signed = FALSE);
# Plot the sample dendrogram and the colors underneath.
plotDendroAndColors(sampleTree, traitColors,
                    groupLabels = names(pheno_sva3), 
                    main = "FEP Sample dendrogram and trait heatmap")
#Output:Clustering dendrogram of samples based on Euclidean distance.


#Constructing a weighted gene network 
#=====================================================================================
#entails choosing of 'soft thresholding power' to which co-expression #similarity is raised to calculate adjacency. 
#The pickSoftThreshold function performs the analysis of network topology. User then chooses a set of candidate powers (the function provides suitable default values), 
#and a set of network indices is returned for inspection as follows:

# Choose a set of soft-thresholding powers
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(t.GX.sva2, powerVector = powers, verbose = 5)
# Plot the results:
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
# Scale-free topology fit index as a function of the soft-thresholding power
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="FEP Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
# this line corresponds to using an R^2 cut-off of h
abline(h=0.90,col="red")
# Mean connectivity as a function of the soft-thresholding power
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#We choose the power 5, which is the lowest power for which the scale-free topology the
#index curve flattens out upon reaching a high value (in this case, roughly 0.90).



# Construct the gene network and identify modules (TOM: Topological overlap Matrix)
#=====================================================================================


net = blockwiseModules(t.GX.sva2, power = 5,
                           TOMType = "unsigned", minModuleSize = 30,
                           reassignThreshold = 0, mergeCutHeight = 0.25,
                           numericLabels = TRUE, pamRespectsDendro = FALSE,
                           randomSeed = 12345,
                           saveTOMs = TRUE,
                           saveTOMFileBase = "FirstEpisodePsychosisTOM",
                           verbose = 3)
#mergeCutHeight = 0.25 corresponds to a correlation of .75. Modules reaching this level of pairwise correlation will be merged

table(net$colors)
#there are 16 modules, labeled 1 through 16 in order of descending size, with sizes ranging from
#1839 to 39 genes. The label 0 is reserved for genes outside of all modules.



#Display hierarchical clustering dendrogram (tree) together with the color assignment 
#=====================================================================================

#open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)


#To change aspects of tree cut, module membership, and module merging criteria, use the  recutBlockwiseTrees function that applies the
#modified criteria without having to recompute the network and the clustering dendrogram -Save a substantial amount of time.


#Keep the module assignment and module eigengene information for subsequent analyses.
#=====================================================================================

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]



#Identify modules significantly associated with the clinical obeseity trait .
#=====================================================================================

#Strategy: Correlate eigengene (ie. the first principal component) of each module with
#trait and look for the most significant associations:

# Define numbers of genes and samples
nGenes = ncol(t.GX.sva2);
nSamples = nrow(t.GX.sva2);
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(t.GX.sva2, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, pheno_sva3, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

#We color code each association by the correlation value:  
sizeGrWindow(10,6)
# Will display correlations and their p-values
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                       signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
# Display the correlation values within a heatmap plot
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(pheno_sva3),#names changed to colnames
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("First-episode psychosis Module-trait relationships"))


#Gene relationship to trait and important modules: Gene Significance and Module Membership
#=====================================================================================

#We quantify trait associations of individual genes by defining Gene Significance GS as
#(the absolute value of) the correlation between the gene and the trait. 
#For each module, we also define a quantitative #measure of module membership MM as the correlation of the module eigengene with each transcript. 
#This allows us to quantify the similarity of all genes on the array to every module.

#Define variable weight containing the weight column of datTrait
obese = data.frame(pheno_sva3$BMI.catgobese); #the index phenotype for this dataset
names(obese) = "obese"
# names (colors) of the modules
modNames = substring(names(MEs), 3)
geneModuleMembership = as.data.frame(cor(t.GX.sva2, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));
names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
geneTraitSignificance = as.data.frame(cor(t.GX.sva2, obese, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));
names(geneTraitSignificance) = paste("GS.", names(obese), sep="");
names(GSPvalue) = paste("p.GS.", names(obese), sep="");


#Intramodular analysis: identifying genes with high GS and MM
#=====================================================================================

#Using GS and MM measures, identify genes that have a high significance for weight as well as high module
#membership in interesting modules. 

#1
#As an example, we look at the brown module that has the highest association
#with weight. We plot a scatterplot of Gene Signicance vs. Module Membership in the brown module:
module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for body weight",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#2
module = "turquoise"
#"purple", "turquoise", "lightcyan", "blue", "brown", "cyan", "yellow", "black", "magenta", 
#"midnightblue", "salmon", "green", "pink", "tan", "greenyellow", "red", "grey"        
column = match(module, modNames);
moduleGenes = moduleColors==module;
sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("FEP: Module Membership in", module, "module"),
                   ylab = "Gene significance for obesity",
                   main = paste("FEP: Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

#Note: Generally correlations between module and trait are weak. However in the prediction aspect we are testing the extent to which multiple weakly correlating genes (modules)
#can accurately predict the trait.


#We now merge this statistical information with gene annotation
#=====================================================================================
#names(datExpr)[moduleColors=="brown"] #This will return probe IDs belonging to the desired module colour. 


#connect probe IDs to gene names and other info
load("C:/Users/spjtcoi/Google Drive/gene_symbols.RData") #load gene_symbol object
probes = colnames(t.GX.sva2)
probes2annot = match(probes, gene_symbols$nuID)
# The following is the number or probes without annotation:
sum(is.na(probes2annot))
# Should return 0


#combine probe ID, gene symbol, Locus Link ID (Entrez code), module color, 
#gene significance for weight, and module membership and p-values in all modules
#=====================================================================================

#Mes0/MEs: The eigengene of a module (ie the first principal component)
#MM.red/MM.blue: Module Membership (ie. correlation between a transcript and the eigengene of a module)
#GS: Gene Significance: The empirical correlation between the gene and trait



#The modules will be ordered by their significance for obessity, with THE MOST IMPORTANT ONES TO THE LEFT.

# Create the starting data frame
geneInfo0 = data.frame(nuID = probes,
                           geneSymbol = gene_symbols$SYMBOL[probes2annot],
                           definition = gene_symbols$DEFINITION[probes2annot],
                           moduleColor = moduleColors,
                           geneTraitSignificance,
                           GSPvalue)
# Order modules by their significance for weight
modOrder = order(-abs(cor(MEs, obese, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]],
                             MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                           paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.obese));
geneInfo = geneInfo0[geneOrder, ]


write.csv(geneInfo, file = "geneInfo.csv")
