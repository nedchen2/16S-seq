#install.packages("BiocManager")
#BiocManager::install("WGCNA")

require(WGCNA)
require(qiime2R)
require(tidyverse)
metadata <- 
data <-read_qza("../results/4.Diversity_ana/core-metrics-results/rarefied_table.qza")$data %>% t()

data <- read.csv("../results/5.Enterotype/genus_abundance_data.csv",check.names = F , row.names = 1) %>% t()

library("vegan")
HellingerData<- decostand(data,method = "hellinger") # Hellinger Transformation (a square root of the relative abundance)

OTUs <- HellingerData %>% as.data.frame()
names(OTUs)

goodSamplesGenes(OTUs, verbose = 3)



#rename 
datExpr0 = OTUs

# =================== detect outliers
sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))

plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)

#  ===================== import the parameter trait

traitData = read.csv("../code/metadata_Microbiome.csv") %>% dplyr::select("SampleID","Nosema_total_per_bee"  ,"Crithidia_total_per_bee", "Apicystis_total_per_bee") 
dim(traitData)
names(traitData)

OTUSamples = rownames(datExpr0);
traitRows = match(OTUSamples, traitData$SampleID);
datTraits = log10(traitData[traitRows, -1]+1);
rownames(datTraits) = traitData[traitRows, 1];
collectGarbage()


# = plot dendrogram


sampleTree2 = hclust(dist(datExpr0), method = "average")

  
traitColors = numbers2colors(datTraits[1:3], signed = FALSE);
  
plotDendroAndColors(sampleTree2, traitColors,
                      groupLabels = names(datTraits[1:3]),
                      main = "Sample dendrogram and trait heatmap")


# =====
options(stringsAsFactors = FALSE);
enableWGCNAThreads()

powers = c(c(1:10), seq(from = 11, to=30, by=1))

sft = pickSoftThreshold(datExpr0, powerVector = powers, verbose = 5, networkType = "signed")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
    xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
    main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers,col="red")

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,col="red")

# ====== power choose to 10

softPower = 10;
adjacency = adjacency(datExpr0, power = softPower, type = "signed")

TOM = TOMsimilarity(adjacency, TOMType = "signed");
dissTOM = 1-TOM

TaxaTree = hclust(as.dist(dissTOM), method = "average")

plot(TaxaTree, xlab="", sub="", main = "Taxa clustering on TOM-based dissimilarity",
     labels = FALSE, hang = 0.04)

minModuleSize = 5

dynamicMods = cutreeDynamic(dendro = TaxaTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);

table(dynamicMods)
  
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)

plotDendroAndColors(TaxaTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Taxa dendrogram and module colors")

#Calculate eigengenes:
  
MEList = moduleEigengenes(datExpr0, colors = dynamicColors)
MEs = MEList$eigengenes
#Calculate dissimilarity of module eigengenes:
  
MEDiss = 1-cor(MEs);
#Cluster module eigengenes:
  
METree = hclust(as.dist(MEDiss), method = "average");
#Plot the result:
  
plot(METree, main = "Clustering of module eigengenes",
     xlab = "", sub = "")



MEDissThres = 0.40
#Plot the cut line into the dendrogram:
  
# --- no module need to be merged

merge = mergeCloseModules(datExpr0, dynamicColors, cutHeight = MEDissThres, verbose = 3)
#The merged module colors:
  
mergedColors = merge$colors;
#Eigengenes of the new merged modules:
  
mergedMEs = merge$newMEs;

moduleColors = mergedColors
  
colorOrder = c("grey", standardColors(50));
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs


save(MEs, moduleLabels, moduleColors, TaxaTree, file = "Monterey-networkConstruction-stepByStep.RData")

# ================== code for 

#Defining numbers of OTUs and samples:
nTaxa = ncol(datExpr0);
nSamples = nrow(datExpr0);
#Recalculate MEs (module eigengenes):
  
MEs0 = moduleEigengenes(datExpr0, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p");
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);
#Now we will visualize it:
  

textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)

pdf("../results/7.Final_graph/Heatmap_module_WGCNA.pdf")

par(mar = c(6, 10, 8, 8));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))

dev.off()

# Each row corresponds to a module eigengene and each column corresponds to an environmental trait or biogeochemical rate (as long as it is continuous–notice that the categorical variables are gray and say NA). Each cell contains the corresponding Pearson correlation coefficient (top number) and a p-value (in parentheses). The table is color-coded by correlation according to the color legend.
# You can see that the Brown module is positively correlated with many indices of upwelling while the Black module is negatively correlated with many indices of upwelling. 
# For this work I was particularly interested in CR and so I focused on modules the positively or negatively correlated with CR. 
# The Red module was negatively associated with CR while the Blue module was positively associated with CR



Nosema = as.data.frame(datTraits$Nosema_total_per_bee);
names(Nosema) = "Nosema"

Nosema = as.data.frame(datTraits$Crithidia_total_per_bee);
names(Nosema) = "Nosema"

modNames = substring(names(MEs), 3)
TaxaModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
names(TaxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
TaxaTraitSignificance = as.data.frame(cor(datExpr0, Nosema, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
names(TaxaTraitSignificance) = paste("GS.", names(Nosema), sep="");
names(GSPvalue) = paste("p.GS.", names(Nosema), sep="");

module = "yellow"
column = match(module, modNames);
moduleTaxa = moduleColors==module;
sizeGrWindow(7, 7)
par(mfrow = c(1,1))

verboseScatterplot(abs(TaxaModuleMembership[moduleTaxa, column]),
                   abs(TaxaTraitSignificance[moduleTaxa, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Taxa significance for Nosema",
                   main = paste("Module membership vs. Taxa significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)


# =========== Crithidia 
Crithidia = as.data.frame(datTraits$Crithidia_total_per_bee);
names(Crithidia ) = "Crithidia "

modNames = substring(names(MEs), 3)
TaxaModuleMembership = as.data.frame(cor(datExpr0, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaModuleMembership), nSamples));
names(TaxaModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");
TaxaTraitSignificance = as.data.frame(cor(datExpr0, Crithidia , use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(TaxaTraitSignificance), nSamples));
names(TaxaTraitSignificance) = paste("GS.", names(Crithidia ), sep="");
names(GSPvalue) = paste("p.GS.", names(Crithidia ), sep="");

module = "grey"
column = match(module, modNames);
moduleTaxa = moduleColors==module;

pdf("../results/7.Final_graph/gray_module_scatter.pdf")
verboseScatterplot(abs(TaxaModuleMembership[, column]),
                   abs(TaxaTraitSignificance[moduleTaxa, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Taxa significance for Crithidia",
                   main = paste("Module membership vs. Taxa significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = "grey7", showPValue = F,lmFnc=T)

dev.off()


# This graph shows you how each taxa (each red dot is an OTU that belongs in the Red module) correlated with 
# 1) the Environmental trait of interest and 2) how important it is to the module. The taxa/OTUs that have high module membership tend to occur 
# whenever the module is represented in the environment and are therefore often connected throughout the samples with other red taxa/OTUs. 
# In this module, these hubs (Red OTUs that occur with other Red OTUs) are also the most important OTUs for predicting parasite loads.


names(datExpr0)
names(datExpr0)[moduleColors=="grey"]

names(datExpr0)[moduleColors=="brown"]

names(datExpr0)[moduleColors=="yellow"]

names(datExpr0)[moduleColors=="purple"]


# ================ for OTU BASED ANNOTATION
annot = read_qza("../results/4.Diversity_ana/taxonomy-corrected.qza")$data %>% dplyr::rename( OTU = Feature.ID )
dim(annot)
names(annot)
probes = names(datExpr0)
probes2annot = match(probes, annot$OTU)


TaxaInfo0 = data.frame(Taxon = probes,
                       TaxaSymbol = annot$OTU[probes2annot],
                       LinkID = annot$Taxon[probes2annot],
                       moduleColor = moduleColors,
                       TaxaTraitSignificance,
                       GSPvalue)

modOrder = order(-abs(cor(MEs, CR, use = "p")));
#Add module membership information in the chosen order:
  
for (mod in 1:ncol(TaxaModuleMembership)){
    oldNames = names(TaxaInfo0)
    TaxaInfo0 = data.frame(TaxaInfo0, TaxaModuleMembership[, modOrder[mod]],
                           MMPvalue[, modOrder[mod]]);
    names(TaxaInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                         paste("p.MM.", modNames[modOrder[mod]], sep=""))
  }
#Order the OTUs in the geneInfo variable first by module color, then by geneTraitSignificance:
  
TaxaOrder = order(TaxaInfo0$moduleColor);
TaxaInfo = TaxaInfo0[TaxaOrder, ]

write.csv(TaxaInfo, file = "../results/6.Network_analysis/WGCNA/TaxaInfo.csv")
