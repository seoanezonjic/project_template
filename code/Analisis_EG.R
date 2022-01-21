# WGCNA

# INSTALACIÓN DE WGCNA

# install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
# BiocManager::install(c("GO.db", "preprocessCore", "impute"))
# 
# install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival")
# )
# #install.packages("~/Downloads/WGCNA_1.67.tar", repos = NULL, lib=.Library)
# install.packages("~/Downloads/WGCNA_1.67.tar", repos = NULL, lib=.Library, type = "source")
# #library(WGCNA) # Error, probably impute is not installed
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("impute")

# CONFIGURACIÓN DE R Y PROCESAMIENTO DE DATOS

getwd();

workingDir = "C:/Users/Sergio/Desktop/Proyecto BS";
setwd(workingDir); 

library(WGCNA)

sars_Cov2 <- read.table(file = "GSE147507.tsv", header = TRUE, sep = "\t")

dim(sars_Cov2)
names(sars_Cov2)

# REAJUSTE DE DATOS

datExpr0 = as.data.frame(t(sars_Cov2[,-1]));
names(datExpr0) = sars_Cov2$X;
rownames(datExpr0) = names(sars_Cov2)[-1];

# COMPROBACIÓN DE DATOS EN BUSCA DE VALORES FALTANTES

gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

if (!gsg$allOK)
{
  
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

abline(h = 150000, col = "red")


clust = cutreeStatic(sampleTree, cutHeight = 150000, minSize = 10)
table(clust)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# CARGA DE DATOS DE RESISTENCIA

traitData = read.csv("resistanceBC.csv");
dim(traitData)
names(traitData)



samples = rownames(datExpr);
traitRows = match(samples, traitData$sample_ID);
datTraits = traitData[traitRows, ];
rownames(datTraits) = traitData[traitRows, 1];
datTraits$sample_ID = NULL

collectGarbage()


sampleTree2 = hclust(dist(datExpr), method = "average")

traitColors = numbers2colors(datTraits, signed = FALSE);

plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = names(datTraits), 
                    main = "Sample dendrogram and trait heatmap")

# RECONFIGURACIÓN DE R

options(stringsAsFactors = FALSE);

allowWGCNAThreads()
enableWGCNAThreads()

# ANÁLISIS DE LA TOPOLOGÍA DE LA RED


powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");

abline(h=0.90,col="red")

plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

net = blockwiseModules(datExpr, power = powers,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "ecoliTOM", 
                       verbose = 3)

table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]

# RELACIONAR MODULOS CON RESISTENCIAS


nGenes = ncol(datExpr);
nSamples = nrow(datExpr);
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)
moduleTraitCor = cor(MEs, datTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)

sizeGrWindow(10,6)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 0.5,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))


resistance = as.data.frame(datTraits$resistance_BC);
names(resistance) = "resistance"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

geneTraitSignificance = as.data.frame(cor(datExpr, resistance, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(resistance), sep="");
names(GSPvalue) = paste("p.GS.", names(resistance), sep="")

module = "brown"
column = match(module, modNames);
moduleGenes = moduleColors==module;

sizeGrWindow(7, 7);
par(mfrow = c(1,1));
verboseScatterplot(abs(geneModuleMembership[moduleGenes, column]),
                   abs(geneTraitSignificance[moduleGenes, 1]),
                   xlab = paste("Module Membership in", module, "module"),
                   ylab = "Gene significance for resistance",
                   main = paste("Module membership vs. gene significance\n"),
                   cex.main = 1.2, cex.lab = 1.2, cex.axis = 1.2, col = module)

# ANÁLISIS FUNCIONAL

GO_analysis <- function(gene) {
  ego <- enrichGO(gene          = gene,
                  OrgDb         = 'org.EcK12.eg.db',
                  ont           = "ALL",
                  pAdjustMethod = "BH")
  
  head(ego)
  
}


library(clusterProfiler)

data(geneList, package="DOSE")
genesFilt <- filter(GSPvalue, p.GS.resistance <= 0.01)
genes <- rownames(genesFilt)
GO_Ecoli <- GO_analysis(genes)

