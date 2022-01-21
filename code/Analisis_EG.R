# WGCNA

# INSTALACIÓN DE WGCNA

install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
BiocManager::install(c("GO.db", "preprocessCore", "impute"))

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

workingDir = ".";
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
