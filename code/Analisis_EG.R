# WGCNA

# INSTALACIÓN DE WGCNA

# Automatico

install.packages("BiocManager")
BiocManager::install("WGCNA")

# Manual (R version 3.0.0 and higher)

# install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
# BiocManager::install(c("GO.db", "preprocessCore", "impute"));

# ENTRADA DE DATOS, LIMPIEZA Y PREPROCESAMIENTO

# CONFIGURACIÓN DE R Y PROCESAMIENTO DE DATOS

getwd();

workingDir = ".";
setwd(workingDir); 

library(WGCNA)

# Dataset de perfiles de expresión génica
sars_Cov2 <- read.table(file = "GSE147507.tsv", header = TRUE, sep = "\t")

# Vistazo a los datos
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

# GUARDADO DE DATOS

save(datExpr, file = "SARS_Cov2_datExpr.RData")

