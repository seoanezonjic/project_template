# WGCNA

# INSTALACIÓN DE WGCNA

# Automatico

#install.packages("BiocManager")
BiocManager::install("WGCNA")

# Manual (R version 3.0.0 and higher)

# install.packages(c("matrixStats", "Hmisc", "splines", "foreach", "doParallel", "fastcluster", "dynamicTreeCut", "survival", "BiocManager"))
# BiocManager::install(c("GO.db", "preprocessCore", "impute"));

# ENTRADA DE DATOS, LIMPIEZA Y PREPROCESAMIENTO

# CONFIGURACIÓN DE R Y PROCESAMIENTO DE DATOS

getwd();

workingDir = "C:/Users/Sergio/Desktop/Proyecto BS";
setwd(workingDir); 

library(WGCNA)

# Dataset de perfiles de expresión génica
sars_Cov2 <- read.table(file = "GEO-GSE147507.tsv", header = TRUE, sep = "\t")

# Vistazo a los datos
dim(sars_Cov2)
names(sars_Cov2)

# REAJUSTE DE DATOS

# Transformación de los datos: genes <- filas, muestras <- columnas
datExpr0 = as.data.frame(t(sars_Cov2[,-1]));
names(datExpr0) = sars_Cov2$X;
rownames(datExpr0) = names(sars_Cov2)[-1];

# COMPROBACIÓN DE DATOS EN BUSCA DE VALORES FALTANTES

# Busqueda de genes con valores perdidos
gsg = goodSamplesGenes(datExpr0, verbose = 3);
gsg$allOK

# Si la última comprobación devuelve TRUE, no se necesita ejecutar el siguiente código. De lo contrario, se eliminan los genes y muestras prescindibles
if (!gsg$allOK)
{
  
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
  
  datExpr0 = datExpr0[gsg$goodSamples, gsg$goodGenes]
}

# Agrupamos las muestras para ver si hay valores atípicos
sampleTree = hclust(dist(datExpr0), method = "average");

sizeGrWindow(12,9)

par(cex = 0.6);
par(mar = c(0,4,2,0))
pdf(file = 'sampleClustering.pdf', width = 12,height = 12)
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, 
     cex.axis = 1.5, cex.main = 2)

abline(h = 150000, col = "red")
dev.off()

# Detectamos algunos valores atípicos. Los eliminamos eligiendo un corte de altura
clust = cutreeStatic(sampleTree, cutHeight = 150000, minSize = 10)
table(clust)

keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# GUARDADO DE DATOS

# La variable datExpr tiene los datos preparados para el análisis de red
save(datExpr, file = "SARS_Cov2_datExpr.RData")
