# WGCNA II

# CONFIGURACIÓN DE R

getwd();

workingDir = ".";
setwd(workingDir); 

library(WGCNA)

options(stringsAsFactors = FALSE);

allowWGCNAThreads() 
enableWGCNAThreads()

lnames = load(file = "SARS_Cov2_datExpr.RData");

lnames

# ANÁLISIS DE LA TOPOLOGÍA DE LA RED

# Selección de umbrales
powers = c(c(1:10), seq(from = 12, to=20, by=2))

sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;

pdf(file = 'independenceScale_meanConnectivity.pdf', width = 8,height = 6)
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
dev.off()

# Transformación del tipo de valores para todas las columnas (Una función no detecta los valores si no son tipo NUMERIC)
test <- datExpr

test <- setNames(data.frame(lapply(test, as.numeric)), 
                 colnames(test))
sum(sapply(test, is.numeric))

# Construcción de la red de genes
net = blockwiseModules(test, power = powers,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "SARS_Cov2TOM", 
                       verbose = 3)

# Identificación de módulos
table(net$colors)

sizeGrWindow(12, 9)
mergedColors = labels2colors(net$colors)

pdf(file = 'clusterDendrogram.pdf', width = 8,height = 6)
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

# Información de genes y módulos
moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]]
