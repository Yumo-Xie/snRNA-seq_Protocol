library(SoupX)
library(Seurat)
library(DoubletFinder)
library(tidyverse)
source("/home/yxie5/Rproject/JinHee_snRNA/Scripts/load10XFlex.R")
source("/home/yxie5/Rproject/JinHee_snRNA/Scripts/run_doubletfinder_custom.R")

# Get the array task ID from command line arguments
args = commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])

dir = "dir_to_files"

epiGenes = c("Muc2", "Muc13", "Muc3", "Epcam", "Krt8")
smcGenes = c("Myh11", "Acta2", "Tagln", "Mylk")
immuGenes = c("Igkc","Igha","Ighm")


sc = load10XFlex(file.path(dir[i], "count"))
useToEst = estimateNonExpressingCells(sc, 
                                      nonExpressedGeneList = list(EPI = epiGenes, SMC = smcGenes, Immu = immuGenes),clusters = F,maximumContamination = 0.3)

sc = calculateContaminationFraction(sc, list(EPI = epiGenes, SMC = smcGenes, Immu = immuGenes), useToEst = useToEst)

out = adjustCounts(sc)
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
print(mostZeroed)

png(paste0(".../Plots/SoupX/", sample[i], ".png"))
plotChangeMap(sc, out, DR = sc$metaData[, c("tSNE1", "tSNE2")], "Muc2")
dev.off()

obj = CreateSeuratObject(out)
obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^mt-")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 5 & nCount_RNA > 200)

doublets = run_doubletfinder_custom(obj)
obj = AddMetaData(obj, doublets, col.name = "doublet_finder")

# Save
qs::qsave(obj, file = paste0("..../obj_", sample[i], "_after_SoupX_doublets.qs"))
