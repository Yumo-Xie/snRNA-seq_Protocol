.libPaths("~/Rproject/JinHee_snRNA/renv/profiles/seurat5/renv/library/linux-ubuntu-bionic/R-4.4/x86_64-pc-linux-gnu")
options(bitmapType = "cairo")

library(SoupX)
library(Seurat)
library(DoubletFinder)
library(tidyverse)
source("/home/yxie5/Rproject/JinHee_snRNA/Scripts/load10XFlex.R")
source("/home/yxie5/Rproject/JinHee_snRNA/Scripts/run_doubletfinder_custom.R")

# Get the array task ID from command line arguments
args = commandArgs(trailingOnly = TRUE)
i = as.numeric(args[1])

dir1 = "/fh/fast/grady_w/SR/ngs/illumina/jkim6/20251014_LH00740_0142_A235C3FLT3/Analysis/cellranger_9.0.1/GEM1_results/outs/per_sample_outs/"
dir2 = "/fh/fast/grady_w/SR/ngs/illumina/jkim6/20251014_LH00740_0142_A235C3FLT3/Analysis/cellranger_9.0.1/GEM2_results/outs/per_sample_outs/"

sample = append(list.files(dir1), list.files(dir2))

dir1 = file.path(dir1, list.files(dir1))
dir2 = file.path(dir2, list.files(dir2))
dir = append((dir1), (dir2))

epiGenes = c("Muc2", "Muc13", "Muc3", "Epcam", "Krt8")
smcGenes = c("Myh11", "Acta2", "Tagln", "Mylk")

# Process only the i-th sample
sc = load10XFlex(file.path(dir[i], "count"))
useToEst = estimateNonExpressingCells(sc, 
                                      nonExpressedGeneList = list(EPI = epiGenes, SMC = smcGenes),clusters = F,maximumContamination = 0.3)

sc = calculateContaminationFraction(sc, list(EPI = epiGenes, SMC = smcGenes), useToEst = useToEst)

out = adjustCounts(sc)
cntSoggy = rowSums(sc$toc > 0)
cntStrained = rowSums(out > 0)
mostZeroed = tail(sort((cntSoggy - cntStrained)/cntSoggy), n = 10)
print(mostZeroed)

png(paste0("/home/yxie5/Rproject/JinHee_snRNA/Plots/SoupX/", sample[i], ".png"))
plotChangeMap(sc, out, DR = sc$metaData[, c("tSNE1", "tSNE2")], "Muc2")
dev.off()

obj = CreateSeuratObject(out)
obj[["percent.mt"]] <- PercentageFeatureSet(object = obj, pattern = "^mt-")
obj <- subset(obj, subset = nFeature_RNA > 200 & percent.mt < 5 & nCount_RNA > 200)

doublets = run_doubletfinder_custom(obj)
obj = AddMetaData(obj, doublets, col.name = "doublet_finder")

# Save individual sample results
qs::qsave(obj, file = paste0("/home/yxie5/Rproject/JinHee_snRNA/Data/Dec13_obj_", sample[i], "_after_SoupX_doublets.qs"))
