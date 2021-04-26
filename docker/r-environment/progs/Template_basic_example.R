library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(reticulate)

args <- commandArgs(trailingOnly = TRUE)

# Input dirs
input_dir <- args[1] #
output_dir <- args[2] #
sample_id <- as.character(args[3]) # sample id

seurattools_dir <- args[4] # Dir with seurattools repo
source(file = paste(seurattools_dir, "Basic/DataLoading_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "Basic/PlotQC_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_violin_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "Basic/Normalize_FindHVG.v3.R", sep = "/"))
source(file = paste(seurattools_dir, "Basic/CellScoringMeta_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_plotViz_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_plotPCA_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "Basic/RunDR_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_dimHeatmap_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_plotTSNE_v3_multi.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_FindAllMarkers_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "scTools/scTool_addmodscores_v3.R", sep = "/"))
source(file = paste(seurattools_dir, "ExtraTools/CreateMultiTabExcels_v3.R", sep = "/"))

wd <- paste(output_dir, sample_id, sep = "/")
dir.create(wd, recursive=TRUE, showWarnings = FALSE)
setwd(wd)

# Analysis params
max.pc.mito <- as.numeric(args[5]) # max percent.mt
min.feature <- as.numeric(args[6]) # min nFeature_RNA
max.feature <- as.numeric(args[7]) # max nFeature_RNA
num_pc <- as.numeric(args[8]) # npcs
cc_rds <- args[9] # CC rds file

use_python(args[10]) # reticulate

mito.prefix = args[11]

# Import data
print(paste(input_dir, sample_id, "outs", "filtered_feature_bc_matrix", sep = "/"))
obj <- DataLoading.v3(fold10x = input_dir, 
                      projname = "Single-Cell", id = gsub("_", "-", sample_id), 
                      mincells = 10, minfeatures = 50, mito.pattern = mito.prefix)

saveRDS(object = obj, file = paste0(sample_id, "_minimal.rds"))

# QC 
PlotQC.v3(seurat.obj = obj, rawQC = TRUE)
scTool_violin_v3(seurat.obj = obj, featuresvect = c("nFeature_RNA", "nCount_RNA"))
scTool_violin_v3(seurat.obj = obj, featuresvect = c("percent.mt", "percent.rp"))

# Filtering
obj <- subset(x = obj, subset = nFeature_RNA > min.feature & nFeature_RNA < max.feature
              & percent.mt < max.pc.mito & percent.rp < 100)

# QC post filtering
PlotQC.v3(seurat.obj = obj, rawQC = FALSE)

# OPTION 1 LOG-NORMALIZATION + SCALING
# Normalize and find HVGs
obj <- NormalizeFindHVG.v3(seurat.obj = obj, top.feat2plot = 50, 
                           n.features = 4000, selmethod = "vst") # find variable genes is embedded here

# cell-cycle scoring
obj <- CellScoringMeta.v3(cellcycle.rds = cc_rds, #"/home/sberetta/regev_lab_cell_cycle_mouse.rds", 
                          seurat.obj = obj, ncolumns = 2)
obj <- ScaleData(object = obj, vars.to.regress = c("CC.Difference", "percent.mt", "nCount_RNA"), 
                 display.progress = T, features = rownames(obj))

obj <- RunPCA(object = obj, features = VariableFeatures(object = obj), npcs = 50)
obj <- RunPCA(object = obj, features = rownames(obj), npcs = 50, reduction.name = "pcall", reduction.key = "pcall_")
obj <- RunPCA(object = obj, features = obj@misc$cc.genes, npcs = 50, reduction.name = "pcacc", reduction.key = "pcacc_")

scTool_PlotViz.v3(seurat.obj = obj, viz.feat = 50, start.pc = 1, end.pc = 2)
scTool_PlotViz.v3(seurat.obj = obj, viz.feat = 50, start.pc = 3, end.pc = 4)
scTool_PlotPCA.v3(seurat.obj = obj, start.pc = 1, end.pc = 2, pca.identity = "Phase")
scTool_PlotPCA.v3(seurat.obj = obj, start.pc = 1, end.pc = 2, pca.identity = "Phase", reduction = "pcall")
scTool_PlotPCA.v3(seurat.obj = obj, start.pc = 1, end.pc = 2, pca.identity = "Phase", reduction = "pcacc")
scTool_PlotPCA.v3(seurat.obj = obj, start.pc = 1, end.pc = 2, pca.identity = "Phase", pca.splitby = "Phase")

# DimHeatmap
scTool_DimHeatmap.v3(seurat.obj = obj, dimheat.n.cells = 500, dimheat.end = 9, dimheat.feat = 20)

# JackStraw
obj <- JackStraw(obj, num.replicate = 100, dims = 50)
obj <- ScoreJackStraw(obj, dims = 1:50)

png(filename = paste0("QC/",obj@meta.data$orig.ident[1],"_Jackstrawplot.png"), units = "in", width = 12, height = 9, res=100)
print(JackStrawPlot(obj, dims = 1:50, reduction = "pca") +
        ggtitle(label = paste0("\nJackStraw plot 50 PCs \n"), subtitle = paste0("Sample ",obj@meta.data$orig.ident[1])) + 
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
        theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
        theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
        theme(legend.text = element_text(size=8, hjust=0.5, color="black")) +
        theme(legend.title = element_blank()) +
        theme(axis.text = element_text(size=12, hjust=0.5, color="black")))
dev.off()

png(filename = paste0("QC/",obj@meta.data$orig.ident[1],"_Elbowplot.png"), units = "in", width = 12, height = 9, res=100)
print(ElbowPlot(obj, ndims = 50) +
        ggtitle(label = paste0("\nElbow plot 50 PCs \n"), subtitle = paste0("Sample ",obj@meta.data$orig.ident[1])) + 
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
        theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
        theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
        theme(legend.text = element_text(size=8, hjust=0.5, color="black")) +
        theme(legend.title = element_blank()) +
        theme(axis.text = element_text(size=12, hjust=0.5, color="black")))
dev.off()

# Find neighbours and clusters
obj <- RunDR.v3(seurat.obj = obj, npcs = num_pc)
obj <- FindNeighbors(object = obj, dims = 1:num_pc, force.recalc = T)
for(i in c(0.6,1.2)){
  obj <- FindClusters(object = obj, algorithm = 1, resolution = i)
  scTool_PlotTSNE.v3.multi(seurat.obj = obj, clustname = paste0("RNA_snn_res.",i))
}

# obj <- scTool_addmodscores.v3(seuratobj = obj, geneset.rds = "/home/mbarcella/Refs/CustomGeneSets.rds", 
#                               runaddmodule = T, reduction.name = "tsne")
# scTool_addmodscores.v3(seuratobj = obj, geneset.rds = "/home/mbarcella/Refs/CustomGeneSets.rds", 
#                        runaddmodule = F, reduction.name = "umap")
# 
# 
obj@misc[["markers"]] <- scTool_FindAllMarkers.v3(seurat.obj = obj, diff.test = c("MAST"),
                                                  clustnames = c("RNA_snn_res.0.6", "RNA_snn_res.1.2"), pval.t = 1e-05)
obj@misc[["markerspos"]] <- scTool_FindAllMarkers.v3(seurat.obj = obj, diff.test = c("MAST"),
                                                     clustnames = c("RNA_snn_res.0.6", "RNA_snn_res.1.2"), pval.t = 1e-05, 
                                                     onlypos = T, outdir = "Markers_onlypos")

saveRDS(object = obj@misc[["markers"]], file = paste0(sample_id, "_allmarkers.rds"))
saveRDS(object = obj@misc[["markerspos"]], file = paste0(sample_id, "_posmarkers.rds"))

# save markers to excel files
for(clusts in c("RNA_snn_res.0.6", "RNA_snn_res.1.2")){
  tag <- gsub("res.", "", gsub("RNA_", "", clusts))
  CreateMultiTabExcels.v3(myinput = paste0(sample_id,"_allmarkers.rds"), reso = clusts, test.used = c("MAST"),
                          tag = tag, seurat.obj = obj)
  CreateMultiTabExcels.v3(myinput = paste0(sample_id,"_posmarkers.rds"), reso = clusts, test.used = c("MAST"),
                          tag = tag, seurat.obj = obj, outfold = "Excel_markers_pos")
  
}

saveRDS(object = obj, file = paste0(sample_id,"_logNormalized.rds"))
# save.image(file = paste0(obj@meta.data$orig.ident[1],"_logNormalized.Rda"))
