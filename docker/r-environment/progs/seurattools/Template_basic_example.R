# Load data

source(file = "/home/mbarcella/scripts/bitbucket/seurattools/Basic/DataLoading_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/Basic/PlotQC_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_violin_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/Basic/Normalize_FindHVG.v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/Basic/CellScoringMeta_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_plotViz_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_plotPCA_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/Basic/RunDR_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_dimHeatmap_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_plotTSNE_v3_multi.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_FindAllMarkers_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/scTools/scTool_addmodscores_v3.R")
source(file = "/home/mbarcella/scripts/bitbucket/seurattools/ExtraTools/CreateMultiTabExcels_v3.R")


library(Seurat)
library(ggplot2)
library(RColorBrewer)

# Import data

obj <- DataLoading.v3(fold10x = "filtered_feature_bc_matrix/", 
                      projname = "AML", id = "M22", 
                      mincells = 10, minfeatures = 50)

saveRDS(object = obj, file = paste0(obj@meta.data$orig.ident[1], "_minimal.rds"))

# QC 

PlotQC.v3(seurat.obj = obj, rawQC = TRUE)

scTool_violin_v3(seurat.obj = obj, featuresvect = c("nFeature_RNA", "nCount_RNA"))
scTool_violin_v3(seurat.obj = obj, featuresvect = c("percent.mt", "percent.rp"))

# Filtering

obj <- subset(x = obj, subset = nFeature_RNA > 200 & nFeature_RNA < 7000 
              & percent.mt < 20 & percent.rp < 100)

# QC post filtering

PlotQC.v3(seurat.obj = obj, rawQC = FALSE)

# OPTION 1 LOG-NORMALIZATION + SCALING

# Normalize and find HVGs

obj <- NormalizeFindHVG.v3(seurat.obj = obj, top.feat2plot = 50, 
                           n.features = 2000, selmethod = "vst") # find variable genes is embedded here
# cell-cycle scoring

obj <- CellScoringMeta.v3(cellcycle.rds = "/DATA_NFS/Metadata/Human/scRNA/regev_lab_cell_cycle_human.rds", 
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

obj <- RunDR.v3(seurat.obj = obj, npcs = 20)
obj <- FindNeighbors(object = obj, dims = 1:20, force.recalc = T)

for(i in c(0.6,1.2)){
  obj <- FindClusters(object = obj, algorithm = 1, resolution = i)
  scTool_PlotTSNE.v3.multi(seurat.obj = obj, clustname = paste0("RNA_snn_res.",i))
}


obj <- scTool_addmodscores.v3(seuratobj = obj, geneset.rds = "/home/mbarcella/Refs/CustomGeneSets.rds", 
                              runaddmodule = T, reduction.name = "tsne")
scTool_addmodscores.v3(seuratobj = obj, geneset.rds = "/home/mbarcella/Refs/CustomGeneSets.rds", 
                       runaddmodule = F, reduction.name = "umap")


obj@misc[["markers"]] <- scTool_FindAllMarkers.v3(seurat.obj = obj, diff.test = c("MAST", "wilcox"),
                                                  clustnames = c("RNA_snn_res.0.6", "RNA_snn_res.1.2"), pval.t = 1e-05)
obj@misc[["markerspos"]] <- scTool_FindAllMarkers.v3(seurat.obj = obj, diff.test = c("MAST", "wilcox"),
                                                     clustnames = c("RNA_snn_res.0.6", "RNA_snn_res.1.2"), pval.t = 1e-05, 
                                                     onlypos = T, outdir = "Markers_onlypos")

saveRDS(object = obj@misc[["markers"]], file = paste0(obj@meta.data$orig.ident[1], "_allmarkers.rds"))
saveRDS(object = obj@misc[["markerspos"]], file = paste0(obj@meta.data$orig.ident[1], "_posmarkers.rds"))

# save markers to excel files

for(clusts in c("RNA_snn_res.0.6", "RNA_snn_res.1.2")){
  
  CreateMultiTabExcels.v3(myinput = paste0(obj@meta.data$orig.ident[1],"_allmarkers.rds"), reso = clusts, test.used = c("MAST", "wilcox"),
                          tag = clusts, seurat.obj = obj)
  CreateMultiTabExcels.v3(myinput = paste0(obj@meta.data$orig.ident[1],"_posmarkers.rds"), reso = clusts, test.used = c("MAST", "wilcox"),
                          tag = clusts, seurat.obj = obj, outfold = "Excel_markers_pos")
  
}

saveRDS(object = obj, file = paste0(obj@meta.data$orig.ident[1],"_logNormalized.rds"))
save.image(file = paste0(obj@meta.data$orig.ident[1],"_logNormalized.Rda"))