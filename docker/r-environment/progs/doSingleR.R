doSingleR <- function(seurat.rds, projname, specie="Human", ofold="SingleR_analysis_noclustering/", 
                      vars2embed = c("RNA_snn_res.0.6","RNA_snn_res.1.2","orig.ident"), 
                      topmains = 30, topsubpops=30, cores, orderbycluster=T){
  
  library(SingleR)
  library(RColorBrewer)
  library(Seurat)
  
  seurat.obj <- readRDS(seurat.rds)

  projfold <- paste0(ofold, "/")
  dir.create(path = projfold, recursive = T)
  
  singler = CreateSinglerObject(counts = seurat.obj@assays$RNA@data, annot = NULL, project.name = projname, 
                                min.genes = 0, technology = "10X", species = specie, citation = "",
                                normalize.gene.length = F, variable.genes = "de",
                                fine.tune = T, do.signatures = T, clusters = NULL, do.main.types = T, 
                                reduce.file.size = T, numCores = cores)
  
  
  ## if using Seurat v3.0 and over use:
  
  singler$seurat = seurat.obj # (optional)
  singler$meta.data$xy = seurat.obj@reductions$tsne@cell.embeddings # the tSNE coordinates
  singler$meta.data$xy.umap = seurat.obj@reductions$umap@cell.embeddings # the UMAP coordinates
  # ## if using a previous Seurat version use:
  # singler$meta.data$xy = seurat.object@dr$tsne@cell.embeddings # the tSNE coordinates
  # singler$meta.data$clusters = obj@meta.data$RNA_snn_res.1.2 # the Seurat clusters (if 'clusters' not provided)
  
  if(isTRUE(orderbycluster)){
    tag <- "noclustering"
  }
  else{
    tag <- "hclustering"
  }
  
  for(i in vars2embed){

    singler$meta.data[i] = seurat.obj@meta.data[i]
    
    png(filename = paste0(projfold, "HM1_top", topmains, "_mains_", i, "_", tag, "_annotation.png"), width = 12, height = 9, units = "in", res = 500)
    SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single.main, top.n= topmains,
                        clusters = singler$meta.data[i], normalize = T, 
                        order.by.clusters = orderbycluster)
    dev.off()
    
    png(filename = paste0(projfold, "HM1_top", topsubpops, "_refined_", i, "_", tag, "_annotation.png"), width = 12, height = 9, units = "in", res = 400)
    SingleR.DrawHeatmap(singler$singler[[1]]$SingleR.single, top.n= topsubpops,
                        clusters = singler$meta.data[i], normalize = T, 
                        order.by.clusters = orderbycluster)
    dev.off()
    
    png(filename = paste0(projfold, "HM2_top", topmains, "_mains_", i, "_", tag, "_annotation.png"), width = 12, height = 9, units = "in", res = 400)
    SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single.main, top.n= topmains,
                        clusters = singler$meta.data[i], normalize = T, 
                        order.by.clusters = orderbycluster)
    dev.off()
    
    png(filename = paste0(projfold, "HM2_top", topsubpops, "_refined_", i, "_", tag, "_annotation.png"), width = 12, height = 9, units = "in", res = 400)
    SingleR.DrawHeatmap(singler$singler[[2]]$SingleR.single, top.n= topsubpops,
                        clusters = singler$meta.data[i], normalize = T, 
                        order.by.clusters = orderbycluster)
    dev.off()
  }
  
  df_singleR_1_main <- as.data.frame(singler$singler[[1]]$SingleR.single.main$labels)
  base::colnames(x = df_singleR_1_main)[1] <- "Celltype.main.1"
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = df_singleR_1_main)
  
  df_singleR_1_refined <- as.data.frame(singler$singler[[1]]$SingleR.single$labels)
  base::colnames(x = df_singleR_1_refined)[1] <- "Celltype.refined.1"
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = df_singleR_1_refined)
  
  df_singleR_2_main <- as.data.frame(singler$singler[[2]]$SingleR.single.main$labels)
  base::colnames(x = df_singleR_2_main)[1] <- "Celltype.main.2"
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = df_singleR_2_main)
  
  df_singleR_2_refined <- as.data.frame(singler$singler[[2]]$SingleR.single$labels)
  base::colnames(x = df_singleR_2_refined)[1] <- "Celltype.refined.2"
  seurat.obj <- AddMetaData(object = seurat.obj, metadata = df_singleR_2_refined)
  
  for(datasets in c("Celltype.main.1", "Celltype.refined.1", "Celltype.main.2", "Celltype.refined.2")){
    colourCount = length(unique(seurat.obj@meta.data[[datasets]]))
    getPalette = colorRampPalette(brewer.pal(n = 9, name = "Set1"))
    for (j in vars2embed) {
      png(filename = paste0(projfold, datasets, "_umap_split_by_", j, ".png"), width = 12, height = 9, units = "in", res = 400)
      print(DimPlot(object = seurat.obj, reduction = "umap", pt.size = 0.5, group.by = datasets, cols = getPalette(colourCount), split.by = j))
      dev.off()
    }
    png(filename = paste0(projfold, datasets,"_umap.png"), width = 12, height = 9, units = "in", res = 400)
    print(DimPlot(object = seurat.obj, reduction = "umap", pt.size = 0.5, group.by = datasets, cols = getPalette(colourCount)))
    dev.off()
  }
  saveRDS(object = singler, file = paste(projfold, "singleR.rds", sep = "/"))
  return(singler)
}

args <- commandArgs(trailingOnly = TRUE)
seurat.rds <- args[1]
sample <- args[2]
species <- args[3]
out_dir <- args[4]
threads <- args[5]

sr <- doSingleR(seurat.rds, sample, species, out_dir, c("RNA_snn_res.0.6", "RNA_snn_res.1.2", "orig.ident"),
                topmains = 30, topsubpops=30, threads, orderbycluster=T)
