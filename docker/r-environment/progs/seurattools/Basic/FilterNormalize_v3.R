FilterNormalizeFindHVG.v3 <- function(seurat.obj, id=NULL, 
                               max.genes=5000, max.pctMT=10, min.UMI=100, max.pctRP=100, 
                               outfold="QC/", n.features=2000, top.feat2plot=20){
  
  require(gridExtra)
  require(grid)
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  # filtering 
  seurat.obj <- subset(x = seurat.obj, subset = nFeature_RNA > min.UMI & nFeature_RNA < max.genes & percent.mt < max.pctMT & percent.rp < max.pctRP)
  
  # Normalization
  seurat.obj <- NormalizeData(object = seurat.obj)
  
  # Find most variable genes (Feature selection)
  
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = "vst", nfeatures = n.features, verbose = T)
  
  png(filename = paste0(outfold, "High_variable_features_", seurat.obj$project.name, "_", id,".png"), 
      units = "in", width = 12, height = 9, res=100)
  plot1 <- VariableFeaturePlot(seurat.obj)
  plot1 <- LabelPoints(plot = plot1, points = top.feat2plot, repel = TRUE)
  plot1
  dev.off()
  
  write.xlsx(x = VariableFeatures(seurat.obj), file = paste0("VariableFeatures_",id,".xlsx"), rownames=T, firstRow=T, firstCol=T, 
             creator="MatteoBarcella", tabColour="#af0f35", rowNames=TRUE, sheetName=paste0(id,"_VariableFeatures"))
  
  
  return(seurat.obj)
}