NormalizeFindHVG.v3 <- function(seurat.obj, id=NULL, outfold="QC/", 
                                n.features=2000, top.feat2plot=20, selmethod = "vst"){
  
  require(gridExtra)
  require(grid)
  library(openxlsx)
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  # Normalization
  
  seurat.obj <- NormalizeData(object = seurat.obj)
  
  # Find most variable genes (Feature selection)
  
  seurat.obj <- FindVariableFeatures(seurat.obj, selection.method = selmethod,  
                                     nfeatures = n.features, verbose = T)
  top <- head(VariableFeatures(seurat.obj), top.feat2plot)
  plot1 <- VariableFeaturePlot(seurat.obj)
  plot2 <- LabelPoints(plot = plot1, points = top, repel = TRUE, xnudge = 0, ynudge=0) + scale_y_log10() +
    ggtitle(label = paste0("Top ", top.feat2plot, " High Variable Genes"), subtitle = paste0("Sample ",id)) + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
    theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
    theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
    theme(legend.title = element_blank()) +
    theme(axis.text = element_text(size=12, hjust=0.5, color="black"))
  
  png(filename = paste0(outfold, id, "_HVG",".png"), 
      units = "in", width = 12, height = 9, res=100)
  print(plot2)
  dev.off()
  
  write.xlsx(x = VariableFeatures(seurat.obj), file = paste0("VariableFeatures_",id,".xlsx"), rownames=T, firstRow=T, firstCol=T, 
             creator="MatteoBarcella", tabColour="#af0f35", rowNames=TRUE, sheetName=paste0(id,"_VariableFeatures"))
  
  
  return(seurat.obj)
}