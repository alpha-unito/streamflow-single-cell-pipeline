scTool_violin_v3 <- function(seurat.obj, featuresvect, groupvar="orig.ident", outfold="QC/Violin/", customcol="#079FAB", id=NULL){
  
  library(Seurat)
  library(RColorBrewer)
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  dir.create(path = outfold, showWarnings = F, recursive = T)
  outname <- paste0(outfold, paste(featuresvect, collapse="_"), "_",seurat.obj@project.name, "_", id, "_by_", groupvar,".png")
  
  ncharts <- length(featuresvect)
  
  png(filename = outname, units = "in", width = 12, height = 9, res=100)
  print(VlnPlot(object = obj, features = featuresvect, cols = customcol, pt.size = 0.2, log = TRUE))
  dev.off()
  
}