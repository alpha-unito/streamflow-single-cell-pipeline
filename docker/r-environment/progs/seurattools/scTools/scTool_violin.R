scTool_violin <- function(seurat.obj, featuresvect, groupvar, setwidth=800, setheight=600, outfold="QC_charts/"){
  
  library(Seurat)
  
  dir.create(path = outfold, showWarnings = F, recursive = T)
  outname <- paste0(outfold, paste(featuresvect, collapse="_"), seurat.obj@project.name, "by", groupvar,".png")
  
  png(filename = outname, width = setwidth, height = setheight)
  print(VlnPlot(object = seurat.obj, features.plot = featuresvect, group.by = groupvar, do.return = TRUE))
  dev.off()
}