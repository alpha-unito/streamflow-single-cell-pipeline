scTool_vizPCA <- function(seurat.obj, n.pcs, ncols, outfold= "QC_charts/"){

  library(Seurat)
  dir.create(path = outfold, showWarnings = F)
  outname <- paste0(outfold,"VizPCA_", seurat.obj@project.name,"_PC1_to_", n.pcs,".png")
  png(filename = outname, res = 120, width = 800, height = 600)
  VizPCA(object = obj, pcs.use = 1:n.pcs, nCol = ncols)
  dev.off()
}

