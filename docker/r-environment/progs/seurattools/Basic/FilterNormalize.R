FilterNormalize <- function(seurat.obj, id=NULL, low.ngenes, high.ngenes, high.pctmito, low.numi, high.numi,
                            scalefactor=10000, norm.method="LogNormalize", 
                            outfold="QC_charts/", bins=20){
  
  require(gridExtra)
  require(grid)
  
  if(is.null(id)){
    
    id <- str_split(string = seurat.obj@cell.names[1], pattern = "_")[[1]][1]
  }
  
  seurat.obj <- FilterCells(object = seurat.obj, subset.names = c("nUMI","nGene", "pct.mito"),
                            low.thresholds = c(low.numi, low.ngenes, -Inf), high.thresholds = c(high.numi, high.ngenes, high.pctmito))
  seurat.obj <- NormalizeData(object = seurat.obj, normalization.method = norm.method, scale.factor = scalefactor)
  seurat.obj <- FindVariableGenes(object = seurat.obj, mean.function = ExpMean, 
                                  dispersion.function = LogVMR, x.low.cutoff = 0.0125, x.high.cutoff = 3, 
                                  y.cutoff = 0.5, do.plot = T, num.bin = bins)
  
  png(filename = paste(outfold, id, "_Mean_Dispersion_plots.png", sep = ""), width = 1200, height = 800)
  par(mfrow = c(1, 2))
  plot(x=seurat.obj@hvg.info$gene.mean, y = seurat.obj@hvg.info$gene.dispersion, 
       main = paste(id, " Mean dispersion plot",sep="" ))
  plot(x=seurat.obj@hvg.info$gene.mean, y = seurat.obj@hvg.info$gene.dispersion.scaled, 
       main = paste(id, " Mean dispersion plot - scaled", sep="" ))
  dev.off()

  return(seurat.obj)
}