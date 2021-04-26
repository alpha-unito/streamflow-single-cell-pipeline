PlotQC <- function(seurat.obj, rawQC=TRUE, outfold = "QC_charts/"){
  
  library(stringr)
  dir.create(path = outfold, showWarnings = FALSE)
  
  id <- str_split(string = seurat.obj@cell.names[1], pattern = "_")[[1]][1]
  
  tag <- "prefiltering"
  
  if(rawQC == FALSE){
    tag <- "postfiltering"
  }

  png(filename = paste(outfold, id, "_pct_mitochondrial_and_nGene_density_",tag,".png",sep = ""),
      width = 800,
      height = 600)
  par(mfrow = c(1, 2))
  plot(density(
    x = seurat.obj@meta.data$pct.mito,
    from = 0,
    to = 30
  ), main = paste(id, " Density % mitochondrial content", sep = ""), lwd=2)
  abline(v = 5, col = "blue", lwd=2, lty=3)
  abline(v = 10, col = "red", lwd=2, lty=3)
  plot(density(x = seurat.obj@meta.data$nGene), main = paste(id," Density nGene", sep = ""), lwd=2)
  dev.off()
  
  png(filename = paste(outfold, id, "_nUMI_vs_pct_mitochondrial_and_nGenes",tag,".png", sep = ""),
      width = 800,
      height = 600)
  par(mfrow = c(1, 2))
  GenePlot(
    object = seurat.obj,
    gene1 = "nUMI",
    gene2 = "pct.mito",
    pch.use = 3,
    col.use = "red"
  )
  GenePlot(
    object = seurat.obj,
    gene1 = "nUMI",
    gene2 = "nGene",
    pch.use = 3,
    col.use = "blue"
  )
  dev.off()
  
}