scTool_pcaplotting <- function(seurat.obj, outfold="QC_charts/PCA/", n.pcs, ...){
  
  library(gridExtra)
  library(grid)
  library(Seurat)
  
  dir.create(path = outfold, showWarnings = F, recursive = T)

  
  png(filename = paste0(outfold,"PCAPlot_with_elbow_",seurat.obj@project.name,".png"), width = 800, height = 800, res = 75)
  a <- PCAPlot(object = seurat.obj, 1, 2) + 
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    
  b <- PCElbowPlot(object = obj, num.pc = n.pcs) + 
    theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
  
  #ggtitle(label = "PCA plot + Elbow", subtitle = seurat.obj@project.name)
  
  grid.arrange(a,b, nrow=2, top = textGrob(label = paste0("PCA plot + Elbow\n", seurat.obj@project.name), just = "center")) 
  
  dev.off()
  
}