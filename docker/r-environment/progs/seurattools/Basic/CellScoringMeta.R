CellScoringMeta <- function(id=NULL, cellcycle.rds ,outfold="QC_charts/", seurat.obj, genevector, ncolumns){
  
  library(reshape2)
  
  if(is.null(id)){
    
    id <- str_split(string = seurat.obj@cell.names[1], pattern = "_")[[1]][1]
  }
  
  # load rds and get genes
  cc.genes <- readRDS(file = cellcycle.rds)
  s.genes <- cc.genes[["s.genes"]]
  g2m.genes <- cc.genes[["g2m.genes"]]
  
  #genevector <- c(s.genes[1:3], g2m.genes[1:3])
  
  # Performing cellcyclescoring
  
  seurat.obj <- CellCycleScoring(object = seurat.obj, s.genes = s.genes, g2m.genes = g2m.genes, set.ident = TRUE)
  seurat.obj@meta.data$CC.Difference <- seurat.obj@meta.data$S.Score - seurat.obj@meta.data$G2M.Score
  png(filename = paste(outfold, id, "_CellScoring_frequency_Table.png", sep = ""), width = 800, height = 600)
  freq_table <- prop.table(x = table(seurat.obj@ident, seurat.obj@meta.data[, "Phase"]))
  barplot(freq_table, main = paste(id," Cell Scoring Frequency Distribution", sep=""))
  dev.off()
  #png(filename = paste(outfold,id,"_RidgePlot_cellcycle_Markers.png", sep = ""), width = 800, height = 600)
  #print(RidgePlot(object = seurat.obj, features.plot = genevector, nCol = ncolumns))
  #dev.off()
  df <- FetchData(object = seurat.obj, vars.all = c("G2M.Score", "S.Score"))
  df.m <- melt(df)
  png(filename = paste(outfold,id,"_CellCycleScore_Density.png", sep = ""), 
      width = 800, height = 600)
  print(ggplot(df.m) + geom_density(aes(x = value, colour = variable)) + labs(x = NULL) + 
    ggtitle(paste(id, " Densities (Kernel Estimator)", sep = ""))) 
  dev.off()
  return(seurat.obj)
}
