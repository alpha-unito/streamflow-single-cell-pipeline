CellScoringMeta.v3 <- function(id=NULL, cellcycle.rds,
                            outfold="QC/", 
                            seurat.obj, ncolumns){
  
  library(reshape2)
  library(RColorBrewer)
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  # load rds
  
  cc.genes <- readRDS(file = cellcycle.rds)
  s.genes <- cc.genes[["s.genes"]]
  g2m.genes <- cc.genes[["g2m.genes"]]
  
  seurat.obj@misc[["cc.genes"]] <- c(s.genes, g2m.genes)
  
  
  # Performing cellcyclescoring
  
  seurat.obj <- CellCycleScoring(object = seurat.obj, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
  seurat.obj@meta.data$CC.Difference <- seurat.obj@meta.data$S.Score - seurat.obj@meta.data$G2M.Score
  
  png(filename = paste(outfold, id, "_CellScoring_frequency_Table.png", sep = ""), units = 'in', width = 12, height = 9, res=100)
  freq_table <- prop.table(x = table(seurat.obj@meta.data$orig.ident, seurat.obj@meta.data[, "Phase"]))
  print(freq_table)
  barplot(freq_table, col = brewer.pal(n = 4, "Set1"), main = paste(id," Cell Scoring Frequency Distribution", sep=""))
  dev.off()
  
  df <- FetchData(object = seurat.obj, vars = c("G2M.Score", "S.Score"))
  df.m <- melt(df)
  
  png(filename = paste(outfold,id,"_CellCycleScore_Density.png", sep = ""),  units = 'in', width = 12, height = 9, res=100)
  print(ggplot(df.m) + geom_density(aes(x = value, fill = variable), alpha=0.4) + labs(x = NULL) + 
    ggtitle(paste(id, " Densities", sep = ""))) 
  dev.off()
  return(seurat.obj)
}
