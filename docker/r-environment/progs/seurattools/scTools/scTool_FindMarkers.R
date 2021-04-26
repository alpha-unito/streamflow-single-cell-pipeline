scTool_FindMarkers <- function(seurat.obj, clustnames, ntop=c(10,50), 
                               log.t = 0.25, pval.t = 0.01,
                               diff.test="wilcox", outdir="Markers", 
                               covs=NULL, mincellsgroup = 10,
                               pseudocounts = 1){
  
  library(Seurat)
  library(RColorBrewer)
  library(dplyr)
  
  all <- list()
  
  for(stat.test in diff.test){
  
    for(clustname in clustnames){
    
    obj <- SetAllIdent(object = seurat.obj, id = clustname)
    all[[clustname]][[stat.test]] <- FindAllMarkers(object = obj, min.pct = 0.1, logfc.threshold = log.t, 
                          min.cells.group = mincellsgroup, test.use = stat.test, return.thresh = pval.t,
                          pseudocount.use = pseudocounts)
    
        for(ngenes in ntop){
          cur.dir <- paste0(outdir, "/", clustname, "/", stat.test, "/top" , ngenes)
          dir.create(path = cur.dir, recursive = T, showWarnings = F)
          top <- all[[clustname]][[stat.test]] %>% group_by(cluster) %>% top_n(ngenes, avg_logFC)
          png(filename = paste0(cur.dir, "/", "HeatMap_",seurat.obj@project.name,"_",clustname, "_pseudo_", pseudocounts, ".png"), width = 1200, height = 800)
          print(DoHeatmap(object = obj, genes.use = top$gene, slim.col.label = TRUE, remove.key = TRUE))
          dev.off()
      
            # cycling over cov for combined heatmaps
          
          if(!is.null(covs)){
            
            for(i in covs){
        
              cur.var <- paste0(clustname, "_", i) # creating on-fly variable
              obj@meta.data[[cur.var]] <- paste0(obj@ident, "_", obj@meta.data[[i]])
              
              png(filename = paste0(cur.dir, "/", "HeatMap_",seurat.obj@project.name,"_",clustname, "_by_", i, ".png"), width = 1200, height = 800)
              print(DoHeatmap(object = obj, genes.use = top$gene, slim.col.label = TRUE, remove.key = TRUE, group.by = cur.var))
              dev.off()
            }
          }
      
        }
    }
  }  
  
  return(all)
}