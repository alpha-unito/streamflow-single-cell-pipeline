scTool_FindAllMarkers.v3 <- function(seurat.obj, clustnames, ntop=c(10,50), 
                               log.t = 0.25, pval.t = 0.01, by.var=NULL,
                               min.pct = 0.1, diff.test="wilcox", outdir="Markers", 
                               covs=NULL, mincellsgroup = 10){
  
  library(Seurat)
  library(RColorBrewer)
  library(dplyr)
  
  all <- list()
  pos <- list()
  
  for(stat.test in diff.test){
    
    for(clustname in clustnames){
      
      obj <- SetIdent(object = seurat.obj, value = clustname)
      all[[clustname]][[stat.test]] <- FindAllMarkers(object = obj, min.pct = min.pct, logfc.threshold = log.t, 
                                                      min.cells.group = mincellsgroup, test.use = stat.test, return.thresh = pval.t)
      pos[[clustname]][[stat.test]] <- FindAllMarkers(object = obj, min.pct = min.pct, logfc.threshold = log.t, only.pos = TRUE,
                                                      min.cells.group = mincellsgroup, test.use = stat.test, return.thresh = pval.t)
      
      for(ngenes in ntop){
        cur.dir <- paste0(outdir, "/", clustname, "/", stat.test, "/top" , ngenes)
        dir.create(path = cur.dir, recursive = T, showWarnings = F)
        top <- pos[[clustname]][[stat.test]] %>% group_by(cluster) %>% top_n(ngenes, avg_logFC)
        png(filename = paste0(cur.dir, "/", "HeatMap_",seurat.obj@project.name,"_",clustname, ".png"), 
            width = 12, height = 9, units = 'in', res=200)
        print(DoHeatmap(object = obj, features = top$gene, label = TRUE,
                        raster = FALSE) + NoLegend())
        dev.off()
        
        if(!is.null(by.var)){
          
          png(filename = paste0(cur.dir, "/", "HeatMap_",seurat.obj@project.name,"_",clustname, "_by_",by.var,".png"), 
              width = 12, height = 9, units = 'in', res=200)
          print(DoHeatmap(object = obj, features = top$gene, label = TRUE, 
                          group.bar = T, group.by = by.var, raster = FALSE) + NoLegend())
          dev.off()
          
        }
      }
    }
  }  
  
  return(list(allmarkers=all, posmarkers=pos))
}