scTool_FindAllMarkers.v3 <- function(seurat.obj, clustnames, ntop=c(10,50), 
                               log.t = 0.25, pval.t = 0.01, by.var=NULL, covs=NULL,
                               min.pct = 0.1, diff.test="wilcox", outdir="Markers", mincellsgroup = 10, onlypos=FALSE){
  
  library(Seurat)
  library(RColorBrewer)
  library(dplyr)
  
  all <- list()
  
  for(stat.test in diff.test){
    
    for(clustname in clustnames){
      
      if(!is.null(covs)){
        
        obj <- SetIdent(object = seurat.obj, value = clustname)
        all[[clustname]][[stat.test]] <- FindAllMarkers(object = obj, min.pct = min.pct, logfc.threshold = log.t, 
                                                        only.pos = onlypos,min.cells.group = mincellsgroup, 
                                                        test.use = stat.test, return.thresh = pval.t, latent.vars = covs)
      }else{
        obj <- SetIdent(object = seurat.obj, value = clustname)
        all[[clustname]][[stat.test]] <- FindAllMarkers(object = obj, min.pct = min.pct, logfc.threshold = log.t, only.pos = onlypos,
                                                        min.cells.group = mincellsgroup, test.use = stat.test, return.thresh = pval.t)
      }
      
      for(ngenes in ntop){
        cur.dir <- paste0(outdir, "/", clustname, "/", stat.test, "/top" , ngenes)
        dir.create(path = cur.dir, recursive = T, showWarnings = F)
        top <- all[[clustname]][[stat.test]] %>% group_by(cluster) %>% top_n(ngenes, avg_logFC)
        png(filename = paste0(cur.dir, "/", "HeatMap_",seurat.obj@meta.data$orig.ident[1],"_",clustname, ".png"), 
            width = 12, height = 9, units = 'in', res=300)
        print(DoHeatmap(object = obj, features = top$gene, label = TRUE, 
                        raster = FALSE) + NoLegend())
        dev.off()
        
        if(!is.null(by.var)){
          
          png(filename = paste0(cur.dir, "/", "HeatMap_",seurat.obj@meta.data$orig.ident[1],"_",clustname, "_by_",by.var,".png"), 
              width = 12, height = 9, units = 'in', res=200)
          print(DoHeatmap(object = obj, features = top$gene, label = TRUE, 
                          group.bar = T, group.by = by.var, raster = FALSE) + NoLegend())
          dev.off()
          
        }
      }
    }
  }  
  
  return(all)
}