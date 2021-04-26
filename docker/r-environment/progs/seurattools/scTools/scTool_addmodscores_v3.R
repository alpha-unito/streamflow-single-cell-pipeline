scTool_addmodscores.v3 <- function(seuratobj, geneset.rds, runaddmodule=TRUE, seed=1234, nctrls=100,
                                outdir="FeaturePlots", imputing=NULL, pt.size=0.2, 
                                reduction.name="tsne", id=NULL, brewpalette="YlGnBu") {
  
  library(RColorBrewer)
  rds.geneset <- readRDS(geneset.rds)
  
  if(is.null(id)){
    
    id <- seuratobj@meta.data$orig.ident[1]
  }

  for (lname in names(rds.geneset)) {
    
    if(isTRUE(runaddmodule)){
      
      seuratobj <-
        AddModuleScore(object = seuratobj, features = rds.geneset[[lname]], seed = seed, ctrl = nctrls)
      colnames(seuratobj@meta.data)[(length(colnames(seuratobj@meta.data)) - (length(names(rds.geneset[[lname]])) -
                                                                                1)):length(colnames(seuratobj@meta.data))] <-
        paste0(lname,"_", names(rds.geneset[[lname]]))
      
    }
    
    # create dirs
    curdir <- paste0(outdir, "/")
    dir.create(path = curdir, showWarnings = FALSE, recursive = T)
    
    # do featplot
    
    if(is.null(imputing)){
      png(filename = paste0(curdir, id, "_", lname, "_featplot_", reduction.name, ".png"), 
          width = 12, height = 9, res=300, units = "in")
      print(FeaturePlot(object = seuratobj, 
                        features=paste0(lname,"_", names(rds.geneset[[lname]])), 
                        reduction = reduction.name, pt.size = pt.size, 
                        cols = brewer.pal(9, brewpalette)))
      dev.off()
    }else{
      png(filename = paste0(curdir, id, "_", lname, "_featplot_",imputing, "_", reduction.name, ".png"), 
          width = 12, height = 9, res=200, units = "in")
      print(FeaturePlot(object = seuratobj, 
                        features=paste0(lname,"_", names(rds.geneset[[lname]])), 
                        reduction = reduction.name, pt.size = pt.size, 
                        cols = brewer.pal(9, brewpalette)))
      dev.off()
      
    }

  }
  
  return(seuratobj)
}
