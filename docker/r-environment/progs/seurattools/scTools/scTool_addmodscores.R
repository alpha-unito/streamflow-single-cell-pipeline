scTool_addmodscores <- function(seuratobj, geneset.rds, runaddmodule=TRUE, seed=1234, nctrls=100,
                                outdir="FeaturePlots", imputing=NULL, reduction.name="tSNE.2d") {
  
  rds.geneset <- readRDS(geneset.rds)

  for (lname in names(rds.geneset)) {
    
    if(isTRUE(runaddmodule)){
      
      seuratobj <-
        AddModuleScore(object = seuratobj, genes.list = rds.geneset[[lname]], random.seed = seed, ctrl.size = nctrls)
      colnames(seuratobj@meta.data)[(length(colnames(seuratobj@meta.data)) - (length(names(rds.geneset[[lname]])) -
                                                                                1)):length(colnames(seuratobj@meta.data))] <-
        paste0(lname,"_", names(rds.geneset[[lname]]))
      
    }
    
    # create dirs
    curdir <- paste0(outdir, "/")
    dir.create(path = curdir, showWarnings = FALSE, recursive = T)
    
    # do featplot
    
    png(filename = paste0(curdir, seuratobj@project.name, "_", lname, "_featplot_", imputing, ".png"), width = 12, height = 9, res=300, units = "in")
    FeaturePlot(object = seuratobj, features.plot=paste0(lname,"_", names(rds.geneset[[lname]])), reduction.use = reduction.name)
    dev.off()
  }
  
  return(seuratobj)
}
