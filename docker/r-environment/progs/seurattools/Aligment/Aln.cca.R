Aln.cca <- function(objlist, gene2use, projectname=NULL, n.cc=30){
  
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(gridExtra)
  library(grid)
  library(clustree)
  library(Radviz)
  library(e1071)
  
  start_time <- Sys.time()
  
  dir.create(path = "QC", showWarnings = F)
  my.QC <- "QC"
  
  if(length(objlist) > 2){
    obj <- RunMultiCCA(object.list = objlist, genes.use = gene2use, num.cc = n.cc)
  }
  else{
    obj <- RunCCA(object = objlist[[1]], object2 = objlist[[2]], genes.use = gene2use, num.cc = n.cc)
  }
 
  obj@project.name <- projectname
  
  # plotting MetageneBicorPlot for choosing n.ccs
  
  #png(filename = paste0(my.QC, "/", obj@project.name,"_MetageneBicorPlot.png"), width = 800, height = 600)
  #MetageneBicorPlot(object = obj, grouping.var = "orig.ident", dims.eval = 1:n.cc, smooth = T, display.progress = T) + ggtitle(label = projectname)
  #dev.off()
  
  # end time
  end_time <- Sys.time()
  total_time <- end_time - start_time
  print(paste0("Running time is: ", total_time))
  
  return(obj)
}
