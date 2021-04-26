# This function produce radviz plots starting from a seurat object and a list of cell programs lists
# For each list a radviz plot is produced

scTool_radvizplots <- function(seuratobj, geneset.rds, clust2fetch, outdir="RadvizPlots", psize=1){
  
  library(Seurat)
  require(Radviz)
  require(RColorBrewer)
  library(e1071)
  
  # load lists
  
  rds.geneset <- readRDS(geneset.rds)
  
  for (i in clust2fetch){
    
    curdir <- paste0(outdir, "/", clust2fetch, "/")
    
    dir.create(path = curdir, showWarnings = FALSE, recursive = T)
    
    seuratobj <- SetAllIdent(object = seuratobj, id = clust2fetch)
    
    # set color palette according to the number of clusters
    
    coul <- NULL
    
    if(length(levels(seuratobj@ident)) < 9){
      coul = brewer.pal(length(levels(seuratobj@ident)),"Set1")
    }
    else{
      coul = brewer.pal(9,"Set1")
      coul <- colorRampPalette(coul)(length(levels(seuratobj@ident)))
    }
    
    
    for (lname in names(rds.geneset)) {
      
      if(length(names(rds.geneset[[lname]])) < 2) next
      
      curdata <- FetchData(object = seuratobj, vars.all = paste0(lname,"_", names(rds.geneset[[lname]])))
      colnames(curdata) <- names(rds.geneset[[lname]])
      curclustannot <- FetchData(object = seuratobj, vars.all = clust2fetch)
      
      # data
      curdata_norm <- apply(curdata,2,do.L,fun =function(x) quantile(x,c(0.025,0.975)))
      curdata_ct.S <- make.S(dimnames(curdata_norm)[[2]])
      curdata_ct.sim <- cosine(curdata_norm)
      
      in.da(curdata_ct.S,curdata_ct.sim) # increases with better projections
      rv.da(curdata_ct.S,curdata_ct.sim) # decreases with better projections
      
      curdata_optim.ct <- do.optim(curdata_ct.S,curdata_ct.sim,iter=100,n=1000)
      curdata_ct.S <- make.S(tail(curdata_optim.ct$best,1)[[1]])
      curdata_ct.rv <- do.radviz(curdata_norm,curdata_ct.S)
      
      legend <- sort(unique(as.numeric(seuratobj@meta.data[[clust2fetch]])), decreasing = F)
      namesG <- table(as.character(seuratobj@meta.data[[clust2fetch]]))
      namesG1 <- as.character(unique(seuratobj@meta.data[[clust2fetch]]))
      
    
      png(filename = paste0(curdir, "RadarPlot_", seuratobj@project.name, "_", lname, "_", clust2fetch ,".png"), width = 1200,height = 800)
      
      plot(curdata_ct.rv, point.shape=18, point.size=psize, point.color = coul[curclustannot[[clust2fetch]]],label.size = 2,label.col= "black",cex=1)
      legend('topright', legend = paste(legend, " (", as.character(namesG[as.character(sort(as.numeric(names(namesG)),decreasing = F))]),")",sep=""),
             col=coul[0:length(legend)],pch=18,bty='n',cex=1.5)
      title(paste0(lname, " scores"))
      dev.off()
      
      
    }
  }
  

}