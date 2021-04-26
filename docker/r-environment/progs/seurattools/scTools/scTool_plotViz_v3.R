scTool_PlotViz.v3 <- function(seurat.obj, outfold = "QC/PCA/", id=NULL, viz.feat=15, 
                              start.pc=1, end.pc=2, reduction="pca", ncols=NULL,
                              pca.identity="orig.ident", pca.ptsize=0.2){
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  dir.create(outfold, showWarnings = FALSE, recursive = T) # create outfolder
  
  png(filename = paste0(outfold, id, "_VizLoadings","_PC_", start.pc, "_to_", end.pc,".png"), units = "in", width = 12, height = 9, res=100)
  print(VizDimLoadings(object = seurat.obj, dims = start.pc:end.pc, reduction = reduction, nfeatures = viz.feat, ncol = ncols) + 
    ggtitle(label = paste0("\nVizplot loadings PC", start.pc, " to PC", end.pc, "\n"), subtitle = paste0("Sample ",id)) + 
    theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
    theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
    theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
    theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
    theme(legend.title = element_blank()) +
    theme(axis.text = element_text(size=12, hjust=0.5, color="black")))
  dev.off()
 
}