scTool_DimHeatmap.v3 <- function(seurat.obj, outfold = "QC/PCA/", id=NULL,
                                 reduction="pca", slot.data="scale.data", 
                                 dimheat.n.cells=500, dimheat.end=9, dimheat.feat=10){
  
  library(RColorBrewer)
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  dir.create(outfold, showWarnings = FALSE, recursive = T) # create outfolder
  
  myname <- paste0(outfold, id, "_Heatmap_PC1_to_PC_",dimheat.end,"_reduction_",reduction,"_Slot_",slot.data,".png")
  
  png(filename = myname, units = "in", width = 12, height = 9, res=100)
  print(DimHeatmap(object = seurat.obj, dims = 1:dimheat.end, cells = dimheat.n.cells, nfeatures = dimheat.feat,
                   reduction = reduction, slot = slot.data, fast = FALSE, raster = FALSE, balanced = TRUE) + 
          ggtitle(label = paste0("\nDimHeatmap PC1 to PC", dimheat.end, " Top ", dimheat.n.cells, " cells \n"), 
                  subtitle = paste0("Sample ",id, "\nslot: ", slot.data, "\nreduction: ",reduction)) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")))
  dev.off()
  
}