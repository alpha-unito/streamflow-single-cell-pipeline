scTool_PlotPCA.v3 <- function(seurat.obj, outfold = "QC/PCA/", id=NULL,
                              start.pc=1, end.pc=2, reduction="pca",
                              pca.identity="orig.ident", pca.ptsize=1,
                              pca.splitby=NULL, pca.shapeby=NULL){
if(is.null(id)){
  
  id <- seurat.obj@meta.data$orig.ident[1]
}

dir.create(outfold, showWarnings = FALSE, recursive = T) # create outfolder

myname <- NULL

if(is.null(pca.shapeby) & is.null(pca.splitby)){
  myname <- paste0(outfold, id, "_PCAplot_",pca.identity,"_PC_", start.pc, "_vs_", end.pc,"_reduction_",reduction,".png")
}
else if(is.null(pca.shapeby)){
  myname <- paste0(outfold, id, "_PCAplot_",pca.identity,"_PC_", start.pc, "_vs_", end.pc,"_splitby_",pca.splitby,"_reduction_",reduction,".png")
}
else if(is.null(pca.splitby)){
  myname <- paste0(outfold, id, "_PCAplot_",pca.identity,"_PC_", start.pc, "_vs_", end.pc,"_shapeby_",pca.shapeby,"_reduction_",reduction,".png")
}

png(filename = myname, units = "in", width = 12, height = 9, res=100)
print(DimPlot(object = seurat.obj, reduction = reduction, group.by = pca.identity, 
                dims = c(start.pc, end.pc), pt.size = pca.ptsize, split.by = pca.splitby, shape.by = pca.shapeby) + 
        ggtitle(label = paste0("\nPCA plot PC", start.pc, " vs PC", end.pc, "\n"), subtitle = paste0("Sample ",id)) + 
        theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
        theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
        theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
        theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
        theme(legend.title = element_blank()) +
        theme(axis.text = element_text(size=12, hjust=0.5, color="black")))
dev.off()

}