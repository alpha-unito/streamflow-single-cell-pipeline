ScaleandRunPCA.v3 <- function(seurat.obj, ncores = 1, out = "PCAs/", 
                              regvarsvector = c("CC.difference", "percent.mt"), 
                              id=NULL, cellcycle.rds, npcs=30) {
    
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
    
    cc.genes <- readRDS(file = cellcycle.rds)
    s.genes <- cc.genes[["s.genes"]]
    g2m.genes <- cc.genes[["g2m.genes"]]
    
    dir.create(out, showWarnings = FALSE) # create outfolder
    
    all.genes <- rownames(seurat.obj)
    seurat.obj <- ScaleData(object = seurat.obj, vars.to.regress = regvarsvector, 
                            display.progress = T, use.umi = T)
    
    # Run PCA using only cell-cycle genes
    
    # seurat.obj <- RunPCA(object = seurat.obj, features = c(s.genes, g2m.genes), 
    #                      do.print = FALSE, verbose = T, npcs = npcs)
    # 
     png(filename = paste(out, id, "_PCA_PCElbow_CellCycle_Var_Genes_",
                         paste(regvarsvector, collapse = "_") , ".png", sep = ""),
         width = 1200, height = 800)
    # 
    # 
    # cc <- PCAPlot(object = seurat.obj) + theme(plot.title = element_text(hjust = 0.5),
    #                                        plot.subtitle = element_text(hjust = 0.5)) + 
    #                                  ggtitle(label = "", subtitle = "PCA Cell Cycle Genes")
    # cc.elbow <- PCElbowPlot(object = seurat.obj, num.pc = npcs) + theme(plot.title = element_text(hjust = 0.5),
    #                                            plot.subtitle = element_text(hjust = 0.5)) +
    #                                  ggtitle(label = "", subtitle = "PCElbow Cell Cycle Genes")
    # 
    # # Run PCA using most variable genes
    # 
     seurat.obj <- RunPCA(object = seurat.obj, features = VariableFeatures(object = seurat.obj), npcs = npcs)
    # 
     vargenes <- PCAPlot(object = seurat.obj) + theme(plot.title = element_text(hjust = 0.5),
                                            plot.subtitle = element_text(hjust = 0.5)) +
                                      ggtitle(label = "", subtitle = "PCA VarGenes Genes")
    # 
     vargenes.elbow <- PCElbowPlot(object = seurat.obj, num.pc = npcs) + theme(plot.title = element_text(hjust = 0.5),
                                                plot.subtitle = element_text(hjust = 0.5)) +
                                      ggtitle(label = "", subtitle = "PCA VarGenes Genes")
    # 
    # 
     grid.arrange(vargenes, vargenes.elbow, nrow = 1, ncol = 2, 
                  top = textGrob(paste(id, " Scaling ", paste(regvarsvector, collapse = ","),sep=""), 
                                 gp=gpar(fontsize=20,font=1)))
     dev.off()
    # 
    return(seurat.obj)
    
  }
