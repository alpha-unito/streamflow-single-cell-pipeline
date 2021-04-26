ScaleandRunPCA <-
  function(seurat.obj, ncores = 1, out = "PCAs/", regvarsvector, id=NULL, cellcycle.rds, npcs=30) {
    
    if(is.null(id)){
      
      id <- str_split(string = seurat.obj@cell.names[1], pattern = "_")[[1]][1]
    }
    
    cc.genes <- readRDS(file = cellcycle.rds)
    s.genes <- cc.genes[["s.genes"]]
    g2m.genes <- cc.genes[["g2m.genes"]]
    
    
    dir.create(out, showWarnings = FALSE) # create outfolder
    
    seurat.obj <-
      ScaleData(
        object = seurat.obj,
        vars.to.regress = regvarsvector,
        display.progress = T,
        do.par = T,
        num.cores = ncores
      )
    
    # Run PCA using only cell-cycle genes
    seurat.obj <-
      RunPCA(
        object = seurat.obj,
        pc.genes = c(s.genes, g2m.genes),
        do.print = FALSE, pcs.compute = npcs
      )
    png(
      filename = paste(
        out,
        id,
        "_PCA_PCElbow_CellCycle_Var_Genes_", paste(regvarsvector, collapse = "_") , ".png", sep = ""
      ),
      width = 1200,
      height = 800
    )
  
    
    cc <- PCAPlot(object = seurat.obj) + theme(plot.title = element_text(hjust = 0.5),
                                           plot.subtitle = element_text(hjust = 0.5)) + 
                                     ggtitle(label = "", subtitle = "PCA Cell Cycle Genes")
    cc.elbow <- PCElbowPlot(object = seurat.obj, num.pc = npcs) + theme(plot.title = element_text(hjust = 0.5),
                                               plot.subtitle = element_text(hjust = 0.5)) +
                                     ggtitle(label = "", subtitle = "PCElbow Cell Cycle Genes")

    # Run PCA using most variable genes
    
    seurat.obj <-
      RunPCA(
        object = seurat.obj,
        pc.genes = seurat.obj@var.genes,
        do.print = FALSE, pcs.compute = npcs
      )
    
    vargenes <-
      PCAPlot(object = seurat.obj) + theme(plot.title = element_text(hjust = 0.5),
                                           plot.subtitle = element_text(hjust = 0.5)) +
                                     ggtitle(label = "", subtitle = "PCA VarGenes Genes")
    vargenes.elbow <-
      PCElbowPlot(object = seurat.obj, num.pc = npcs) + theme(plot.title = element_text(hjust = 0.5),
                                               plot.subtitle = element_text(hjust = 0.5)) +
                                     ggtitle(label = "", subtitle = "PCA VarGenes Genes")
    grid.arrange(
      cc,
      cc.elbow,
      vargenes,
      vargenes.elbow,
      nrow = 2,
      ncol = 2,
      top = textGrob(paste(id, " Scaling ", paste(regvarsvector, collapse = ","),sep=""), gp=gpar(fontsize=20,font=1)) 
    )
    dev.off()
    return(seurat.obj)
    
  }
