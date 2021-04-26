Aln.alnsb_dr.R <- function(seurat.obj, n.cc, groupvar='orig.ident', clustprefix="celltype", 
                           res.start=0.6, res.final=2, res.step=0.2){
  
  library(grid)
  library(gridExtra)
  # Function to perform Alignsubspace function followed by different dimensional reduction calculation and 
  # cluster finding
  
  # This function takes as input several parameters:
  
  # 1. seurat object to align
  # 2. group variable by which the aligning is performed
  # 3. number of CC dimensions to use (inspect MetageneBicorPlot to choose a value)
  # 4. cluster prefix that preceed different resolutions
  # 5. resolutions settings
  
  
  obj <- AlignSubspace(object = seurat.obj, reduction.type = "cca", 
                         grouping.var = groupvar, 
                         dims.align = 1:n.cc)
  
  obj <- RunTSNE(obj, reduction.use = "cca.aligned", dims.use = 1:n.cc, do.fast = T, seed.use = 123, dim.embed = 2)
  obj <- RunUMAP(object = obj, dims.use = 1:n.cc, reduction.use = "cca.aligned", seed.use = 123, max.dim = 2)
  #obj <- RunPHATE(object = obj)
  
  # FindClusters 
  
  umap.f <- "DimReductions/UMAP/"
  dir.create(path = umap.f, showWarnings = F, recursive = T)
  tsne.f <- "DimReductions/tSNE/"
  dir.create(path = tsne.f, showWarnings = F, recursive = T)
  
  resolution <- seq(from=0.6, to=2, by=0.2)
  for (i in 1:length(resolution)){
    obj <- FindClusters(obj, reduction.type = "cca.aligned", resolution = resolution[i], 
                        dims.use = 1:n.cc, force.recalc = TRUE)
    
    # Plotting tSNE
    png(filename = paste0(tsne.f,"/",obj@project.name, "_tSNE_",resolution[i],"_aligned.png"), width = 12, height = 9, res=300, units = "in")
    print(TSNEPlot(obj, pt.size = 2.5, do.label = T, label.size = 12) + ggtitle(label = paste0("tSNE ", obj@project.name), subtitle = paste0("Aligned", clustprefix, "\nres: ", resolution[i])) + 
        theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5)))
    dev.off()
    
    # Plotting UMAP
    png(filename = paste0(umap.f,"/",obj@project.name, "_UMAP_",resolution[i],"_aligned.png"), width = 12, height = 9, res=300, units = "in")
    print(DimPlot(object = obj, reduction.use = 'umap', pt.size = 2.5, do.label = T, label.size = 12) + ggtitle(label = paste0("UMAP ", obj@project.name), subtitle = paste0("Aligned", clustprefix, "\nres: ", resolution[i])) + 
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
    dev.off()
    
    # Plotting tSNE + tSNE.orig ident side by side
    
    png(filename = paste0(tsne.f,"/",obj@project.name, "_tSNE_",resolution[i],"_with_orig.ident_aligned.png"), width = 12, height = 9, res=300, units = "in")
    a <- TSNEPlot(obj, pt.size = 2.5, do.label = T, label.size = 12) + ggtitle(label = paste0("tSNE ", obj@project.name), subtitle = paste0("Aligned ", clustprefix, "\nres: ", resolution[i])) + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5))
    b <- TSNEPlot(obj, pt.size = 2.5, do.label = T, label.size = 12, group.by="orig.ident") + ggtitle(label = paste0("tSNE ", obj@project.name), subtitle = paste0("Aligned original identity ", "\nres: ", resolution[i])) + 
      theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5))
    grid.arrange(a, b, ncol=2)
    dev.off()
    
    
    # Plotting UMAP + UMAP.orig ident side by side
    
    png(filename = paste0(umap.f,"/",obj@project.name, "_UMAP_",resolution[i],"_with_orig.ident_aligned.png"), width = 12, height = 9, res=300, units = "in")
    a <- DimPlot(object = obj, reduction.use = 'umap', pt.size = 2.5, do.label = T, label.size = 12) + ggtitle(label = paste0("UMAP ", obj@project.name), subtitle = paste0("Aligned", clustprefix, "\nres: ", resolution[i])) + 
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    b <- DimPlot(object = obj, reduction.use = 'umap', pt.size = 2.5, do.label = T, label.size = 12, group.by="orig.ident") + ggtitle(label = paste0("UMAP ", obj@project.name), subtitle = paste0("Aligned original identity ", "\nres: ", resolution[i])) + 
      theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5))
    grid.arrange(a, b, ncol=2)
    dev.off()
    
    obj <- StashIdent(obj, save.name = paste(clustprefix,resolution[i], sep=""))
    
  }
  
  return(obj)
  
}

