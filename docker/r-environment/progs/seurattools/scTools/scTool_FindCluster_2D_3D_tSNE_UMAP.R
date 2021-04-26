scTool_FindCluster_2D_3D_tSNE_UMAP <- function(seuratobj, rerun.dr = F , n.cc = 30, id=NULL, 
                                               umap.f=NULL, tsne.f=NULL, plot3D=FALSE, feature.3D=NULL,
                                               redu.type="pca", var2plot=NULL, onlyvarplot=FALSE, labsize =10,
                                               pointsize=1.5, onlyplot=FALSE, n.cores=4, rerun.dr.redu.type = NULL,
                                               resolution=c(0.6,0.8,1,1.2,1.4,1.6,1.8,2), mypalette="Set1"){
  
  library(plotly)
  library(reshape2)
  library(doParallel)
  
  obj <- seuratobj
  
  # setting id of sample
  
  if(is.null(id)){
    
    id <- str_split(string = obj@cell.names[1], pattern = "_")[[1]][1]
  }
  
  # if tSNE and UMAP was not calculated on 3 dimensions rerun the calculation
  
  if(is.null(umap.f) & is.null(tsne.f)){
    
    umap.f <- paste0(normalizePath("."),"/DimReductions/UMAP/")
    dir.create(path = umap.f, showWarnings = F, recursive = T)
    dir.create(path = paste0(umap.f, "2D/"), showWarnings = F, recursive = T)
    dir.create(path = paste0(umap.f, "3D/"), showWarnings = F, recursive = T)
    tsne.f <- paste0(normalizePath("."),"/DimReductions/tSNE/")
    dir.create(path = tsne.f, showWarnings = F, recursive = T)
    dir.create(path = paste0(tsne.f, "2D/"), showWarnings = F, recursive = T)
    dir.create(path = paste0(tsne.f, "3D/"), showWarnings = F, recursive = T)
  }
  
  # RuntSNE and RunUMAP and plot 2D and 3D
  
  if(rerun.dr == TRUE){

    if(is.null(rerun.dr.redu.type)){
      
      obj <- RunTSNE(object = obj, dims.use = 1:n.cc, do.fast = T, seed.use = 123, reduction.name = "tSNE.2d", reduction.key = "tSNE.2d_")
      obj <- RunUMAP(object = obj, dims.use = 1:n.cc, seed.use = 123, max.dim = 2, reduction.name = "UMAP.2d", reduction.key = "UMAP.2d_")
      obj <- RunTSNE(object = obj, dims.use = 1:n.cc, do.fast = T, seed.use = 123, dim.embed = 3, reduction.name = "tSNE.3d", reduction.key = "tSNE.3d_")
      obj <- RunUMAP(object = obj, dims.use = 1:n.cc, seed.use = 123, max.dim = 3, reduction.name = "UMAP.3d", reduction.key = "UMAP.3d_")
      
    }
    else{
    
      obj <- RunTSNE(object = obj, dims.use = 1:n.cc, do.fast = T, seed.use = 123, reduction.name = "tSNE.2d", reduction.key = "tSNE.2d_", reduction.use = rerun.dr.redu.type)
      obj <- RunUMAP(object = obj, dims.use = 1:n.cc, seed.use = 123, max.dim = 2, reduction.name = "UMAP.2d", reduction.key = "UMAP.2d_", reduction.use = rerun.dr.redu.type)
      obj <- RunTSNE(object = obj, dims.use = 1:n.cc, do.fast = T, seed.use = 123, dim.embed = 3, reduction.name = "tSNE.3d", reduction.key = "tSNE.3d_", reduction.use = rerun.dr.redu.type)
      obj <- RunUMAP(object = obj, dims.use = 1:n.cc, seed.use = 123, max.dim = 3, reduction.name = "UMAP.3d", reduction.key = "UMAP.3d_", reduction.use = rerun.dr.redu.type)
    }
    

    # store cell embeddings
    
    tsne.cell_emb_2d <- as.data.frame(obj@dr$tSNE.2d@cell.embeddings)
    umap.cell_emb_2d <- as.data.frame(obj@dr$UMAP.2d@cell.embeddings)
    tsne.cell_emb_3d <- as.data.frame(obj@dr$tSNE.3d@cell.embeddings)
    umap.cell_emb_3d <- as.data.frame(obj@dr$UMAP.3d@cell.embeddings)

    # add embeddings to metadata 
    
    obj <- AddMetaData(object = obj, metadata = tsne.cell_emb_2d)
    obj <- AddMetaData(object = obj, metadata = umap.cell_emb_2d)
    
    obj <- AddMetaData(object = obj, metadata = tsne.cell_emb_3d)
    obj <- AddMetaData(object = obj, metadata = umap.cell_emb_3d)
  }
  

  
  # Run findcluster and tSNE / UMAP plotting
  
  registerDoParallel(cores=n.cores)
  
  for (i in resolution){
    if(onlyplot == FALSE){
      
      obj <- FindClusters(obj, reduction.type = redu.type, resolution = i, 
                          dims.use = 1:n.cc, force.recalc = TRUE, save.SNN = T, 
                          random.seed = 123, print.output = FALSE)
      
     # obj <- ValidateClusters(obj, pc.use = 1:n.cc)
      
      freq_table <- prop.table(x = table(obj@meta.data[,paste0("res.",i)], obj@meta.data[, "Phase"]), margin = 2)
      freq_table <- melt(data = freq_table)
      freq_table$Var1 <-  as.factor(freq_table$Var1)
      freq_table$Var2 <-  as.factor(freq_table$Var2)
      
      barplot.path <- paste0(normalizePath("."),"/QC_charts/Barplots/res.",i, "/")
      dir.create(path = barplot.path, showWarnings = F, recursive = T)
      
      png(filename = paste0(barplot.path ,obj@project.name, "barplot_",i, ".png"), width = 800, height = 600)
      print(ggplot(data = freq_table, mapping = aes(x = Var2, y = value, fill=Var1)) + geom_bar(stat = "identity") +
              ggtitle(label = paste0(id, " Barplot Cell Cycle"), subtitle = paste0("Resolution: ", i)))
      dev.off()

    # Plotting tSNE 2D
    
    png(filename = paste0(tsne.f,"2D/",obj@project.name, "_tSNE_",i, ".png"), units = "in", width = 12, height = 9, res = 300)
    print(DimPlot(obj, pt.size = pointsize, do.label = T, label.size = labsize, reduction.use="tSNE.2d") + ggtitle(label = paste0("tSNE ", obj@project.name), subtitle = paste0("Resolution: ", i)) + 
            theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5)))
    dev.off()
    
    # Plotting UMAP
    png(filename = paste0(umap.f,"2D/",obj@project.name, "_UMAP_",i, ".png"), units = "in", width = 12, height = 9, res = 300)
    print(DimPlot(object = obj, reduction.use="UMAP.2d", pt.size = pointsize, do.label = T, label.size = labsize) + ggtitle(label = paste0("UMAP ", obj@project.name), subtitle = paste0("Resolution: ", i)) + 
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
    dev.off()
    
    
    if(plot3D == TRUE){
      
      if(is.null(feature.3D)){
        cur.res <- paste0("res.",i)
      }
      else{
        cur.res <- feature.3D      
        }
      
      
      # 3d plot tSNE
      
      p <- plot_ly(obj@meta.data, x = ~tSNE.3d_1, y = ~tSNE.3d_2, z = ~tSNE.3d_3, color =~ obj@meta.data[[cur.res]], colors = "Paired") %>%
        add_markers(size=3) %>%
        layout(title=paste0("3D tSNE ",id," res:", i), scene = list(xaxis = list(title = 'tSNE_1'),
                                                                        yaxis = list(title = 'tSNE_2'),
                                                                        zaxis = list(title = 'tSNE_3')))
      
      htmlwidgets::saveWidget(widget = as_widget(p), file = paste0(tsne.f,"3D/", obj@project.name, "_tSNE_",i,".html")) 
      
      # 3d plot UMAP
      
      p <- plot_ly(obj@meta.data, x = ~UMAP.3d_1, y = ~UMAP.3d_2, z = ~UMAP.3d_3, color = ~ obj@meta.data[[cur.res]], colors = "Paired") %>%
        add_markers(size=3) %>%
        layout(title=paste0("3D UMAP ",id," res:", i), scene = list(xaxis = list(title = 'UMAP1'),
                                                                        yaxis = list(title = 'UMAP2'),
                                                                        zaxis = list(title = 'UMAP3')))
      
      htmlwidgets::saveWidget(as_widget(p), paste0(umap.f,"3D/", obj@project.name, "_UMAP_",i,".html"))
      
    }
    
    registerDoSEQ()
    
    }
    
    if(onlyplot == TRUE & onlyvarplot == FALSE){
      
    obj <- SetAllIdent(object = obj, id = paste0("res.",i))
      
    #if(onlyvarplot == FALSE){
    
    # Plotting tSNE 2D
    
    png(filename = paste0(tsne.f,"2D/",obj@project.name, "_tSNE_",i, ".png"), units = "in", width = 12, height = 9, res=300)
    print(DimPlot(obj, pt.size = pointsize, do.label = T, label.size = labsize, reduction.use="tSNE.2d") + ggtitle(label = paste0("tSNE ", obj@project.name), subtitle = paste0("Resolution: ", i)) + 
            theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5)))
    dev.off()
    
    # Plotting UMAP
    png(filename = paste0(umap.f,"2D/",obj@project.name, "_UMAP_",i, ".png"), units = "in", width = 12, height = 9, res=300)
    print(DimPlot(object = obj, reduction.use="UMAP.2d", pt.size = pointsize, do.label = T, label.size = labsize) + ggtitle(label = paste0("UMAP ", obj@project.name), subtitle = paste0("Resolution: ", i)) + 
            theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
    dev.off()
    
      if(plot3D == TRUE){
      
        if(is.null(feature.3D)){
        cur.res <- paste0("res.",i)
        }
        else{
        cur.res <- feature.3D      
        }
      
      
      # 3d plot tSNE
      
        p <- plot_ly(obj@meta.data, x = ~tSNE.3d_1, y = ~tSNE.3d_2, z = ~tSNE.3d_3, color =~ obj@meta.data[[cur.res]], colors = "Paired") %>%
          add_markers(size=3) %>%
          layout(title=paste0("3D tSNE ",id," res:", i), scene = list(xaxis = list(title = 'tSNE_1'),
                                                                      yaxis = list(title = 'tSNE_2'),
                                                                      zaxis = list(title = 'tSNE_3')))
        
        htmlwidgets::saveWidget(widget = as_widget(p), file = paste0(tsne.f,"3D/", obj@project.name, "_tSNE_",i,".html")) 
        
        # 3d plot UMAP
        
        p <- plot_ly(obj@meta.data, x = ~UMAP.3d_1, y = ~UMAP.3d_2, z = ~UMAP.3d_3, color = ~ obj@meta.data[[cur.res]], colors = "Paired") %>%
          add_markers(size=3) %>%
          layout(title=paste0("3D UMAP ",id," res:", i), scene = list(xaxis = list(title = 'UMAP1'),
                                                                      yaxis = list(title = 'UMAP2'),
                                                                      zaxis = list(title = 'UMAP3')))
        
        htmlwidgets::saveWidget(as_widget(p), paste0(umap.f,"3D/", obj@project.name, "_UMAP_",i,".html"))
      
      }
    }
    }
    if(onlyplot == TRUE & onlyvarplot == TRUE & (!is.null(var2plot))){
        
        obj <- SetAllIdent(object = obj, id = var2plot)
        
        png(filename = paste0(tsne.f,"2D/",obj@project.name, "_tSNE_by_", var2plot, ".png"), units = "in", width = 12, height = 9, res=300)
        print(DimPlot(obj, pt.size = pointsize, do.label = T, label.size = labsize, group.by = var2plot, reduction.use="tSNE.2d") + 
                ggtitle(label = paste0("tSNE ", obj@project.name), subtitle = paste0("By: ", var2plot)) + 
                theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5)))
        dev.off()
        
        
        png(filename = paste0(umap.f,"2D/",obj@project.name, "_UMAP_by_", var2plot, ".png"), units = "in", width = 12, height = 9, res=300)
        print(DimPlot(object = obj, reduction.use="UMAP.2d", pt.size = pointsize, do.label = T, label.size = labsize, group.by = var2plot) + 
                ggtitle(label = paste0("UMAP ", obj@project.name), subtitle = paste0("By: ", var2plot)) + 
                theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
        dev.off()
        
        
        # 3d plot tSNE
        
        p <- plot_ly(obj@meta.data, x = ~tSNE.3d_1, y = ~tSNE.3d_2, z = ~tSNE.3d_3, color =~ obj@meta.data[[var2plot]], colors = mypalette) %>%
          add_markers(size=3) %>%
          layout(title=paste0("3D tSNE ",id," by:", var2plot), scene = list(xaxis = list(title = 'tSNE_1'),
                                                                      yaxis = list(title = 'tSNE_2'),
                                                                      zaxis = list(title = 'tSNE_3')))
        
        htmlwidgets::saveWidget(widget = as_widget(p), file = paste0(tsne.f,"3D/", obj@project.name, "_tSNE_",var2plot,".html")) 
        
        # 3d plot UMAP
        
        p <- plot_ly(obj@meta.data, x = ~UMAP.3d_1, y = ~UMAP.3d_2, z = ~UMAP.3d_3, color = ~ obj@meta.data[[var2plot]], colors = mypalette) %>%
          add_markers(size=3) %>%
          layout(title=paste0("3D UMAP ",id," by:", var2plot), scene = list(xaxis = list(title = 'UMAP1'),
                                                                      yaxis = list(title = 'UMAP2'),
                                                                      zaxis = list(title = 'UMAP3')))
        
        htmlwidgets::saveWidget(as_widget(p), paste0(umap.f,"3D/", obj@project.name, "_UMAP_",var2plot,".html"))
        
  
  }
  return(obj)
}



