scTool_PlotTSNE.v3.multi <- function(seurat.obj, id=NULL, 
                                     umap.f=NULL, tsne.f=NULL, clustname="RNA_snn_res.1", 
                                     pointsize=1, labsize=10){
  
  library(RColorBrewer)
  library(plotly)
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
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
  
# Plotting tSNE 2D

png(filename = paste0(tsne.f,"2D/",id, "_tSNE_",clustname, ".png"), units = "in", width = 12, height = 9, res = 200)
print(DimPlot(obj, pt.size = pointsize, label = T, label.size = labsize, reduction="tsne", group.by = clustname) + 
        ggtitle(label = paste0("tSNE ", obj@project.name, " Sample ", id), subtitle = paste0("Cluster: ",clustname)) + 
        theme(plot.title = element_text(hjust = 0.5, size = 12), plot.subtitle = element_text(hjust = 0.5)))
dev.off()

# Plotting UMAP
png(filename = paste0(umap.f,"2D/",id, "_UMAP_",clustname, ".png"), units = "in", width = 12, height = 9, res = 200)
print(DimPlot(object = obj, reduction="umap", pt.size = pointsize, label = T, label.size = labsize, group.by = clustname) +
      ggtitle(label = paste0("UMAP ", obj@project.name, " Sample ", id), subtitle = paste0("Cluster: ", clustname)) +
        theme(plot.title = element_text(hjust = 0.5), plot.subtitle = element_text(hjust = 0.5)))
dev.off()


# 3d plot tSNE
  
p <- plot_ly(obj@meta.data, x = ~tSNE3d_1, y = ~tSNE3d_2, z = ~tSNE3d_3, 
             color =~ obj@meta.data[[clustname]], colors = "Paired") %>%
    add_markers(size=3) %>%
    layout(title=paste0("3D tSNE ",id," cluster:", clustname), scene = list(xaxis = list(title = 'tSNE_1'),
                                                                yaxis = list(title = 'tSNE_2'),
                                                                zaxis = list(title = 'tSNE_3')))
  
htmlwidgets::saveWidget(widget = as_widget(p), file = paste0(tsne.f,"3D/", id, "_tSNE_",clustname,".html")) 
  
# 3d plot UMAP

p <- plot_ly(obj@meta.data, x = ~UMAP3d_1, y = ~UMAP3d_2, z = ~UMAP3d_3, 
             color = ~ obj@meta.data[[clustname]], colors = "Paired") %>%
    add_markers(size=3) %>%
    layout(title=paste0("3D UMAP ",id," cluster:", clustname), scene = list(xaxis = list(title = 'UMAP1'),
                                                                yaxis = list(title = 'UMAP2'),
                                                                zaxis = list(title = 'UMAP3')))

htmlwidgets::saveWidget(as_widget(p), paste0(umap.f,"3D/", id, "_UMAP_",clustname,".html"))

}

