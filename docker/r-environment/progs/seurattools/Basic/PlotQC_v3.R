PlotQC.v3 <- function(seurat.obj, rawQC=TRUE, outfold = "QC/Raw/"){
  
  library(stringr)
  library(ggplot2)
  
  id <- seurat.obj@meta.data$orig.ident[1]

  if(rawQC == TRUE){
    tag <- "raw"
    dir.create(path = outfold, showWarnings = FALSE, recursive = T) # create output directory
  }
  
  if(rawQC == FALSE){
    tag <- "postfiltering"
    outfold <- "QC/PostFiltering/"
    dir.create(path = outfold, showWarnings = FALSE, recursive = T) # create output directory
  }

  # % MITO density plot
  
  png(filename = paste(outfold, id,"_pctMT_density_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(percent.mt)) +
    geom_density() +
      ggtitle(label = "Density plot %Mitochondrial genes", subtitle = paste0("Sample ",id," RawData")) + 
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
      theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
      theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
      theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
      theme(legend.title = element_blank()) +
      theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
    scale_x_log10())
  dev.off()
  
  # % MITO density plot
  
  png(filename = paste(outfold, id,"_pctRPs_density_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(percent.rp)) +
          geom_density() +
          ggtitle(label = "Density plot %Ribosomal proteins", subtitle = paste0("Sample ",id," RawData")) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
          scale_x_log10())
  dev.off()
  
  # Scatterplot with 2D density estimation for comparing MT vs RPP percentages
  
  png(filename = paste(outfold, id,"_pctMT_vs_pctRP_scatter_with2D_density_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(percent.rp, percent.mt)) +
    geom_point()+
    geom_density_2d() +
    ggtitle(label = "Scatter with 2D density (%RP vs %MT)", subtitle = paste0("Sample ",id," RawData")) + 
      theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
      theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
      theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
      theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
      theme(legend.title = element_blank()) +
      theme(axis.text = element_text(size=12, hjust=0.5, color="black"))
    )
  dev.off()
  
  # % nUMI vs nGene scatterplot
  
  png(filename = paste(outfold, id,"_nUMI_vs_nGenes_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(nCount_RNA, nFeature_RNA)) +
          geom_point() +
          stat_smooth(method = "gam", formula = y ~ s(x), size = 1) +
          ggtitle(label = "Scatterplot nUMIs vs nGenes", subtitle = paste0("Sample ",id," RawData")) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
          scale_x_log10())
  dev.off()
  
  # % nUMI vs %MT scatterplot
  
  png(filename = paste(outfold, id,"_nUMI_vs_pctMT_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(nCount_RNA, percent.mt)) +
          geom_point() +
          stat_smooth(method = "gam", formula = y ~ s(x), size = 1) +
          ggtitle(label = "Scatterplot nUMIs vs %MT", subtitle = paste0("Sample ",id," RawData")) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
          scale_x_log10())
  dev.off()
  
  # % nUMI vs %RP scatterplot
  
  png(filename = paste(outfold, id,"_nUMI_vs_pctRP_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(nCount_RNA, percent.rp)) +
          geom_point() +
          stat_smooth(method = "gam", formula = y ~ s(x), size = 1) +
          ggtitle(label = "Scatterplot nUMIs vs %RP", subtitle = paste0("Sample ",id," RawData")) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
          scale_x_log10())
  dev.off()
  
  # % nGene vs %MT scatterplot
  
  png(filename = paste(outfold, id,"_nGenes_vs_pctMT_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(nFeature_RNA, percent.mt)) +
          geom_point() +
          stat_smooth(method = "gam", formula = y ~ s(x), size = 1) +
          ggtitle(label = "Scatterplot nGenes vs %MT", subtitle = paste0("Sample ",id," RawData")) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
          scale_x_log10())
  dev.off()
  
  # % nGene vs %RP scatterplot
  
  png(filename = paste(outfold, id,"_nGenes_vs_pctRP_",tag,".png",sep = ""), 
      units = 'in', width = 12, height = 9, res=200)
  print(ggplot(seurat.obj@meta.data, aes(nFeature_RNA, percent.rp)) +
          geom_point() +
          stat_smooth(method = "gam", formula = y ~ s(x), size = 1) +
          ggtitle(label = "Scatterplot nGenes vs %RP", subtitle = paste0("Sample ",id," RawData")) + 
          theme(plot.title = element_text(size = 16, face = "bold", hjust = 0.5)) +
          theme(plot.subtitle=element_text(size=12, hjust=0.5, face="italic", color="black")) +
          theme(axis.title = element_text(size=12, hjust=0.5, face="bold", color="black")) +
          theme(legend.text = element_text(size=12, hjust=0.5, color="black")) +
          theme(legend.title = element_blank()) +
          theme(axis.text = element_text(size=12, hjust=0.5, color="black")) +
          scale_x_log10())
  dev.off()
  
}