

CompareClust <- function(obj1.rds,
                         obj2.rds, 
                         obj1_clusters=c("SNN_UNSUP_Cluster_0.6", "SNN_UNSUP_Cluster_1.2"), 
                         lfc_thresh = 0.25, 
                         p_thresh = 0.01, 
                         obj2_clusters = "SNN_UNSUP_Cluster_1.2",
                         outname="ClustCompare.png", seed = 1234, 
                         nctrls = 100){
  
  library(Seurat)
  library(ggplot2)
  library(ggrepel)
  library(gplots)
  library(data.table)
  
  obj1 <- readRDS(obj1.rds)
  obj2 <- readRDS(obj2.rds)
  
  myvect <- NULL
  for(i in obj2_clusters){
    
    myvect <- c(myvect, i)
    
  }
  
  # cycle over n.clusters configurations:
  # 1. for each cluster resolution find the markers
  # 2. with each list of genes intersect the second object and add a module score to it based on such list of genes
  # 3. Calculate the value of the added module score by each cluster defined in each obj2 configuration
  # 4. Stores all mean values in a dataframe
  # 5. The final data frame will incorporate on the rows the clusters in obj1 and in the cols the clusters in obj2.
  # 6. The value on n-row and m-col will be filled by the average value of the module score in m cluster in obj2 
  #    build with the genes that mark the cluster n in object1
  
  for (clustconf in 1:length(obj1_clusters)) {
    
    # set cluster configuration in order to perform FindAllMarkers specific to that resolution
    obj1 <- SetAllIdent(object = obj1, id = obj1_clusters[[clustconf]])
    # retrieve the list of all markers that determine each cluster
    obj1.markers <- FindAllMarkers(object = obj1)
    obj1.markers.copy <<- obj1.markers
    # filter out genes based on padj and lfc
    if(lfc_thresh > 0){
      obj1.markers <- obj1.markers[obj1.markers$p_val_adj < p_thresh & obj1.markers$avg_logFC > lfc_thresh, ]
    }
    else{
      obj1.markers <- obj1.markers[obj1.markers$p_val_adj < p_thresh & obj1.markers$avg_logFC < lfc_thresh, ]
    }
    
    for(obj1.curclust in levels(obj1.markers$cluster)){
      
      # select only genes obtained by FindAllMarkers in obj1 that are present in obj2
      cur.genes <- intersect(obj1.markers[obj1.markers$cluster == obj1.curclust, ][["gene"]], 
                             rownames(obj2@scale.data))
      obj2 <- AddModuleScore(obj2, genes.list = list(cur.genes), random.seed = seed, ctrl.size = nctrls,
                             enrich.name = paste("obj1", "n", obj1.curclust, obj1_clusters[[clustconf]], sep="_"))
      print(paste("obj1", "n", obj1.curclust, obj1_clusters[[clustconf]], sep="_"))
    }
    
  }

  # retrieve all 
  extended_metadata <- obj2@meta.data[,c(myvect ,
                            grep(pattern = "^obj1_n", x = colnames(obj2@meta.data), perl = T, value = T))]
  
  print(grep(pattern = "^obj1_n", x = colnames(obj2@meta.data), perl = T, value = T))
  
   colnames(extended_metadata) <- gsub(pattern = "$0.61", replacement = "0.6", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$0.81", replacement = "0.8", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$11", replacement = "1", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$1.21", replacement = "1.2", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$1.41", replacement = "1.4", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$1.61", replacement = "1.6", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$1.81", replacement = "1.8", x = colnames(extended_metadata), perl = T)
   colnames(extended_metadata) <- gsub(pattern = "$21", replacement = "2", x = colnames(extended_metadata), perl = T)
   
   

   extended_metadata_cells <<- extended_metadata
   extended_metadata_cells$cnames <<- rownames(extended_metadata_cells)
   
   meta.dt <- data.table(extended_metadata)


  meta.dt.melted <- melt.data.table(data = meta.dt, measure.vars = as.vector(myvect))

  meta.dt.melted$celltype_clust <- paste("obj2", meta.dt.melted$value, meta.dt.melted$variable,  sep = "_n_")

  meta.dt.melted$value <- NULL
  meta.dt.melted$variable<- NULL

  meta.dt.melted$celltype_clust <- as.factor(meta.dt.melted$celltype_clust)
  
  vectlist <- NULL
  
  for(i in 1:(ncol(meta.dt.melted)-1)){
    
    colname <- colnames(meta.dt.melted)[i]
    vectlist[[colname]] <- meta.dt.melted[, median(get(colname)), by = celltype_clust]
    colnames(vectlist[[colname]])[2] <- colname
    vectlist[[colname]] <- t(vectlist[[colname]])
  }
  combined <- do.call("rbind", vectlist) 
  colnames(combined) <- combined[1,]
  combined.unique <- unique(combined)
  
  combined.unique <- combined.unique[2:nrow(combined.unique),]
  
  combined.unique.mt <- as.matrix(apply(combined.unique, 2, as.numeric))
  rownames(combined.unique.mt) <- row.names(combined.unique)
  
  library(gplots)
  png(filename = outname, width = 1200, height = 800)
  heatmap.2(x = combined.unique.mt, scale = 'none', density.info = 'none', 
                margins=c(25,25), trace = 'none', 
                col=bluered(75))
  dev.off()
  
  return(list(combined.unique.mt, obj2))

}


