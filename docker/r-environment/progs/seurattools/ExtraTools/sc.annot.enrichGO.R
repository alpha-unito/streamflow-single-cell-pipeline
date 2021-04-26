# script for annotatating markers and list of genes

sc.annot.enrichGO <- function(seurat.rds, markers.rds, specie="human", 
                              res="res.0.6", test="MAST", go="BP"){
  
  library(clusterProfiler)
  library(DOSE)
  library(org.Hs.eg.db)
  library(org.Mm.eg.db)
  
  db <- NULL
  
  if (specie == "human"){
    
    db <- "org.Hs.eg.db"
  }
  if (specie == "mouse"){
    db <- "org.Mm.eg.db"
  }
  
  markers <- readRDS(markers.rds)
  seurat.obj <- readRDS(seurat.rds)
  
  tot.clusts <- levels(markers[[res]][[test]]$cluster)
  tmp.markers <- markers[[res]][[test]]
  
   ego <- list()

   for(i in tot.clusts){
     markers <- tmp.markers[tmp.markers$cluster == i,]
     myset <- bitr(geneID = markers$gene, fromType="SYMBOL",
                   toType=c("ENSEMBL", "ENTREZID"),
                   OrgDb=db)
     myuni <- bitr(geneID = row.names(seurat.obj@data), fromType="SYMBOL",
                   toType=c("ENSEMBL", "ENTREZID"),
                   OrgDb=db)

     ego[[paste0("Cluster_n_",i)]] <- as.data.frame(enrichGO(gene = myset$ENTREZID,
                                   universe = myuni$ENTREZID,
                                   OrgDb    = db,
                                   ont      = go,
                                   pAdjustMethod = "BH",
                                   pvalueCutoff  = 0.01,
                                   qvalueCutoff  = 0.05,
                                   readable = TRUE))

   }

   ego$detest <- test
   ego$resolution <- res
   return(ego)

}
