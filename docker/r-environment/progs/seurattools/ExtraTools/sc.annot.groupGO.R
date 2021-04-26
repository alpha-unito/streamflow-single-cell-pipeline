# script for annotatating markers and list of genes

sc.annot.groupGO <- function(markers.rds, specie="human", res="res.0.6", test="MAST", go="BP", golevel=3){
  
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
  tot.clusts <- levels(markers[[res]][[test]]$cluster)
  tmp.markers <- markers[[res]][[test]]
  
  ggo <- list()
  
  for(i in tot.clusts){
    markers <- tmp.markers[tmp.markers$cluster == i,]
    myset <- bitr(geneID = markers$gene, fromType="SYMBOL", 
                  toType=c("ENSEMBL", "ENTREZID"), 
                  OrgDb=db)
    ggo[[paste0("Cluster_n_",i)]] <- as.data.frame(groupGO(gene = myset$ENTREZID,
                                 OrgDb    = db,
                                 ont      = go,
                                 level    = golevel,
                                 readable = TRUE))
  }
    ggo$detest <- test
    ggo$resolution <- res
    return(ggo)
}


