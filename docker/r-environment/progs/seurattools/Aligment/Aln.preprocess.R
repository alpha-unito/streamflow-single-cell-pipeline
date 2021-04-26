Aln.preprocess <- function(rdsfiles, rdsnames, tophvg_pct, imputing.method="noimputing", metadesign=NULL, idvarinmeta="Sample"){
  
  listobjs <- list() # variable to store objects
  mostvar <- list() # store most var genes 

  # Setting tophvg value if less than 10  

  if(tophvg_pct < 10){
    tophvg_pct <- paste0(0,tophvg_pct)
  }
 
  # Cycle over rds files, takes the most var genes across them and build merged object that
  # will undergo to the alignment procedure

  for(i in 1:length(rdsfiles)){
    
    listobjs[[rdsnames[[i]]]] <- readRDS(file = rdsfiles[[i]]) # read rds
    
    if(!is.null(metadesign)){
      
      design <- read.table(file = metadesign, header = T, sep = ",")
      
      tmp.metadata <- merge.data.frame(x = listobjs[[rdsnames[[i]]]]@meta.data, y = design, by.x = 'orig.ident', by.y = idvarinmeta)
      rownames(tmp.metadata) <- rownames(listobjs[[rdsnames[[i]]]]@meta.data)
      listobjs[[rdsnames[[i]]]]@meta.data <- tmp.metadata
    }
    
    if(imputing.method == "alra"){
      
      # loading functions to perform alra imputation
      
      source(file = "https://raw.githubusercontent.com/KlugerLab/ALRA/master/alra.R")
      source(file = "https://raw.githubusercontent.com/KlugerLab/ALRA/master/alraSeurat2.R")
      
       sink(file = 'Alra_report.log', append=TRUE)
       print(paste0("Object: ", rdsnames[[i]]))
       listobjs[[rdsnames[[i]]]] <- alraSeurat2(obj = listobjs[[rdsnames[[i]]]])
       sink()
    }

    #listobjs[[rdsnames[[i]]]]@meta.data$cond <- rdsnames[[i]]  # store in the cond variable the original id
    
    # store in the mostvar list the most variable genes according to selected threshold tophvg_pct
    # the threshold is the quantile: if tophvg_pct variable is 10, the firs 10% of the genes are used.
    
    mostvar[[rdsnames[[i]]]] <- head(rownames(listobjs[[rdsnames[[i]]]]@hvg.info), 
                                     round(quantile(x = seq(1:length(rownames(listobjs[[rdsnames[[i]]]]@hvg.info))), 
                                                    probs = as.numeric(paste0(0,".",tophvg_pct))), digits = 0))
    
  }
  
  # compute the uniquely genes among the lists of genes (from different objects)
  # basically we remove duplicated genes that are present in more than 1 list
  
  union_hgvs <- unique(unlist(mostvar)) 
  
  # Select genes that are most variable in each sample and that are in common across samples
  
  genes.use <- intersect(union_hgvs, rownames(listobjs[[1]]@scale.data))
  
  for(k in 2:length(listobjs)){
    genes.use <- intersect(genes.use, rownames(listobjs[[k]]@scale.data)) # cycle over the list of objects
  }
  
  # the function returns a list with 3 elements:
  # a list with loaded objects, a list with the most variable genes for each object, a vector of genes 
  # to be used in the Canonical Correlation Analysis (CCA) in seurat.

  return(list(objects=listobjs, hgvinfos = mostvar, gene2use = genes.use)) 
  
}





