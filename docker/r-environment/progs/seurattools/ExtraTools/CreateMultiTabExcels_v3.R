# Testing multitab excel creation

CreateMultiTabExcels.v3 <- function(myinput, id=NULL, seurat.obj, 
                                 reso=c("RNA_snn_res.0.6"),
                                 test.used="MAST", tag="snn_0.6",
                                 outfold="Excel_markers"){
  
  library(openxlsx)
  library(Seurat)
  
  if(is.null(id)){
    
    id <- seurat.obj@meta.data$orig.ident[1]
  }
  
  dir.create(path = outfold, showWarnings = F)
  
  my.rds <- readRDS(file = myinput)

  for(k in test.used){
    
    of <- paste0(outfold, "/", id, "_", k, "_", tag ,".xlsx")
    OUT <- createWorkbook()
    
    for(i in reso)
    {
      for(clust in 0:(length(unique(my.rds[[i]][[k]]$cluster))-1)){
        
        mydf <- my.rds[[i]][[k]]
        mydf <- mydf[mydf$cluster == clust,]
        sname <-paste(id, "_", tag, "_c", clust, sep = "")
        addWorksheet(OUT, sname)
        writeData(OUT, sheet = sname, x = mydf)
      }
    }
    saveWorkbook(OUT,of)
  }

}

