# Testing multitab excel creation

CreateMultiTabExcels <- function(myinput, id, 
                                 reso=c("res.0.6"),
                                 file.is.present = FALSE,
                                 test.used="MAST", tag="res.0.6",
                                 outfold="Excel_markers"){
  
  library(openxlsx)
  library(Seurat)
  
  dir.create(path = outfold, showWarnings = F)
  
  if(isTRUE(file.is.present)){
    file.remove(paste0(outfold, "/", id, "_", test.used, "_", tag ,".xlsx"))
  }
  
  my.rds <- readRDS(file = myinput)
  of <- paste0(outfold, "/", id, "_", test.used, "_", tag ,".xlsx")
  OUT <- createWorkbook()
  for(i in reso)
  {
    for(clust in 0:(length(unique(my.rds[[i]][[test.used]]$cluster))-1)){
      
      mydf <- my.rds[[i]][[test.used]]
      mydf <- mydf[mydf$cluster == clust,]
      sname<-paste(id,"_",i,"_clust_",clust,sep="")
      addWorksheet(OUT, sname)
      writeData(OUT, sheet = sname, x = mydf)
    }
  }
  saveWorkbook(OUT,of)
}

