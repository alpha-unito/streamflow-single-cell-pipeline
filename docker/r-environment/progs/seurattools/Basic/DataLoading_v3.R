# This function is able to load data from 10x folder

DataLoading.v3 <- function(fold10x, projname, id, mincells, minfeatures, mito.pattern="^MT-"){
  
  # load libraries
  
  library(Seurat)
  library(openxlsx)
  
  # load data from 10X folder and collect basic stats
  
  scdata <- Read10X(data.dir = fold10x) # read data
  colnames(x = scdata) <- paste(id, colnames(x = scdata), sep = '_') # add id tag to cellnames
  
  # create seurat object based on two main parameters mincells and mingenes
  
  scdata <- CreateSeuratObject(counts = scdata, min.cells = mincells, 
                               min.features = minfeatures, 
                               names.field = 1,
                               names.delim = "_",
                               project = projname) 
  
  # retrieve mitocondrial genes 

  scdata[["percent.mt"]] <- PercentageFeatureSet(scdata, pattern = mito.pattern)
  scdata[["percent.rp"]] <- PercentageFeatureSet(scdata, pattern = "^RPL|^RPS")
  
  write.xlsx(x = scdata@meta.data, file = paste0("Standard_metadata_",id,".xlsx"), rownames=T, firstRow=T, firstCol=T, 
             creator="MatteoBarcella", tabColour="#af0f35", rowNames=TRUE, sheetName=paste0(id,"_RawMetadata"))
  
  return(scdata)
  
}