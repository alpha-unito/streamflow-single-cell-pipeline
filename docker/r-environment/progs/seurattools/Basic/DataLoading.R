# This function is able to load data from 10x folder

DataLoading <- function(fold10x, projname, id, mincells, mingenes){
  
  # load libraries
  
  library(Seurat)
  
  # load data from 10X folder and collect basic stats
  
  scdata <- Read10X(data.dir = fold10x) # read data
  colnames(x = scdata) <- paste(id, colnames(x = scdata), sep = '_') # add id tag to cellnames
  
  # create seurat object based on two main parameters mincells and mingenes
  
  scdata <- CreateSeuratObject(raw.data = scdata, min.cells = mincells, min.genes = mingenes, project = paste0(projname,"_",id)) 
  
  # retrieve mitocondrial genes 
  mito.genes <- grep(pattern = "^MT-|^mt-", x = rownames(x = scdata@data), value = TRUE) 
  
  # Calculate real mitocondrial percentage (not ratio) and do the plot
  percent.mito <- Matrix::colSums(scdata@raw.data[mito.genes, ])/Matrix::colSums(scdata@raw.data) * 100 
  scdata <- AddMetaData(object = scdata, metadata = percent.mito, col.name = "pct.mito") # adding % mito as metadata
  
  return(scdata)
  
}