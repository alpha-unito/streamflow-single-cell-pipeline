RunMerge <- function(objlist){
  
  library(Seurat)
  
  MergedObj <- MergeSeurat(object1 = objlist[[1]], object2 = objlist[[2]])
  for(i in 3:length(objlist)){
    
    MergedObj <- MergeSeurat(object1 = MergedObj, object2 = objlist[[i]])
  }
  return(MergedObj)
}