scTool_saveAnalysis <- function(markers.obj, seurat.obj, outfold="Analysis_objects/"){
  
  dir.create(path = outfold, showWarnings = F, recursive = T)

  saveRDS(object = markers.obj, file = paste0(outfold,seurat.obj@project.name,"_all_markers_1.rds"))
  saveRDS(object = seurat.obj, file = paste0(outfold,seurat.obj@project.name,".rds"))
  save.image(file = paste0(outfold,seurat.obj@project.name,".RData"))
  savehistory(file = paste0(outfold,seurat.obj@project.name,".Rhistory"))
  
}