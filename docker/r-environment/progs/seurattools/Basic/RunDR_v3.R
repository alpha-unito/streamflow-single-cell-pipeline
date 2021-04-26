RunDR.v3 <- function(seurat.obj, id=NULL, umap.f=NULL, npcs=30,
                     tsne.f=NULL, reductiontype="pca", perplexity=30){
  

  # RuntSNE and RunUMAP and plot 2D and 3D
      
  seurat.obj <- RunTSNE(object = seurat.obj, do.fast = T, seed.use = 123, reduction = reductiontype, perplexity = perplexity)
  seurat.obj <- RunUMAP(object = seurat.obj, seed.use = 123, reduction = reductiontype, dims=1:npcs)
  seurat.obj <- RunTSNE(object = seurat.obj, do.fast = T, seed.use = 123, reduction.name = "tSNE3d", 
                        reduction.key = "tSNE3d_", reduction = reductiontype, dim.embed = 3, perplexity = perplexity)
  seurat.obj <- RunUMAP(object = seurat.obj, seed.use = 123, n.components = 3, 
                        reduction.name = "UMAP3d", reduction.key = "UMAP3d_", 
                        reduction = reductiontype, dims=1:npcs)

    # store cell embeddings
    
    tsne.cell_emb_2d <- as.data.frame(seurat.obj@reductions$tsne@cell.embeddings)
    umap.cell_emb_2d <- as.data.frame(seurat.obj@reductions$umap@cell.embeddings)
    tsne.cell_emb_3d <- as.data.frame(seurat.obj@reductions$tSNE3d@cell.embeddings)
    umap.cell_emb_3d <- as.data.frame(seurat.obj@reductions$UMAP3d@cell.embeddings)
    
    # add embeddings to metadata 
    
    seurat.obj <- AddMetaData(object = seurat.obj, metadata = tsne.cell_emb_2d)
    seurat.obj <- AddMetaData(object = seurat.obj, metadata = umap.cell_emb_2d)
    
    seurat.obj <- AddMetaData(object = seurat.obj, metadata = tsne.cell_emb_3d)
    seurat.obj <- AddMetaData(object = seurat.obj, metadata = umap.cell_emb_3d)
    
    return(seurat.obj)
  }
  