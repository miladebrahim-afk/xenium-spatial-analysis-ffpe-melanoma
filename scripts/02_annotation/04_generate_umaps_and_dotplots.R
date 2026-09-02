# Figures 1B, 1C, 1E: final UMAPs and marker DotPlot.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))
source(file.path('R','constants.R'))
source(file.path('R','analysis_helpers.R'))
suppressPackageStartupMessages({library(Seurat);library(ggplot2);library(grid)})

obj <- load_merged_object()
require_metadata(obj,c('celltype','major_celltype','NF1'))
out <- ensure_dir(file.path(FIGURES_DIR,'figure1_annotation'))

# Prefer the exact UMAP retained in merged_obj_filtered. If it is absent, rebuild
# only the embedding from the best retained reduction; this fallback will not be
# numerically identical to a historically saved UMAP.
if (!'umap' %in% Reductions(obj)) {
  red <- if ('harmony' %in% Reductions(obj)) 'harmony' else if ('pca' %in% Reductions(obj)) 'pca' else stop('No UMAP/Harmony/PCA reduction found.')
  nd <- min(N_PCS,ncol(Embeddings(obj,red)))
  set.seed(42)
  obj <- RunUMAP(obj,reduction=red,dims=seq_len(nd),seed.use=42,verbose=FALSE)
  saveRDS(obj,ANNOTATED_OBJECT_FILE)
}

add_mini_axes <- function(p,object,reduction='umap',axis_length=1) {
  emb <- Embeddings(object,reduction); xmin<-min(emb[,1]); ymin<-min(emb[,2])
  p + theme(axis.line=element_blank(),axis.text=element_blank(),axis.ticks=element_blank(),axis.title=element_blank()) +
    annotate('segment',x=xmin,xend=xmin+axis_length,y=ymin,yend=ymin,linewidth=.5,arrow=arrow(length=unit(.18,'cm'),type='closed')) +
    annotate('segment',x=xmin,xend=xmin,y=ymin,yend=ymin+axis_length,linewidth=.5,arrow=arrow(length=unit(.18,'cm'),type='closed'))
}

p1 <- DimPlot(obj,reduction='umap',group.by='major_celltype',label=TRUE,label.size=4,repel=TRUE,shuffle=TRUE) + theme(legend.position='none')
p1 <- add_mini_axes(p1,obj)
save_gg(p1,file.path(out,'Fig1B_major_celltype_UMAP.pdf'),6,6)

p2 <- DimPlot(obj,reduction='umap',group.by='celltype',label=TRUE,label.size=3,repel=TRUE,shuffle=TRUE) + theme(legend.position='none')
p2 <- add_mini_axes(p2,obj)
save_gg(p2,file.path(out,'Fig1E_celltype_UMAP.pdf'),7,7)

p3 <- DimPlot(obj,reduction='umap',group.by='celltype',split.by='NF1',label=TRUE,label.size=3,repel=TRUE) + theme(legend.position='none')
save_gg(p3,file.path(out,'UMAP_celltypes_split_by_NF1.pdf'),12,6)

obj <- normalize_rna_if_needed(obj)
markers <- intersect(CANONICAL_MARKERS,rownames(obj[['RNA']]))
pdot <- DotPlot(obj,features=markers,group.by='celltype',assay='RNA',dot.scale=6) +
  RotatedAxis() + theme(axis.text.x=element_text(angle=45,hjust=1))
save_gg(pdot,file.path(out,'Fig1C_celltype_marker_dotplot.pdf'),15,8)
