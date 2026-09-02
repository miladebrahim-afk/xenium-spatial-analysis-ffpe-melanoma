# Additional spatial focus panels used for Figure 3B/E.
source(file.path('config','project_paths.R')); source(file.path('R','project_io.R'))
suppressPackageStartupMessages({library(Seurat);library(ggplot2)})
xlist <- load_xenium_list(); out <- ensure_dir(file.path(FIGURES_DIR,'spatial_focus'))
for (i in seq_along(xlist)) {
  obj <- xlist[[i]]; pid <- paste(unique(as.character(obj$patient_id)),collapse='_')
  obj$focus_group <- 'Other'
  obj$focus_group[grepl('Melanoma',obj$celltype,ignore.case=TRUE)] <- 'Melanoma'
  obj$focus_group[grepl('CD8',obj$celltype,ignore.case=TRUE)] <- 'CD8'
  obj$focus_group[grepl('CAF',obj$celltype,ignore.case=TRUE)] <- 'CAF'
  p <- ImageDimPlot(obj,group.by='focus_group',size=.8)+ggtitle(pid)
  ggsave(file.path(out,paste0('Melanoma_CD8_CAF_',gsub('[^A-Za-z0-9_-]','_',pid),'.pdf')),p,width=7,height=6)
}
