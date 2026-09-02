# Patient-level major and detailed cell-type composition (supporting Fig. 1 / Fig. S2).
# Figure 3C CD8/Treg values are generated separately in
# scripts/05_spatial_distance/03_tcell_abundance_and_ratio.R.
source(file.path('config','project_paths.R'))
source(file.path('R','project_io.R'))
source(file.path('R','analysis_helpers.R'))
suppressPackageStartupMessages({library(dplyr);library(ggplot2);library(ggpubr)})
obj <- load_merged_object(); md <- obj@meta.data
out <- ensure_dir(file.path(RESULTS_DIR,'annotation')); figout <- ensure_dir(file.path(FIGURES_DIR,'composition'))

celltype_comp <- md |> count(patient_id,NF1,celltype,name='n') |> group_by(patient_id) |> mutate(fraction=n/sum(n)) |> ungroup()
major_comp <- md |> count(patient_id,NF1,major_celltype,name='n') |> group_by(patient_id) |> mutate(fraction=n/sum(n)) |> ungroup()
write.csv(celltype_comp,file.path(out,'celltype_fraction_per_patient.csv'),row.names=FALSE)
write.csv(major_comp,file.path(out,'major_celltype_fraction_per_patient.csv'),row.names=FALSE)

p1 <- ggplot(major_comp,aes(NF1,fraction,fill=NF1))+geom_boxplot(outlier.shape=NA,alpha=.8)+geom_jitter(width=.15,size=1.5,alpha=.7)+
  stat_compare_means(method='wilcox.test',label='p.format')+facet_wrap(~major_celltype,scales='free_y')+theme_bw(base_size=11)+theme(legend.position='none')+labs(x='NF1 status',y='Fraction of cells per patient')
save_gg(p1,file.path(figout,'major_celltype_composition_NF1.pdf'),10,6)

p2 <- ggplot(celltype_comp,aes(NF1,fraction,fill=NF1))+geom_boxplot(outlier.shape=NA,alpha=.8)+geom_jitter(width=.15,size=1.1,alpha=.7)+
  stat_compare_means(method='wilcox.test',label='p.format')+facet_wrap(~celltype,scales='free_y',ncol=5)+theme_bw(base_size=9)+theme(legend.position='none')+labs(x='NF1 status',y='Fraction of cells per patient')
save_gg(p2,file.path(figout,'celltype_composition_NF1.pdf'),13,9)
