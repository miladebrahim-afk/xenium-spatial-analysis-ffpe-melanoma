# Figure S4A: spatial-distance QC by NF1 genotype
suppressPackageStartupMessages({
  library(Seurat); library(dplyr); library(dbscan); library(ggplot2); library(patchwork)
})

xlist_file <- file.path("data","processed","xenium_list_with_meta_niche_published_reference.rds")
out_dir <- file.path("results","spatial_cd8","distance_qc")
dir.create(out_dir, recursive=TRUE, showWarnings=FALSE)

melanoma_labels <- c("Melanoma NC","Melanoma Melanocytic","Melanoma Intermediate",
                     "Melanoma Mesenchymal","Melanoma Proliferative")

norm_nf1 <- function(x) {
  x <- as.character(x)
  y <- ifelse(x %in% c("Mut","NF1Mut","NF1_Mut","NF1 Mut"),"Mut",
       ifelse(x %in% c("WT","NF1WT","NF1_WT","NF1 WT"),"WT",NA))
  if (anyNA(y)) stop("Unexpected NF1 values: ", paste(unique(x[is.na(y)]), collapse=", "))
  y
}

one <- function(x, nm, i) {
  z <- unique(as.character(x[!is.na(x)]))
  if (length(z) != 1) stop("Expected one ", nm, " in core ", i)
  z
}

fmt_p <- function(p) ifelse(p < 0.001, "p < 0.001", paste0("p = ", formatC(p, format="f", digits=3)))

if (!file.exists(xlist_file)) stop("Missing: ", xlist_file)
cat("Loading Xenium tissue-core objects...\n")
xenium_list <- readRDS(xlist_file)
if (length(xenium_list) != 39) stop("Expected 39 cores; found ", length(xenium_list))

rows <- vector("list", length(xenium_list))
cat("\nCalculating per-core QC metrics...\n")

for (i in seq_along(xenium_list)) {
  obj <- xenium_list[[i]]; md <- obj@meta.data
  sample_id <- one(md$patient_id,"patient_id",i)
  nf1 <- one(norm_nf1(md$NF1),"NF1",i)

  coords <- GetTissueCoordinates(obj)
  coords$celltype <- md$celltype[match(coords$cell, rownames(md))]
  coords <- coords %>% filter(!is.na(x), !is.na(y), !is.na(celltype))

  n_total <- nrow(coords)
  area <- (max(coords$x)-min(coords$x)) * (max(coords$y)-min(coords$y))
  if (!is.finite(area) || area <= 0) area <- NA_real_
  density <- ifelse(is.na(area), NA_real_, n_total/area)

  mel <- coords %>% filter(celltype %in% melanoma_labels)
  avg_nbr <- NA_real_
  if (nrow(mel) >= 2) {
    nn <- dbscan::frNN(as.matrix(mel[,c("x","y")]), eps=40, sort=FALSE)
    avg_nbr <- mean(lengths(nn$id)-1, na.rm=TRUE)
  }

  rows[[i]] <- data.frame(sample_index=i, sample_id=sample_id, NF1=nf1,
                          n_total=n_total, analyzed_area=area, total_density=density,
                          avg_melanoma_neighbors=avg_nbr, n_melanoma_cells=nrow(mel))
  cat(sprintf("  core %02d | %-10s | NF1 %-3s | cells=%7d | melanoma=%6d\n",
              i,sample_id,nf1,n_total,nrow(mel)))
  rm(obj,md,coords,mel); gc(verbose=FALSE)
}

core_qc <- bind_rows(rows) %>% mutate(NF1=factor(NF1, levels=c("Mut","WT")))
write.csv(core_qc, file.path(out_dir,"distance_qc_core_metrics.csv"), row.names=FALSE)

metrics <- c("n_total","analyzed_area","total_density","avg_melanoma_neighbors")
labels <- c("Total cells analyzed","Estimated analyzed tissue area","Cellular density",
            "Average number of melanoma neighbors")

stats <- bind_rows(lapply(seq_along(metrics), function(j) {
  m <- metrics[j]; d <- core_qc[!is.na(core_qc[[m]]),]
  p <- suppressWarnings(wilcox.test(d[[m]] ~ d$NF1, exact=FALSE)$p.value)
  data.frame(metric=m,label=labels[j], n_mut=sum(d$NF1=="Mut"), n_wt=sum(d$NF1=="WT"),
             median_mut=median(d[[m]][d$NF1=="Mut"]),
             median_wt=median(d[[m]][d$NF1=="WT"]), wilcox_p=p)
}))
write.csv(stats, file.path(out_dir,"distance_qc_wilcoxon.csv"), row.names=FALSE)

cat("\nFigure S4A Wilcoxon QC:\n"); print(stats, row.names=FALSE)

make_plot <- function(metric, ylab, p) {
  d <- core_qc[!is.na(core_qc[[metric]]),]
  rng <- range(d[[metric]], na.rm=TRUE); ypos <- rng[2] + 0.08*diff(rng)
  ggplot(d, aes(x=NF1,y=.data[[metric]],fill=NF1)) +
    geom_boxplot(width=.55,outlier.shape=NA,alpha=.7) +
    geom_jitter(width=.12,size=2,alpha=.7) +
    annotate("text",x=1.5,y=ypos,label=fmt_p(p),size=3.5) +
    expand_limits(y=ypos) + labs(x=NULL,y=ylab) + theme_classic(base_size=12) +
    theme(legend.position="none",axis.text.x=element_text(face="bold"),
          axis.title.y=element_text(face="bold"))
}

plots <- lapply(seq_along(metrics), function(j)
  make_plot(metrics[j], labels[j], stats$wilcox_p[j]))
p <- plots[[1]] + plots[[2]] + plots[[3]] + plots[[4]] + plot_layout(nrow=1)

ggsave(file.path(out_dir,"FigureS4A_distance_QC_four_panels.pdf"), p,
       width=16,height=4.5,useDingbats=FALSE)
ggsave(file.path(out_dir,"FigureS4A_distance_QC_four_panels.png"), p,
       width=16,height=4.5,dpi=300)

cat("\nDistance-QC summary:\n")
cat("  tissue cores: ",nrow(core_qc),"\n",sep="")
cat("  NF1Mut cores: ",sum(core_qc$NF1=="Mut"),"\n",sep="")
cat("  NF1WT cores: ",sum(core_qc$NF1=="WT"),"\n",sep="")
cat("\nDONE: Figure S4A distance QC complete.\n")
