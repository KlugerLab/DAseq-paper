### DA-seq on COVID-19 data
### Original paper: https://www.nature.com/articles/s41587-020-0602-4
### This script reproduces analysis presented in Figure 2

library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)

source("convenience.R")

## Set Python and GPU
python2use <- "/data/henry/henry_env/venv/bin/python"
GPU <- 3


##=============================================##
## Data

## Load data

main_S <- readRDS(url("https://ndownloader.figshare.com/files/22927382"))
DefaultAssay(main_S) <- "integrated"
head(main_S@meta.data)
table(main_S@meta.data$severity)

# check commanes
names(main_S@commands)



## Prepare object according to commands

# analysis
main_S <- ScaleData(main_S)
main_S <- RunPCA(
  main_S, npcs = 90, verbose = F
)
main_S <- runFItSNE(
  main_S, dims.use = 1:90, seed.use = 3, 
  ann_not_vptree = FALSE, nthreads = 12
)
TSNEPlot(main_S, group.by = "severity")
TSNEPlot(main_S, group.by = "celltype", label = T)



## Separate immune cells from patient samples

# get immune cells
sort(unique(main_S@meta.data$celltype))
epi_type <- c("Basal","Ciliated","Ciliated-diff","FOXN4","Ionocyte","IRC",
              "Secretory","Secretory-diff","Squamous","outliers_epithelial","unknown_epithelial")
immune_type <- setdiff(sort(unique(main_S@meta.data$celltype)), epi_type)

immune_S <- subset(x = main_S, cells = which(main_S@meta.data$celltype %in% immune_type))

# remove control cells
immune_S <- subset(immune_S, cells = which(immune_S@meta.data$severity != "control"))
immune_S@meta.data$severity <- factor(as.character(immune_S@meta.data$severity))

TSNEPlot(immune_S, group.by = "severity")
TSNEPlot(immune_S, group.by = "celltype", label = T)


# assign a number label for each cell type
immune_S@meta.data$celltype_num <- as.numeric(factor(immune_S@meta.data$celltype))





##=============================================##
## DA-seq on immune cells

## Get sample labels
table(main_S@meta.data$patient)
table(main_S@meta.data$sample)
moderate_labels <- unique(main_S@meta.data$sample[main_S@meta.data$severity == "moderate"])
critical_labels <- unique(main_S@meta.data$sample[main_S@meta.data$severity == "critical"])



## DA cells

da_cells_immune <- getDAcells(
  X = immune_S@reductions$pca@cell.embeddings, 
  cell.labels = immune_S@meta.data$sample, 
  labels.1 = moderate_labels, labels.2 = critical_labels, 
  k.vector = seq(200,3200,200), n.runs = 5, 
  plot.embedding = immune_S@reductions$tsne@cell.embeddings, 
)

da_cells_immune <- updateDAcells(
  da_cells_immune, pred.thres = c(0.05,0.95),
  plot.embedding = immune_S@reductions$tsne@cell.embeddings, size = 0.01
)
da_cells_immune$pred.plot
da_cells_immune$da.cells.plot



## DA regions

da_regions_immune <- getDAregion(
  X = immune_S@reductions$pca@cell.embeddings, da.cells = da_cells_immune, 
  cell.labels = immune_S@meta.data$sample, 
  labels.1 = moderate_labels, labels.2 = critical_labels, 
  resolution = 0.01, 
  plot.embedding = immune_S@reductions$tsne@cell.embeddings, size = 0.1
)
da_regions_immune$da.region.plot
da_regions_immune$DA.stat

n_da_immune <- length(unique(da_regions_immune$da.region.label)) - 1



## Local markers

# set cluster overlap for each DA
da_clusters_immune <- c(
  "1" = "Neu", "2" = "nrMa", "4" = "CTL", "5" = "Neu"
)

# Seurat negbinom to find local markers
immune_S <- addDAslot(immune_S, da.regions = da_regions_immune)
local_markers_immune <- list()
for(i in 4:n_da_immune){
  if(!as.character(i) %in% names(da_clusters_immune)){next()}
  local_markers_immune[[as.character(i)]] <- SeuratLocalMarkers(
    object = immune_S, da.region.to.run = i, cell.label.slot = "celltype", cell.label.to.run = da_clusters_immune[as.character(i)], 
    assay = "RNA", test.use = "negbinom", min.diff.pct = 0.09, only.pos = T
  )
  local_markers_immune[[as.character(i)]]$pct.diff <- local_markers_immune[[as.character(i)]]$pct.1 - 
    local_markers_immune[[as.character(i)]]$pct.2
}





##=============================================##
## Data from Liao et al.
## https://www.nature.com/articles/s41591-020-0901-9

## Download
liao_S <- readRDS(url("http://cells.ucsc.edu/covid19-balf/nCoV.rds"))


## Set clusters
epi_clusters <- c(13,16,25,28,31)

macro_clusters <- as.character(c(0,1,2,3,4,5,7,8,10,11,12,18,21,22,23,26))
t_clusters <- as.character(c(6,9,14))
cd8t_clusters <- as.character(c(9,14))


## Calculate module score
liao_immune_S <- subset(liao_S, cells = which(!liao_S@meta.data$cluster %in% epi_clusters & liao_immune_S$group != "HC"))
liao_immune_S@meta.data$severity <- gsub("O", "moderate", liao_immune_S@meta.data$group1)
liao_immune_S@meta.data$severity <- gsub("S/C", "critical", liao_immune_S@meta.data$severity)

da_gene_modules <- lapply(local_markers_immune, FUN = function(x){
  rownames(x)[1:min(100,nrow(x))]
})
names(da_gene_modules) <- paste0("DA", names(da_gene_modules))

liao_immune_S <- AddModuleScore(
  object = liao_immune_S, features = da_gene_modules, assay = "RNA", name = names(da_gene_modules)
)
for(i in 1:4){
  colnames(liao_immune_S@meta.data)[grep(names(da_gene_modules)[i],colnames(liao_immune_S@meta.data))] <- 
    names(da_gene_modules)[i]
}





##=============================================##
## Generate plots

library(scales)
da_cols <- hue_pal()(n_da_immune)
da_order <- order(da_regions_immune$da.region.label)

tsne_embedding <- immune_S@reductions$tsne@cell.embeddings

## TSNE plots

gg1 <- plotCellLabel(
  tsne_embedding, label = immune_S@meta.data$severity, size = 0.01, do.label = F
) + theme_tsne
ggsave(gg1, filename = "figs/covidChua_a.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg1, legend.position = "right"), 
       filename = "figs/covidChua_a_legend.pdf", width = 0.7, height = 0.3, dpi = 1200)


gg2 <- plotCellLabel(
  tsne_embedding, label = factor(immune_S@meta.data$celltype_num), 
  size = 0.01, label.size = 2
) + scale_color_hue(labels = paste(sort(unique(immune_S@meta.data$celltype_num)), immune_type, sep = "-")) + theme_tsne
ggsave(gg2, filename = "figs/covidChua_b.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg2), filename = "figs/covidChua_b_legend.pdf", width = 1, height = 1.6, dpi = 1200)


gg3 <- da_cells_immune$pred.plot + theme_tsne
ggsave(gg3, filename = "figs/covidChua_c.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg3, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm")), 
       filename = "figs/covidChua_c_legend.pdf", height = 30, width = 15, units = "mm", dpi = 1200)


gg4 <- plotCellLabel(
  tsne_embedding[da_order,], label = as.character(da_regions_immune$da.region.label[da_order]), 
  size = 0.01, label.size = 2, label.plot = as.character(c(1:n_da_immune))
) + scale_color_manual(
  values = c("gray", da_cols), breaks = c(1:n_da_immune), labels = paste0("DA",c(1:n_da_immune))
) + theme_tsne
ggsave(gg4, filename = "figs/covidChua_d.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg4), filename = "figs/covidChua_d_legend.pdf", width = 0.5, height = 0.8, dpi = 1200)



## Dot plot

marker_genes <- list(
  "1" = c("CXCR4","CD63","CD48","IRAK3"),
  "2" = c("RGL1","MAFB","MAF","CD163"),
  "4" = c("IFNG","CCR6","TOX2","IL26"),
  "5" = c("IL1RN","SOCS3","PTGS2","IL1B")
)

DefaultAssay(immune_S) <- "RNA"
gg5 <- DotPlot(
  immune_S, features = unlist(marker_genes), cols = c("gray","blue"), group.by = "da"
) + theme_dot + RotatedAxis() + guides(
  color = guide_colorbar(title = "Average Expression", order = 2), 
  size = guide_legend(title = "Percent Expressed", order = 1)
)
ggsave(gg5, filename = "figs/covidChua_e.pdf", width = 90, height = 40, units = "mm", dpi = 1200)
ggsave(
  g_legend(gg5, legend.key.height = unit(0.15,"cm"), legend.key.width = unit(0.2,"cm"), legend.spacing = unit(0.25, 'cm')), 
  filename = "figs/covidChua_e_legend.pdf", height = 50, width = 30, units = "mm", dpi = 1200
)
DefaultAssay(immune_S) <- "integrated"



## Feature plots

# DA1
sgg1 <- c(
  list(plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions_immune$da.region.label[da_order]), size = 0.01, do.label = F, 
    cell.col = c("gray",da_cols[1],"gray","gray","gray","gray")
  ) + ggtitle("DA1") + theme_tsne), 
  lapply(c("CD48","CD63","CXCR4"), FUN = function(x){
    plotCellScore(
      tsne_embedding, immune_S@assays$RNA@data[x,], cell.col = c("gray","blue"), size = 0.01
    ) + ggtitle(x) + theme_tsne
  })
)
ggsave(
  plot_grid(plotlist = sgg1, nrow = 1), 
  filename = "figs/covidChua_s_a.png", height = 45, width = 160, units = "mm", dpi = 1200
)

# DA5
sgg2 <- c(
  list(plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions_immune$da.region.label[da_order]), size = 0.01, do.label = F, 
    cell.col = c("gray","gray","gray","gray","gray",da_cols[5])
  ) + ggtitle("DA5") + theme_tsne), 
  lapply(c("IL1RN","SOCS3","PTGS2"), FUN = function(x){
    plotCellScore(
      tsne_embedding, immune_S@assays$RNA@data[x,], cell.col = c("gray","blue"), size = 0.01
    ) + ggtitle(x) + theme_tsne
  })
)
ggsave(
  plot_grid(plotlist = sgg2, nrow = 1), 
  filename = "figs/covidChua_s_b.png", height = 45, width = 160, units = "mm", dpi = 1200
)

# DA2
sgg3 <- c(
  list(plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions_immune$da.region.label[da_order]), size = 0.01, do.label = F, 
    cell.col = c("gray","gray",da_cols[2],"gray","gray","gray")
  ) + ggtitle("DA2") + theme_tsne), 
  lapply(c("RGL1","MAFB"), FUN = function(x){
    plotCellScore(
      tsne_embedding, immune_S@assays$RNA@data[x,], cell.col = c("gray","blue"), size = 0.01
    ) + ggtitle(x) + theme_tsne
  })
)
ggsave(
  plot_grid(plotlist = sgg3, nrow = 1), 
  filename = "figs/covidChua_s_c.png", height = 45, width = 120, units = "mm", dpi = 1200
)

# DA4
sgg4 <- c(
  list(plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions_immune$da.region.label[da_order]), size = 0.01, do.label = F, 
    cell.col = c("gray","gray","gray","gray",da_cols[4],"gray")
  ) + ggtitle("DA4") + theme_tsne), 
  lapply("IFNG", FUN = function(x){
    plotCellScore(
      tsne_embedding, immune_S@assays$RNA@data[x,], cell.col = c("gray","blue"), size = 0.01
    ) + ggtitle(x) + theme_tsne
  })
)
ggsave(
  plot_grid(plotlist = sgg4, nrow = 1), 
  filename = "figs/covidChua_s_d.png", height = 45, width = 80, units = "mm", dpi = 1200
)



## Violin plots of gene module scores for Liao's data

gg6 <- list(
  # neutrophil, DA1&DA5, C15
  VlnPlot(
    subset(liao_immune_S, cells = which(liao_immune_S$cluster == "15")),
    features = "DA1", group.by = "severity", pt.size = 0
  ) + theme_tsne,
  VlnPlot(
    subset(liao_immune_S, cells = which(liao_immune_S$cluster == "15")),
    features = "DA5", group.by = "severity", pt.size = 0
  ) + theme_tsne,
  # macrophage, DA2
  VlnPlot(
    subset(liao_immune_S, cells = which(liao_immune_S$cluster %in% macro_clusters)), 
    features = "DA2", group.by = "severity", pt.size = 0
  ) + theme_tsne,
  # CTL, DA4
  VlnPlot(
    subset(liao_immune_S, cells = which(liao_immune_S$cluster %in% cd8t_clusters)), 
    features = "DA4", group.by = "severity", pt.size = 0
  ) + theme_tsne
)
ggsave(
  plot_grid(plotlist = gg6, nrow = 2), 
  filename = "figs/covidChua_s_e.pdf", height = 45, width = 45, units = "mm", dpi = 1200
)




