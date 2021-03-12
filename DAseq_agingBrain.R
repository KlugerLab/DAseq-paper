### DA-seq on aging mouse brain data
### Original paper: https://www.nature.com/articles/s41593-019-0491-3
### This script reproduces analysis presented in Supp Figure S5,S6

library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)

source("convenience.R")

## Set Python and GPU
python2use <- "/data/henry/henry_env/venv/bin/python"
GPU <- 2

## Set path for FIt-SNE R wrapper
fitsneR <- "~/git/FIt-SNE/fast_tsne.R"


##=============================================##
## Data prep

## Download data

if(!dir.exists("./data/")){
  dir.create("./data/")
}

# Please go to https://singlecell.broadinstitute.org/single_cell/study/SCP263/aging-mouse-brain#/ to download
# processed data.

# Make sure downloaded data is stored under ./data/ directory.



## Load data

# exp mat
data_exp <- read.table(
  "./data/expression_Aging_mouse_brain_portal_data_updated.txt", sep = "\t", header = T, 
  row.names = 1, stringsAsFactors = F
)

data_exp <- Matrix(as.matrix(data_exp), sparse = T)


# meta data
data_meta <- read.table(
  "./data/meta_Aging_mouse_brain_portal_data.txt", sep = "\t", header = T, stringsAsFactors = F
)
data_meta <- data_meta[-1,]
data_meta$cell_type <- gsub("NEUR_mature","mNEUR",data_meta$cell_type)
data_meta$cell_type <- gsub("NEUR_immature","ImmN",data_meta$cell_type)

table(data_meta$cell_type)
table(data_meta$all_cells_by_age)
rownames(data_meta) <- data_meta[,1]

celltype2label <- c(
  "1-OPC","2-OLG","3-OEG","4-NSC","5-ARP","6-ASC","7-NRP","8-ImmN","9-mNEUR","10-NendC",
  "11-EPC","12-HypEPC","13-TNC","14-CPC","15-EC","16-PC","17-VSMC","18-Hb_VC","19-VLMC","20-ABC",
  "21-MG","22-MNC","23-MAC","24-DC","25-NEUT"
)
names(celltype2label) <- sapply(celltype2label, function(x) unlist(strsplit(x,"-"))[2])



## Seurat

# create object
data_S <- CreateSeuratObject(
  counts = data_exp, names.field = 6, names.delim = "_"
)
table(data_S@meta.data$orig.ident)


# add metadata
data_S@meta.data$cell_type <- data_meta[colnames(data_S),"cell_type"]
data_S@meta.data$age <- data_meta[colnames(data_S),"all_cells_by_age"]
data_S@meta.data$condition <- data_S@meta.data$age
data_S@meta.data$condition[data_S@meta.data$age == "2-3mo"] <- "young"
data_S@meta.data$condition[data_S@meta.data$age == "21-22mo"] <- "old"
table(data_S@meta.data$condition)

data_S@meta.data$cell_type_label <- factor(celltype2label[data_S@meta.data$cell_type], levels = celltype2label)
data_S@meta.data$cell_type_num <- as.numeric(data_S@meta.data$cell_type_label)


# analysis
data_S <- ScaleData(data_S)
data_S <- FindVariableFeatures(
  data_S, selection.method = "mvp", mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf)
)
data_S <- RunPCA(
  data_S, npcs = 20, verbose = F
)
data_S <- runFItSNE(
  data_S, dims.use = 1:20, seed.use = 3, fast.R.path = fitsneR, 
  ann_not_vptree = FALSE, nthreads = 12
)

TSNEPlot(data_S, group.by = "age", pt.size = 0.1)
TSNEPlot(data_S, group.by = "cell_type_label", label = T, pt.size = 0.1)





##=============================================##
## DA-Seq


## labels

labels_old <- unique(as.character(data_S@meta.data$orig.ident[data_S@meta.data$age == "21-22mo"]))
labels_young <- unique(as.character(data_S@meta.data$orig.ident[data_S@meta.data$age == "2-3mo"]))



## DA cells

da_cells <- getDAcells(
  X = data_S@reductions$pca@cell.embeddings, 
  cell.labels = as.character(data_S@meta.data$orig.ident), 
  labels.1 = labels_young, labels.2 = labels_old, 
  k.vector = round(seq(50,1000,length.out = 10)), 
  plot.embedding = data_S@reductions$tsne@cell.embeddings
)

da_cells <- updateDAcells(
  da_cells, 
  plot.embedding = data_S@reductions$tsne@cell.embeddings, size = 0.01
)

da_cells$pred.plot
da_cells$da.cells.plot



## DA regions

da_regions <- getDAregion(
  X = data_S@reductions$pca@cell.embeddings, da.cells = da_cells, 
  cell.labels = as.character(data_S$orig.ident), 
  labels.1 = labels_young, labels.2 = labels_old, 
  resolution = 0.05, min.cell = 50, 
  plot.embedding = data_S@reductions$tsne@cell.embeddings, size = 0.01, do.label = T
)

da_regions$da.region.plot
da_regions$DA.stat

n_da <- length(unique(da_regions$da.region.label)) - 1

xlims <- c(min(data_S[["tsne"]]@cell.embeddings[data_S$cell_type == "MG",1]), 
           max(data_S[["tsne"]]@cell.embeddings[data_S$cell_type == "MG",1]))
ylims <- c(min(data_S[["tsne"]]@cell.embeddings[data_S$cell_type == "MG",2]), 
           max(data_S[["tsne"]]@cell.embeddings[data_S$cell_type == "MG",2]))
gridExtra::grid.arrange(
  TSNEPlot(data_S, cells = which(da_regions$da.region.label == 1&data_S$cell_type == "MG"), 
           group.by = "age", pt.size = 0.01) + xlim(xlims) + ylim(ylims),
  TSNEPlot(data_S, cells = which(da_regions$da.region.label != 1&data_S$cell_type == "MG"), 
           group.by = "age", pt.size = 0.01) + xlim(xlims) + ylim(ylims), 
  nrow=1
)



## DA markers

STG_markers <- STGmarkerFinder(
  X = as.matrix(data_S@assays$RNA@data), 
  da.regions = da_regions, 
  lambda = 1.5, n.runs = 5, return.model = T, 
  python.use = python2use, GPU = GPU
)



## Compare DA2 and DA9 (within cluster 11/cell type MG)

data_S <- addDAslot(data_S, da.regions = da_regions, set.ident = T)
Seurat_local_marker <- FindMarkers(
  data_S, ident.1 = 2, ident.2 = 9
)





##=============================================##
## Young mice analysis

young_S <- subset(data_S, cells = which(data_S@meta.data$condition == "young"))
young_S <- ScaleData(young_S)
young_S <- FindVariableFeatures(
  young_S, selection.method = "mvp", mean.cutoff = c(0.0125,3), dispersion.cutoff = c(0.5,Inf)
)
young_S <- RunPCA(
  young_S, npcs = 20, verbose = F
)

young_idx1 <- sample(c(1:8), size = 4)
young_idx2 <- setdiff(c(1:8),young_idx1)
young_S@meta.data$condition[young_S@meta.data$orig.ident%in%labels_young[young_idx1]] <- "0"
young_S@meta.data$condition[young_S@meta.data$orig.ident%in%labels_young[young_idx2]] <- "1"
TSNEPlot(young_S, group.by = "condition")
da_cells_young <- getDAcells(
  X = young_S@reductions$pca@cell.embeddings,
  cell.labels = as.character(young_S@meta.data$orig.ident),
  labels.1 = labels_young[young_idx1], labels.2 = labels_young[young_idx2],
  k.vector = seq(50,500,50), pred.thres = c(0,0), 
  plot.embedding = young_S@reductions$tsne@cell.embeddings
)
da_cells_young <- updateDAcells(
  da_cells_young, pred.thres = c(0,0), 
  plot.embedding = young_S@reductions$tsne@cell.embeddings, size = 0.01
)

da_regions_young <- getDAregion(
  X = young_S@reductions$pca@cell.embeddings, da.cells = da_cells_young, 
  cell.labels = as.character(young_S@meta.data$orig.ident),
  labels.1 = labels_young[young_idx1], labels.2 = labels_young[young_idx2], min.cell = 50, 
  plot.embedding = young_S@reductions$tsne@cell.embeddings
)





##=============================================##
## Generate plots

library(scales)
da_cols <- hue_pal()(n_da)
da_order <- order(da_regions$da.region.label)
idx_local <- intersect(which(da_regions$da.region.label %in% c(0,2,9)),which(data_S@meta.data$cell_type == "MG"))
da_order_local <- order(da_regions$da.region.label[idx_local])

tsne_embedding <- data_S@reductions$tsne@cell.embeddings

## TSNE plots

gg1 <- plotCellLabel(
  tsne_embedding, label = data_S@meta.data$condition, size = 0.01, do.label = F
) + theme_tsne
ggsave(gg1, filename = "figs/agingBrain_a.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg1, legend.position = "right"), 
       filename = "figs/agingBrain_a_legend.pdf", width = 0.5, height = 0.3, dpi = 1200)


gg2 <- plotCellLabel(
  tsne_embedding, label = factor(data_S@meta.data$cell_type_num), 
  size = 0.01, label.size = 2
) + scale_color_hue(labels = levels(data_S@meta.data$cell_type_label)) + theme_tsne
ggsave(gg2, filename = "figs/agingBrain_b.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg2), filename = "figs/agingBrain_b_legend.pdf", width = 1.6, height = 1.6, dpi = 1200)


gg3 <- da_cells$pred.plot + theme_tsne
ggsave(gg3, filename = "figs/agingBrain_c.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg3, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm"), legend.title = element_blank()), 
       filename = "figs/agingBrain_c_legend.pdf", height = 30, width = 15, units = "mm", dpi = 1200)


gg4 <- plotCellLabel(
  tsne_embedding[da_order,], label = as.factor(da_regions$da.region.label[da_order]), 
  size = 0.01, label.size = 2, label.plot = as.character(c(1:n_da))
) + scale_color_manual(
  values = c("gray", da_cols), breaks = c(1:n_da), labels = paste0("DA",c(1:n_da))
) + theme_tsne
ggsave(gg4, filename = "figs/agingBrain_d.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg4), filename = "figs/agingBrain_d_legend.pdf", width = 0.5, height = 2.1, dpi = 1200)
gg4sub <- plotCellLabel(
  tsne_embedding[idx_local,][da_order_local,], 
  factor(da_regions$da.region.label[idx_local][da_order_local]),
  size = 0.1, do.label = F, cell.col = c("gray", da_cols[c(2,9)])
) + theme_tsne
ggsave(gg4sub, filename = "figs/agingBrain_d_sub.png", width = 25, height = 25, units = "mm", dpi = 1200)



## Rand plot

gg1 <- ggplot() +
  geom_point(data = data.frame(
    order = seq(1,ncol(data_S),length.out = length(unlist(da_cells$rand.pred))), random = sort(unlist(da_cells$rand.pred))
  ), aes(order, random), col = "gray", size = 0.1, alpha = 0.5) +
  geom_point(data = data.frame(
    order = c(1:ncol(data_S)), da = sort(da_cells$da.pred)
  ), aes(order,da), col = "black", size = 0.1, alpha = 0.75) +
  geom_hline(yintercept = min(unlist(da_cells$rand.pred)), size = 0.1) +
  geom_hline(yintercept = max(unlist(da_cells$rand.pred)), size = 0.1) +
  ylim(-1,1) + theme_tsne + theme(axis.text.y = element_text(size=6), axis.ticks.y = element_line(size = 0.2))
ggsave(gg1, filename = "figs/agingBrain_rand.png", width = 50, height = 50, units = "mm", dpi = 1200)



## Dot plot

# marker genes to plot
marker_genes <- list(
  "1" = "Tmem119",
  "2" = "Sox11",
  "3" = "Pdgfra",
  "4" = "Cdk1",
  "5" = "Rassf10",
  "6" = "Cdca3"
)


# add STG info
STG.marker.info <- do.call(rbind, lapply(STG_markers$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(unlist(marker_genes))))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(unlist(marker_genes))
rownames(STG.marker.info) <- c(1:(n_da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]
STG.marker.info.m$value <- 7.5 * STG.marker.info.m$value


# generate dot plot
gg5 <- DotPlot(
  data_S, features = unlist(marker_genes), cols = c("gray","blue"), group.by = "da"
) + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value)) + theme_dot + RotatedAxis()
ggsave(gg5, filename = "figs/agingBrain_e.pdf", width = 50, height = 40, units = "mm", dpi = 1200)
ggsave(
  g_legend(gg5, legend.key.height = unit(0.15,"cm"), legend.key.width = unit(0.2,"cm")), 
  filename = "figs/agingBrain_e_legend.pdf", height = 40, width = 30, units = "mm", dpi = 1200
)



## Feature plots

# DA1, MG
gg6 <- plotCellScore(
  tsne_embedding[STG_local_marker$model$cells,], STG_local_marker$model$pred, 
  cell.col = c("gray","blue"), size = 0.5
) + theme_tsne
ggsave(gg6, filename = "figs/agingBrain_f.png", width = 50, height = 50, units = "mm", dpi = 1200)

# DA5, Rassf10
sgg1 <- list(
  plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions$da.region.label[da_order]), size = 0.1, do.label = F, 
    cell.col = c("gray","gray","gray","gray","gray",da_cols[5],"gray")
  ) + ggtitle("DA5") + theme_tsne,
  plotCellScore(
    tsne_embedding, data_S@assays$RNA@data["Rassf10",], cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("Rassf10") + theme_tsne
)
ggsave(
  plot_grid(plotlist = sgg1, nrow = 1), 
  filename = "figs/agingBrain_s_a.png", height = 45, width = 80, units = "mm", dpi = 1200
)

# DA1, MG
sgg2 <- list(
  plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions$da.region.label[da_order]), size = 0.1, do.label = F, 
    cell.col = c("gray",da_cols[1],"gray","gray","gray","gray","gray")
  ) + ggtitle("DA1") + theme_tsne,
  plotCellScore(
    tsne_embedding, data_S@assays$RNA@data["Tmem119",], cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("Tmem119") + theme_tsne
)
ggsave(
  plot_grid(plotlist = sgg2, nrow = 1), 
  filename = "figs/agingBrain_s_b.png", height = 45, width = 80, units = "mm", dpi = 1200
)
sgg3 <- list(
  plotCellLabel(
    tsne_embedding[STG_local_marker$model$cells,][da_order_local,], 
    factor(da_regions$da.region.label[match(STG_local_marker$model$cells,colnames(data_S))][da_order_local]),
    cell.col = c("gray",da_cols[1]), size = 0.5, do.label = F
  ) + ggtitle("DA1") + theme_tsne,
  plotCellScore(
    tsne_embedding[STG_local_marker$model$cells,], 
    data_S@assays$RNA@data["Lyz2",STG_local_marker$model$cells], 
    cell.col = c("gray","blue"), size = 0.5
  ) + ggtitle("Lyz2") + theme_tsne,
  plotCellScore(
    tsne_embedding[STG_local_marker$model$cells,], 
    data_S@assays$RNA@data["Aldoc",STG_local_marker$model$cells], 
    cell.col = c("gray","blue"), size = 0.5
  ) + ggtitle("Aldoc") + theme_tsne
)
ggsave(
  plot_grid(plotlist = sgg3, nrow = 1), filename = "figs/agingBrain_s_c.png", 
  width = 120, height = 45, units = "mm", dpi = 1200
)

ggsave(g_legend(
  gg6, 
  legend.text = element_blank(), legend.title = element_blank(), 
  legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm")
), filename = "figs/agingBrain_f_legend.pdf", height = 30, width = 10, units = "mm", dpi = 1200)



## Young mice plots

gg1 <- plotCellLabel(
  young_S@reductions$tsne@cell.embeddings, label = young_S@meta.data$orig.ident, size = 0.01, do.label = F
) + theme_tsne
ggsave(gg1, filename = "figs/agingBrain_young_1.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg1, legend.position = "right"), 
       filename = "figs/agingBrain_a_legend.pdf", width = 0.5, height = 0.3, dpi = 1200)

gg1 <- ggplot() +
  geom_point(data = data.frame(
    order = seq(1,ncol(young_S),length.out = length(unlist(da_cells_young$rand.pred))), random = sort(unlist(da_cells_young$rand.pred))
  ), aes(order, random), col = "gray", size = 0.1, alpha = 0.5) +
  geom_point(data = data.frame(
    order = c(1:ncol(young_S)), da = sort(da_cells_young$da.pred)
  ), aes(order,da), col = "black", size = 0.1, alpha = 0.75) +
  geom_hline(yintercept = min(unlist(da_cells_young$rand.pred)), size = 0.1) +
  geom_hline(yintercept = max(unlist(da_cells_young$rand.pred)), size = 0.1) +
  ylim(-1,1) + theme_tsne + theme(axis.text.y = element_text(size=6), axis.ticks.y = element_line(size = 0.2))
ggsave(gg1, filename = "figs/agingBrain_young_2.png", width = 50, height = 50, units = "mm", dpi = 1200)

gg1 <- da_cells_young$da.cells.plot + theme_tsne
ggsave(gg1, filename = "figs/agingBrain_young_3.png", width = 50, height = 50, units = "mm", dpi = 1200)



