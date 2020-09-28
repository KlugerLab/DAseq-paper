### DA-seq on melanoma immune cells data
### Original paper: https://www.sciencedirect.com/science/article/pii/S0092867418313941
### This script reproduces analysis presented in Figure 3

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

## Set path for FIt-SNE R wrapper
fitsneR <- "~/git/FIt-SNE/fast_tsne.R"


##=============================================##
## Data prep

## Load data

if(!dir.exists("./data/")){
  dir.create("./data/")
}

# Expression matrix
download.file(
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  "./data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz"
)

data_exp <- read.table(
  "./data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, row.names = 1, stringsAsFactors = F, skip = 2
)
data_exp <- data_exp[,-16292]

data_colInfo <- read.table(
  "./data/GSE120575_Sade_Feldman_melanoma_single_cells_TPM_GEO.txt.gz",
  sep = "\t", header = F, stringsAsFactors = F, nrows = 2
)
data_colInfo <- data_colInfo[,-1]

colnames(data_exp) <- data_colInfo[1,]
data_exp <- Matrix(as.matrix(data_exp), sparse = T)


# Patient info
download.file(
  "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE120nnn/GSE120575/suppl/GSE120575_patient_ID_single_cells.txt.gz",
  "./data/GSE120575_patient_ID_single_cells.txt.gz"
)
patient_info <- read.table(
  "./data/GSE120575_patient_ID_single_cells.txt.gz", 
  sep = "\t", header = T, stringsAsFactors = F, skip = 19, nrows = 16291
)
mean(colnames(data_exp) == patient_info$title)
rownames(patient_info) <- patient_info$title


# Cluster info
download.file(
  "https://raw.githubusercontent.com/KlugerLab/DAseq-paper/master/data/melanoma_cluster_info",
  "./data/melanoma_cluster_info"
)
cluster_info <- read.table(
  "./data/melanoma_cluster_info", sep = "\t", header = T, stringsAsFactors = F
)
rownames(cluster_info) <- cluster_info$Cell.Name



## Seurat

data_S <- CreateSeuratObject(
  counts = data_exp, project = "melanoma.immune"
)


# set metadata for each cell
data_S@meta.data$condition <- patient_info[colnames(data_S), "characteristics..response"]
data_S@meta.data$lesion <- patient_info[colnames(data_S), 
                                        "characteristics..patinet.ID..Pre.baseline..Post..on.treatment."]
data_S@meta.data$cluster <- paste0("G", cluster_info$Cluster.number)
data_S@meta.data$cluster <- factor(data_S@meta.data$cluster, levels = paste0("G", c(1:11)))
cluster2celltype <- c(
  "G1"="G1-B cells", "G2"="G2-Plasma cells", "G3"="G3-Monocytes/Macrophages", "G4"="G4-Dendritic cells",
  "G5"="G5-Lymphocytes", "G6"="G6-Exhausted CD8+ T cells", "G7"="G7-Regulatory T cells", 
  "G8"="G8-Cytotoxicity Lymphocytes", "G9"="G9-Exhausted/HS CD8+ T cells", "G10"="G10-Memory T cells",
  "G11"="G11-Lymphocytes exhausted/cell cycle"
)
data_S@meta.data$cell_type <- cluster2celltype[as.character(data_S@meta.data$cluster)]
data_S@meta.data$cell_type <- factor(data_S@meta.data$cell_type, levels = c(
  cluster2celltype
))


# scale data
data_S <- ScaleData(data_S)

# calculate gene variance to set variable genes
gene_var <- apply(data_exp, 1, var)
VariableFeatures(data_S) <- names(gene_var)[gene_var > 6]


# dimension reduction
data_S <- RunPCA(data_S, npcs = 10, verbose = F)

data_S <- runFItSNE(
  data_S, dims.use = 1:10, seed.use = 3, fast.R.path = fitsneR, 
  ann_not_vptree = FALSE, nthreads = 12
)
TSNEPlot(data_S, group.by = "condition", pt.size = 0.5)
TSNEPlot(data_S, group.by = "cluster", label = T, pt.size = 0.5) + 
  scale_color_hue(labels = cluster2celltype)

pca_embedding <- data_S@reductions$pca@cell.embeddings
tsne_embedding <- data_S@reductions$tsne@cell.embeddings





##=============================================##
## DA-Seq

# sample labels for responders and non-responders
labels_res <- X.label.info[X.label.info$condition == "R", "label"]
labels_nonres <- X.label.info[X.label.info$condition == "NR", "label"]



## DA cells

da_cells <- getDAcells(
  X = pca_embedding, 
  cell.labels = data_S@meta.data$lesion, 
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  k.vector = seq(50, 500, 50), 
  plot.embedding = tsne_embedding
)

da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(0.05,0.9), 
  do.plot = T, plot.embedding = tsne_embedding, size = 0.1
)

da_cells$pred.plot
da_cells$da.cells.plot



## DA regions

da_regions <- getDAregion(
  X = pca_embedding, 
  da.cells = da_cells, 
  cell.labels = data_S@meta.data$lesion,
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  resolution = 0.01, 
  plot.embedding = tsne_embedding, size = 0.1
)

da_regions$da.region.plot
da_regions$DA.stat

n_da <- length(unique(da_regions$da.region.label)) - 1



## DA markers

# Seurat
data_S <- addDAslot(data_S, da.regions = da_regions, da.slot = "da")
Seurat_markers <- SeuratMarkerFinder(
  data_S, da.slot = "da", test.use = "negbinom", only.pos = T
)
FeaturePlot(
  data_S, features = c("VCAM1","TYROBP","MS4A1","KLRK1","LEF1"), cols = c("gray","red")
)


# STG
STG_markers <- STGmarkerFinder(
  X = as.matrix(data_S@assays$RNA@data), 
  da.regions = da_regions, 
  lambda = 1.5, n.runs = 5, return.model = T, 
  python.use = python2use, GPU = GPU
)

plotCellScore(
  tsne_embedding, score = STG_markers$model$`4`$pred, cell.col = c("gray","blue")
)



## Local markers

da_clusters <- c("4"="G5")

# Seurat negbinom to find local markers
local_markers <- list()
for(i in 1:n_da){
  if(!as.character(i) %in% names(da_clusters)){next()}
  local_markers[[as.character(i)]] <- SeuratLocalMarkers(
    object = data_S, da.region.to.run = i, cell.label.slot = "cluster", cell.label.to.run = da_clusters[as.character(i)], 
    assay = "RNA", test.use = "negbinom", min.diff.pct = 0.09, only.pos = T
  )
  local_markers[[as.character(i)]]$pct.diff <- local_markers[[as.character(i)]]$pct.1 - 
    local_markers[[as.character(i)]]$pct.2
}





##=============================================##
## Generate plots

library(scales)
da_cols <- hue_pal()(n_da)
da_order <- order(da_regions$da.region.label)

## TSNE plots

gg1 <- plotCellLabel(tsne_embedding, label = data_S@meta.data$condition, size = 0.1, do.label = F) + theme_tsne
ggsave(gg1, filename = "figs/melanoma_a.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg1, legend.position = "top"), 
       filename = "figs/melanoma_a_legend.pdf", width = 2, height = 0.25, dpi = 1200)


gg2 <- plotCellLabel(tsne_embedding, label = data_S@meta.data$cluster, size = 0.1, label.size = 3) + 
  scale_color_hue(labels = cluster2celltype) + theme_tsne
ggsave(gg2, filename = "figs/melanoma_b.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg2), filename = "figs/melanoma_b_legend.pdf", width = 2, height = 1.5, dpi = 1200)


gg3 <- da_cells$pred.plot + theme_tsne
ggsave(gg3, filename = "figs/melanoma_c.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg3, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm")), 
       filename = "figs/melanoma_c_legend.pdf", height = 30, width = 15, units = "mm", dpi = 1200)


gg4 <- plotCellLabel(
  tsne_embedding[da_order,], label = as.character(da_regions$da.region.label[da_order]), 
  size = 0.1, label.size = 3, label.plot = as.character(c(1:n_da))
) + scale_color_manual(
  values = c("gray", da_cols), breaks = c(1:n_da), labels = paste0("DA",c(1:n_da))
) + theme_tsne
ggsave(gg4, filename = "figs/melanoma_d.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg4), filename = "figs/melanoma_d_legend.pdf", width = 0.5, height = 0.7, dpi = 1200, useDingbats=F)



## Dot plot

# markers
marker_genes <- list(
  "1" = c("VCAM1","LAG3","CD27","CD38"),
  "2" = c("CD14","FCER1G","VCAN","LYZ"),
  "3" = c("CD19","MS4A1","IGHM","CD79A"),
  "4" = c("IL7R","KLRK1","CD8A","CCL5"),
  "5" = c("LEF1","TCF7","CCR7","SELL")
)

# add STG information
STG.marker.info <- do.call(rbind, lapply(STG_markers$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(unlist(marker_genes))))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(unlist(marker_genes))
rownames(STG.marker.info) <- c(1:(n_da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]

# generate dot plot
gg5 <- DotPlot(
  data_S, features = unlist(marker_genes), cols = c("gray","blue"), group.by = "da"
) + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value)) + theme_dot + RotatedAxis()
ggsave(gg5, filename = "figs/melanoma_e.pdf", width = 90, height = 40, units = "mm", dpi = 1200)
ggsave(
  g_legend(gg5, legend.key.height = unit(0.15,"cm"), legend.key.width = unit(0.2,"cm")), 
  filename = "figs/melanoma_e_legend.pdf", height = 40, width = 30, units = "mm", dpi = 1200
)

# dot plot for DA4 and G5
i <- 4
data_S@meta.data$da.local <- "0"
data_S@meta.data$da.local[data_S@meta.data$cluster == da_clusters[as.character(i)]] <- 
  da_clusters[as.character(i)]
data_S@meta.data$da.local[data_S@meta.data$da == i] <- paste0("DA",i)
data_S@meta.data$da.local <- factor(data_S@meta.data$da.local, levels = c("0",da_clusters[as.character(i)],paste0("DA",i)))
gg6 <- DotPlot(
  data_S, features = c("CTLA4","FCRL3","KLRC4","CD84"), group.by = "da.local", cols = c("gray","blue")
) + theme_dot + RotatedAxis()
ggsave(gg6, filename = paste0("figs/melanoma_g_DA",i,".pdf"), width = 25, height = 30, units = "mm", dpi = 1200)
data_S@meta.data$da.local <- NULL



## Feature plots

# DA1, VCAM1
sgg1 <- list(
  plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions$da.region.label[da_order]), size = 0.1, do.label = F, 
    cell.col = c("gray",da_cols[1],"gray","gray","gray","gray")
  ) + ggtitle("DA1") + theme_tsne, 
  plotCellScore(
    tsne_embedding, data_S@assays$RNA@data["VCAM1",], cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("VCAM1") + theme_tsne
)
ggsave(
  plot_grid(plotlist = sgg1, nrow = 1), 
  filename = "figs/melanoma_s_a.png", height = 45, width = 80, units = "mm", dpi = 1200
)

# DA4, KRLK1
sgg2 <- list(
  plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions$da.region.label[da_order]), size = 0.01, do.label = F, 
    cell.col = c("gray","gray","gray","gray",da_cols[4],"gray")
  ) + ggtitle("DA4") + theme_tsne,
  plotCellScore(
    tsne_embedding, data_S@assays$RNA@data["KLRK1",], cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("KLRK1") + theme_tsne,
  plotCellScore(
    tsne_embedding, STG_markers$model$`4`$pred, cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("STG_DA4") + theme_tsne
)
ggsave(
  plot_grid(plotlist = sgg2, nrow = 1), 
  filename = "figs/melanoma_s_b.png", height = 45, width = 120, units = "mm", dpi = 1200
)

# DA5, LEF1
sgg3 <- list(
  plotCellLabel(
    tsne_embedding[da_order,], as.factor(da_regions$da.region.label[da_order]), size = 0.1, do.label = F, 
    cell.col = c("gray","gray","gray","gray","gray",da_cols[5])
  ) + ggtitle("DA5") + theme_tsne,
  plotCellScore(
    tsne_embedding, data_S@assays$RNA@data["LEF1",], cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("LEF1") + theme_tsne,
  plotCellScore(
    tsne_embedding, STG_markers$model$`5`$pred, cell.col = c("gray","blue"), size = 0.1
  ) + ggtitle("STG_DA5") + theme_tsne
)
ggsave(
  plot_grid(plotlist = sgg3, nrow = 1), 
  filename = "figs/melanoma_s_c.png", height = 45, width = 120, units = "mm", dpi = 1200
)

# legend
ggsave(g_legend(
  sgg1[[2]], 
  legend.text = element_blank(), legend.title = element_blank(), 
  legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm")
), filename = "figs/melanoma_s_legend.pdf", height = 30, width = 10, units = "mm", dpi = 1200)




