### DA-seq on melanoma immune cells data
### Original paper: https://www.sciencedirect.com/science/article/pii/S0092867418313941
### This script reproduces analysis presented in Figure 2

library(Seurat) # Version 2.3.4
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)


## Set Python and GPU
python2use <- "/home/henry/henry_env/venv/bin/python"
GPU <- 3


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
  raw.data = data_exp, project = "melanoma.immune"
)


# set metadata for each cell
data_S@meta.data$condition <- patient_info[data_S@cell.names, "characteristics..response"]
data_S@meta.data$lesion <- patient_info[data_S@cell.names, 
                                        "characteristics..patinet.ID..Pre.baseline..Post..on.treatment."]
data_S@meta.data$cluster <- cluster_info$Cluster.number

data_S <- ScaleData(data_S)


# calculate gene variance to set variable genes
gene_var <- apply(data_exp, 1, var)
data_S@var.genes <- names(gene_var)[gene_var > 6]


# dimension reduction
data_S <- RunPCA(data_S, pcs.compute = 10, do.print = F)

data_S <- RunTSNE(data_S, dims.use = 1:10)
TSNEPlot(data_S, group.by = "condition", pt.size = 0.5)
TSNEPlot(data_S, group.by = "cluster", do.label = T, pt.size = 0.5)

pca_embedding <- data_S@dr$pca@cell.embeddings
tsne_embedding <- data_S@dr$tsne@cell.embeddings





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
  plot.embedding = tsne_embedding, 
  python.use = python2use, GPU = GPU
)

da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(0.075,0.925), 
  do.plot = T, plot.embedding = tsne_embedding
)

# Fig. 2C
da_cells$pred.plot
da_cells$da.cells.plot



## DA regions

da_regions <- getDAregion(
  X = pca_embedding, 
  cell.idx = da_cells$da.cell.idx, 
  k = 5, alpha = 0.25, iter.max = 30, 
  cell.labels = data_S@meta.data$lesion,
  labels.1 = labels_res, 
  labels.2 = labels_nonres, 
  plot.embedding = tsne_embedding
)

# Fig. 2D
da_regions$da.region.plot
da_regions$DA.stat

n.da <- length(unique(da_regions$cluster.res)) - 1



## DA markers

STG.genes <- STGmarkerFinder(
  X = as.matrix(data_S@data), 
  cell.idx = da_cells$da.cell.idx, 
  da.region.label = da_regions$cluster.res, 
  lambda = 1.2, n.runs = 5, return.model = T, 
  python.use = python2use, GPU = GPU
)



## Dot plot

# markers
marker_genes <- c(
  "CD19","MS4A1","IGHM","CD79A",
  "VCAM1","LAG3","CD27","CD38",
  "CD14","CSF3R","VCAN","LYZ",
  "IL7R","TCF7","CD8A","CCL5",
  "CCR7","LEF1","SELL"
)


# add metadata
data_S@meta.data$da <- 0
data_S@meta.data$da[da_cells$da.cell.idx] <- da_regions$cluster.res


# add STG information
STG.marker.info <- do.call(rbind, lapply(STG.genes$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(marker_genes)))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(marker_genes)
rownames(STG.marker.info) <- c(1:(n.da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]
STG.marker.info.m$value <- 10 * STG.marker.info.m$value


# generate dot plot
ggdot <- DotPlot(
  data_S, genes.plot = marker_genes, cols = c("gray","red"), group.by = "da", 
  x.lab.rot = T, do.return = T
)

# Fig. 2E
ggdot + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value))




