### DA-seq on aging mouse brain data
### Original paper: https://www.nature.com/articles/s41593-019-0491-3
### This script reproduces analysis presented in Figure 4

library(Seurat) # Version 2.3.4
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)
library(viridis)


## Set Python and GPU
python2use <- "/home/henry/henry_env/venv/bin/python"
GPU <- 3


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

table(data_meta$cell_type)
table(data_meta$all_cells_by_age)
rownames(data_meta) <- data_meta[,1]



## Seurat

# create object
data_S <- CreateSeuratObject(
  raw.data = data_exp, names.field = 6, names.delim = "_"
)
table(data_S@meta.data$orig.ident)


# add metadata
data_S@meta.data$cell_type <- data_meta[data_S@cell.names,"cell_type"]
data_S@meta.data$age <- data_meta[data_S@cell.names,"all_cells_by_age"]
data_S@meta.data$label <- data_S@meta.data$age
data_S@meta.data$label[data_S@meta.data$age == "2-3mo"] <- "young"
data_S@meta.data$label[data_S@meta.data$age == "21-22mo"] <- "old"
table(data_S@meta.data$label)

data_S@meta.data$cell_type_num <- as.factor(as.numeric(as.factor(data_S@meta.data$cell_type)))


# analysis
data_S <- ScaleData(data_S)
data_S <- FindVariableGenes(data_S, x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = F)
data_S <- RunPCA(
  data_S, pcs.compute = 20, do.print = F
)
data_S <- RunTSNE(
  data_S, tsne.method = "Rtsne", seed.use = 3, dims.use = 1:20
)

# Fig. 4A
TSNEPlot(data_S, group.by = "age", pt.size = 0.1)
# Fig. 4B
TSNEPlot(data_S, group.by = "cell_type", do.label = T, pt.size = 0.1)





##=============================================##
## DA-Seq


## labels

labels_old <- unique(as.character(data_S@meta.data$orig.ident[data_S@meta.data$age == "21-22mo"]))
labels_young <- unique(as.character(data_S@meta.data$orig.ident[data_S@meta.data$age == "2-3mo"]))



## DA cells

da.cells <- getDAcells(
  X = data_S@dr$pca@cell.embeddings, 
  cell.labels = as.character(data_S@meta.data$orig.ident), 
  labels.1 = labels_young, labels.2 = labels_old, 
  k.vector = seq(50,800,50),
  plot.embedding = data_S@dr$tsne@cell.embeddings,
  python.use = python2use, GPU = GPU
)

da.cells <- updateDAcells(
  da.cells, pred.thres = c(0.05,1),
  do.plot = T, plot.embedding = data_S@dr$tsne@cell.embeddings
)

# Fig. 4C
da.cells$pred.plot
da.cells$da.cells.plot



## DA regions

da.regions <- getDAregion(
  X = data_S@dr$tsne@cell.embeddings, 
  cell.labels = as.character(data_S@meta.data$orig.ident), 
  labels.1 = labels_young, labels.2 = labels_old, 
  cell.idx = da.cells$da.cell.idx, 
  k = 3, alpha = 0.37, restr.fact = 12, 
  do.plot = T, plot.embedding = data_S@dr$tsne@cell.embeddings
)

# Fig. 4D
da.regions$da.region.plot
da.regions$DA.stat



## DA markers

STG_markers <- STGmarkerFinder(
  X = as.matrix(data_S@data), 
  cell.idx = da.cells$da.cell.idx, 
  da.region.label = da.regions$cluster.res, 
  lambda = 1.5, n.runs = 10, return.model = T, 
  python.use = python2use, GPU = GPU 
)



## Dot plot

# add metadata
data_S@meta.data$da <- 0
data_S@meta.data$da[da.cells$da.cell.idx] <- da.regions$cluster.res


# marker genes to plot
goi_list <- list(
  "3" = c("Cdk1","Sox11"),
  "1" = "Pdgfra",
  "2" = "Tmem119"
)
marker_genes <- unlist(goi_list[order(names(goi_list))])


# add STG info
n.da <- length(unique(da.regions$cluster.res)) - 1
STG.marker.info <- do.call(rbind, lapply(STG_markers$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(marker_genes)))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(marker_genes)
rownames(STG.marker.info) <- c(1:(n.da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]
STG.marker.info.m$value <- 7.5 * STG.marker.info.m$value


# generate dot plot
ggdot <- DotPlot(
  data_S, genes.plot = marker_genes, cols = c("gray","red"), group.by = "da", 
  x.lab.rot = T, do.return = T
)

# Fig. 4E
ggdot + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value))



## Run STG on DA region 2 (within cluster 11)

# add meta data
data_S@meta.data$cluster_da <- paste(data_S@meta.data$cell_type_num, data_S@meta.data$da, sep = "_")

# subset Seurat object to keep cells only from cluster 11
data_S_sub <- SubsetData(data_S, subset.name = "cell_type_num", accept.value = 11)


# run STG to separate 11_2
STG_local_marker <- runSTG(
  X = data_S_sub@data, X.labels = data_S_sub@meta.data$cluster_da, 
  label.1 = "11_2", label.2 = "11_0", 
  lambda = 3, n.runs = 1, return.model = T,
  python.use = python2use, GPU = GPU
)


# add STG predictions to Seurat object
data_S_sub@meta.data$STG <- STG_local_marker$model$pred


# Fig. 4F
FeaturePlot(data_S_sub, features.plot = "STG", pt.size = 0.1, do.return = T, no.legend = F)[[1]] + 
  scale_color_viridis_c() + ggtitle(NULL) + 
  xlim(c(min(data_S@dr$tsne@cell.embeddings[,1]),max(data_S@dr$tsne@cell.embeddings[,1]))) + 
  ylim(c(min(data_S@dr$tsne@cell.embeddings[,2]),max(data_S@dr$tsne@cell.embeddings[,2])))


