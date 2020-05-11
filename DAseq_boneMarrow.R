### DA-seq on bone marrow data
### Original paper: https://www.nature.com/articles/s41586-019-1104-8
### This script reproduces analysis presented in Figure 5

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
## Functions

# Combine data from multiple files
combineData <- function(data1, data2, name1 = NULL, name2 = NULL, delim = "_"){
  # read in data
  # exp_pre <- read.table(data1, header = T,sep = "\t",stringsAsFactors = F, row.names = 1)
  # exp_post <- read.table(data2, header = T,sep = "\t",stringsAsFactors = F, row.names = 1)
  exp_pre <- data1
  exp_post <- data2
  if(!is.null(name1)){
    colnames(exp_pre) <- paste(colnames(exp_pre), name1, sep = delim)
  }
  if(!is.null(name2)){
    colnames(exp_post) <- paste(colnames(exp_post), name2, sep = delim)
  }
  
  # find all genes
  allgene <- unique(c(rownames(exp_pre), rownames(exp_post)))
  missG_pre <- setdiff(allgene, rownames(exp_pre))
  missG_post <- setdiff(allgene, rownames(exp_post))
  
  # create matrix with 0 for missing genes
  miss_pre <- matrix(0, ncol = dim(exp_pre)[2], nrow = length(missG_pre))
  rownames(miss_pre) <- missG_pre
  colnames(miss_pre) <- colnames(exp_pre)
  miss_post <- matrix(0, ncol = dim(exp_post)[2], nrow = length(missG_post))
  rownames(miss_post) <- missG_post
  colnames(miss_post) <- colnames(exp_post)
  
  # bind data
  new_pre <- rbind(exp_pre, miss_pre)
  new_post <- rbind(exp_post, miss_post)
  new_pre <- new_pre[order(rownames(new_pre)),]
  new_post <- new_post[order(rownames(new_post)),]
  print(mean(rownames(new_pre) == rownames(new_post)))
  data_new <- cbind(new_pre, new_post)
  return(data_new)
}





##=============================================##
## Data prep

## Load data

if(!dir.exists("./data/")){
  dir.create("./data/")
}

# Gene expression files
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE108892&format=file", 
  "./data/GSE108892_RAW.tar"
)
system(
  "tar -xvf ./data/GSE108892_RAW.tar -C ./data/"
)


# Read data from 5FU and niche
sample_names <- c(
  "GSM2915575_niche-col23", "GSM2915576_niche-lepr", "GSM2915577_niche-vecad",
  "GSM2915578_5FU-cntrl", "GSM2915579_5FU-cntrl-col23", 
  "GSM2915580_5FU-treat-lepr", "GSM2915581_5FU-treat-vecad", "GSM2915582_5FU-treat-col23"
)
n_sample <- length(sample_names)

data_list <- list()
for(i in 1:n_sample){
  data_list[[i]] <- Matrix(as.matrix(read.table(
    paste("./data/", sample_names[i], ".txt.gz", sep = ""),
    sep = "\t", header = T, row.names = 1, stringsAsFactors = F
  )), sparse = T)
}
names(data_list) <- sample_names


# meta data file (provided by first authors of the original paper)
download.file(
  "https://raw.githubusercontent.com/KlugerLab/DAseq-paper/master/data/boneMarrow_meta.csv",
  "./data/boneMarrow_meta.csv"
)

cell_meta <- read.table("./data/boneMarrow_meta.csv", sep = ",", header = T, row.names = 1)

# get sample names
sample_ids <- c("col23","lepr","vecad","CTRL-mix","CTRL-col23","5FU-mix","5FU-vecad","5FU-col23")


# rename cell names
for(i in 1:n_sample){
  colnames(data_list[[i]]) <- sapply(colnames(data_list[[i]]), function(x){
    x.s <- unlist(strsplit(x, split = ".", fixed = T))
    x.l <- length(x.s)
    return(x.s[x.l])
  })
  colnames(data_list[[i]]) <- paste(sample_ids[i], colnames(data_list[[i]]), sep = ":")
}



## Filter & Merge data

# take only cells in the metadata
data_filter_list <- list()
for(i in 1:n_sample){
  data_filter_list[[i]] <- data_list[[i]][,colnames(data_list[[i]]) %in% rownames(cell_meta)]
}
names(data_filter_list) <- sample_ids


# niche data
data_niche <- combineData(
  data1 = combineData(data1 = data_filter_list[[1]], data2 = data_filter_list[[2]]),
  data2 = data_filter_list[[3]]
)

# 5-FU data
data_fu <- combineData(data1 = data_filter_list[[4]], data2 = data_filter_list[[5]])
for(i in 6:8){
  data_fu <- combineData(data1 = data_fu, data2 = data_filter_list[[i]])
}



## Seurat

# niche
niche_S <- CreateSeuratObject(data_niche, project = "niche", names.field = 1, names.delim = ":")
table(niche_S@meta.data$orig.ident)
niche_S@meta.data$batch <- "non-pbs"

niche_S <- NormalizeData(niche_S)
niche_S <- ScaleData(niche_S)
niche_S <- FindVariableGenes(niche_S, do.plot = F)


# 5-FU
fu_S <- CreateSeuratObject(data_fu, project = "fu", names.field = 1, names.delim = ":")
table(fu_S@meta.data$orig.ident)
fu_S@meta.data$batch <- "pbs"

fu_S <- NormalizeData(fu_S)
fu_S <- ScaleData(fu_S)
fu_S <- FindVariableGenes(fu_S, do.plot = F)


# CCA
genes.cca <- union(head(rownames(niche_S@hvg.info), 2000), head(rownames(fu_S@hvg.info), 2000))
genes.cca <- intersect(genes.cca, rownames(niche_S@data))
genes.cca <- intersect(genes.cca, rownames(fu_S@data))

data_S <- RunCCA(
  object = niche_S, object2 = fu_S, genes.use = genes.cca, num.cc = 30
)
data_S <- AlignSubspace(
  data_S, reduction.type = "cca", grouping.var = "batch", dims.align = 1:30
)


# add condition
data_S@meta.data$condition <- cell_meta[rownames(data_S@meta.data),"ident2"]


# t-SNE
data_S <- RunTSNE(
  data_S, reduction.use = "cca.aligned", dims.use = 1:30, seed.use = 3
)
# Fig. 5A
TSNEPlot(data_S, group.by = "condition", pt.size = 0.1)


# add cluster information
data_S@meta.data$cluster <- cell_meta[rownames(data_S@meta.data),"ident11named"]
data_S@meta.data$cluster <- gsub("[[:digit:]]{2}-","",data_S@meta.data$cluster)
table(data_S@meta.data$cluster)
# Fig. 5B
TSNEPlot(data_S, group.by = "cluster", pt.size = 0.1, do.label = T)





##=============================================##
## DA-Seq

## DA cells

da.cells <- getDAcells(
  X = data_S@dr$cca.aligned@cell.embeddings, 
  cell.labels = as.character(data_S@meta.data$condition),
  labels.1 = "CTRL", labels.2 = "FU", 
  k.vector = seq(50,600,50), 
  plot.embedding = data_S@dr$tsne@cell.embeddings,
  python.use = python2use, GPU = GPU
)

da.cells <- updateDAcells(
  da.cells, pred.thres = c(0,0.95),
  do.plot = T, plot.embedding = data_S@dr$tsne@cell.embeddings
)
# Fig. 5C
da.cells$pred.plot
da.cells$da.cells.plot



## DA regions

da.regions <- getDAregion(
  X = data_S@dr$cca.aligned@cell.embeddings[,1:20], 
  cell.labels = as.character(data_S@meta.data$condition),
  labels.1 = "CTRL", labels.2 = "FU", 
  cell.idx = da.cells$da.cell.idx, 
  k = 3, alpha = 0.7, restr.fact = 12, 
  plot.embedding = data_S@dr$tsne@cell.embeddings
)

# Fig. 5D
da.regions$da.region.plot
da.regions$DA.stat

n.da <- length(unique(da.regions$cluster.res)) - 1



## DA Markers

STG_markers <- STGmarkerFinder(
  X = as.matrix(data_S@data), 
  cell.idx = da.cells$da.cell.idx, 
  da.region.label = da.regions$cluster.res, 
  lambda = 1.2, n.runs = 5, 
  python.use = python2use, GPU = GPU
)

goi_list <- list(
  "1" = c("Birc5","Ccnb2","Cenpa","Cenpf"),
  "2" = c("Lpl"), 
  "3" = c("Gas6","Hp")
)
marker_genes <- unlist(goi_list[order(names(goi_list))])



## Dot plot

# add metadata
data_S@meta.data$da <- 0
data_S@meta.data$da[da.cells$da.cell.idx] <- da.regions$cluster.res


# add STG information
STG.marker.info <- do.call(rbind, lapply(STG_markers$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(marker_genes)))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(marker_genes)
rownames(STG.marker.info) <- c(1:(n.da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]
STG.marker.info.m$value <- 10 * STG.marker.info.m$value


# dot plot
ggdot <- DotPlot(
  data_S, genes.plot = marker_genes, cols = c("gray","red"), group.by = "da", 
  x.lab.rot = T, do.return = T
)

# Fig. 5E
ggdot + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value))



