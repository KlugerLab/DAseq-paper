### DA-seq on mouse skin data
### Original paper: https://www.sciencedirect.com/science/article/pii/S1534580718309882
### This script reproduces analysis presented in Figure 3

library(Seurat) #V3
library(DAseq)
library(Matrix)
library(reshape2)
library(ggplot2)
library(cowplot)

source('convenience.R')

## Set Python and GPU
python2use <- "/data/henry/henry_env/venv/bin/python"
GPU <- 2

## Set path for FIt-SNE R wrapper
fitsneR <- "~/git/FIt-SNE/fast_tsne.R"


##=============================================##
## Functions

# Read in 10X data in .mtx format and add col and row names
read10X <- function(folder = NULL, matFile = NULL, geneFile = NULL, cellFile = NULL, 
                    suffix = "", sep = "_", gz = T){
  if(!is.null(folder)){
    if(gz){
      matFile <- paste(folder, "/matrix.mtx.gz", sep = "")
      geneFile <- paste(folder, "/genes.tsv.gz", sep = "")
      cellFile <- paste(folder, "/barcodes.tsv.gz", sep = "")
    } else {
      matFile <- paste(folder, "/matrix.mtx", sep = "")
      geneFile <- paste(folder, "/genes.tsv", sep = "")
      cellFile <- paste(folder, "/barcodes.tsv", sep = "")
    }
  }
  
  geneNames <- read.table(geneFile, header = F, sep = "\t", as.is = T)[,2]
  cellNames <- paste(read.table(cellFile, header = F, sep = "\t", as.is = T)[,1], suffix, sep = sep)
  
  # add suffix to duplicate gene names
  if(max(table(geneNames)) > 1){
    for(dupgene in names(which(table(geneNames) != 1))){
      geneidx <- which(geneNames == dupgene)
      for(ii in 2:length(geneidx)){
        geneNames[geneidx[ii]] <- paste(dupgene, ii-1, sep = ".")
      }
    }
  }
  
  rawMat <- readMM(matFile)
  rownames(rawMat) <- geneNames
  colnames(rawMat) <- cellNames
  
  return(rawMat)
}



##=============================================##
## Data prep

## Load data

if(!dir.exists("./data/")){
  dir.create("./data/")
}


# Gene expression files
download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE122043&format=file", 
  destfile = "./data/GSE122043_RAW.tar"
)
system(
  "tar -xvf ./data/GSE122043_RAW.tar -C ./data/"
)

download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE122nnn/GSE122043/suppl/GSE122043_genes.tsv.gz",
  "./data/GSE122043_genes.tsv.gz"
)


# Read data of E13.5 and E14.5 samples
sample_names <- c(
  "GSM3453535_e13Control","GSM3453536_e13Control_replicate",
  "GSM3453537_e14Control","GSM3453538_e14Control_replicate"
)
n_sample <- length(sample_names)

data_list <- list()
for(i in 1:n_sample){
  data_list[[i]] <- read10X(
    matFile = paste("./data/", sample_names[i], "_matrix.mtx.gz", sep = ""), 
    cellFile = paste("./data/", sample_names[i], "_barcodes.tsv.gz", sep = ""),
    geneFile = "./data/GSE122043_genes.tsv.gz", 
    suffix = sample_names[i], sep = "-"
  )
}
names(data_list) <- sample_names



## Seurat

# Seurat object for each sample
data_S_list <- lapply(data_list, function(x){
  x_S <- CreateSeuratObject(counts =  x, min.features = 1000, names.delim = "-", names.field = 3)
  x_S <- subset(x_S, subset = nCount_RNA > 2500 & nCount_RNA < 50000)
  x_S <- NormalizeData(x_S)
  x_S <- ScaleData(x_S)
  x_S <- FindVariableFeatures(x_S, selection.method = "mvp", do.plot = F)
  x_S <- RunPCA(
    x_S, features = rownames(subset(x_S@assays$RNA@meta.features, mvp.dispersion > 0.8)),
    npcs = 10, verbose = F
  )
  x_S <- RunTSNE(x_S, dims = 1:10)
  x_S <- FindNeighbors(x_S, dims = 1:10, verbose = F)
  x_S <- FindClusters(x_S, resolution = 0.1, verbose = F)
  
  return(x_S)
})
names(data_S_list) <- sample_names

lapply(data_S_list, function(x){
  table(x@meta.data$orig.ident)
})
lapply(data_S_list, function(x) DotPlot(x, features = c("Col1a1","Krt14","Krt10"), cols = c("gray","red")))

# select only dermal cells based on Col1a1
dermal_cluster <- lapply(data_S_list, function(x){
  gene.ratio.cluster <- by(x@assays$RNA@data["Col1a1",], INDICES = x@active.ident, 
                           FUN = function(xx) mean(xx > 0))
  gene.ratio.cluster.2 <- by(x@assays$RNA@data["Krt14",] + x@assays$RNA@data["Krt10",], 
                             INDICES = x@active.ident, 
                             FUN = function(xx) mean(xx > 0))
  return(names(gene.ratio.cluster)[gene.ratio.cluster > 0.8 & gene.ratio.cluster.2 < 0.5])
})
data_derm_S_list <- list()
for(i in 1:n_sample){
  data_derm_S_list[[i]] <- subset(
    data_S_list[[i]], cells = which(data_S_list[[i]]@active.ident %in% dermal_cluster[[i]])
  )
}
names(data_derm_S_list) <- sample_names



## Merge data

data_anchors <- FindIntegrationAnchors(object.list = data_derm_S_list)
data_S <- IntegrateData(data_anchors, normalization.method = "LogNormalize")
data_S <- ScaleData(data_S)
data_S <- RunPCA(data_S, verbose = F)
plot(data_S@reductions$pca@stdev)

data_S@meta.data$time <- gsub("Control","",sapply(data_S@meta.data$orig.ident, FUN = function(x){
  unlist(strsplit(x, split = "_", fixed = T))[2]
}))
data_S@meta.data$time <- factor(data_S@meta.data$time, levels = c("e14","e13"))

data_S <- runFItSNE(
  data_S, dims.use = 1:40, seed.use = 3, fast.R.path = fitsneR, 
  ann_not_vptree = FALSE, nthreads = 12
)
TSNEPlot(data_S, group.by = "time")
FeaturePlot(data_S, features = paste("rna",c("Col1a1","Krt14","Krt10"),sep="_"), cols = c("gray","red"))
FeaturePlot(data_S, features = paste("rna",c("Dkk1","Lef1","Ptch1","Bmp4"),sep="_"), cols = c("gray","red"))





##=============================================##
## DA-seq

## DA cells

da_cells <- getDAcells(
  X = data_S@reductions$pca@cell.embeddings[,1:40], k.vector = seq(50,500,50), 
  cell.labels = data_S@meta.data$orig.ident,
  labels.1 = c("GSM3453535_e13Control","GSM3453536_e13Control_replicate"),
  labels.2 = c("GSM3453537_e14Control","GSM3453538_e14Control_replicate"),
  plot.embedding = data_S@reductions$tsne@cell.embeddings
)

da_cells <- updateDAcells(
  X = da_cells, pred.thres = c(-0.8,0.8), 
  plot.embedding = data_S@reductions$tsne@cell.embeddings, size = 0.1
)
da_cells$pred.plot
da_cells$da.cells.plot



## DA regions

da_regions <- getDAregion(
  X = data_S@reductions$pca@cell.embeddings[,1:40], da.cells = da_cells, 
  cell.labels = data_S@meta.data$orig.ident,
  labels.1 = c("GSM3453535_e13Control","GSM3453536_e13Control_replicate"),
  labels.2 = c("GSM3453537_e14Control","GSM3453538_e14Control_replicate"), 
  resolution = 0.05, min.cell = 50, 
  plot.embedding = data_S@reductions$tsne@cell.embeddings
)
da_regions$da.region.plot
n_da <- length(unique(da_regions$da.region.label)) - 1



## DA-score per pairwise comparison

getDAscoreOnly <- function(cell.labels, cell.idx, labels.1, labels.2){
  labels.1 <- labels.1[labels.1 %in% cell.labels]
  labels.2 <- labels.2[labels.2 %in% cell.labels]
  
  idx.label <- cell.labels[cell.idx]
  ratio.1 <- sum(idx.label %in% labels.1) / sum(cell.labels %in% labels.1)
  ratio.2 <- sum(idx.label %in% labels.2) / sum(cell.labels %in% labels.2)
  ratio.diff <- (ratio.2 - ratio.1) / (ratio.2 + ratio.1)
  
  return(ratio.diff)
}

for(i in 1:n_da){
  cat("DA", i, ": ",
    getDAscoreOnly(
      cell.labels = data_S@meta.data$orig.ident,
      cell.idx = which(da_regions$da.region.label == i),
      labels.1 = "GSM3453535_e13Control", labels.2 = "GSM3453537_e14Control"
    ), ", ",
    getDAscoreOnly(
      cell.labels = data_S@meta.data$orig.ident,
      cell.idx = which(da_regions$da.region.label == i),
      labels.1 = "GSM3453535_e13Control", labels.2 = "GSM3453538_e14Control_replicate"
    ), ", ",
    getDAscoreOnly(
      cell.labels = data_S@meta.data$orig.ident,
      cell.idx = which(da_regions$da.region.label == i),
      labels.1 = "GSM3453536_e13Control_replicate", labels.2 = "GSM3453537_e14Control"
    ), ", ",
    getDAscoreOnly(
      cell.labels = data_S@meta.data$orig.ident,
      cell.idx = which(da_regions$da.region.label == i),
      labels.1 = "GSM3453536_e13Control_replicate", labels.2 = "GSM3453538_e14Control_replicate"
    ), "\n", sep = ""
  )
}



## DA markers

# Seurat
data_S <- addDAslot(data_S, da.regions = da_regions, da.slot = "da")
Seurat_markers <- SeuratMarkerFinder(
  data_S, da.slot = "da", assay = "RNA", test.use = "negbinom", only.pos = T
)


# STG
STG_markers <- STGmarkerFinder(
  X = as.matrix(data_S@assays$RNA@data), 
  da.regions = da_regions, 
  lambda = 1.5, n.runs = 5, return.model = T, 
  python.use = python2use, GPU = GPU
)





##=============================================##
## Data from Fan et al.
### https://www.sciencedirect.com/science/article/pii/S1534580718306804

download.file(
  "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE102086&format=file",
  "data/GSE102086_RAW.tar"
)
system(
  "tar -xvf ./data/GSE102086_RAW.tar -C ./data/"
)
download.file(
  "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE102nnn/GSE102086/suppl/GSE102086_genes.tsv.gz",
  "data/GSE102086_genes.tsv.gz"
)


# Read data
sample_names_fan <- c(
  "GSM2723623_E13_WT","GSM2723627_E15_WT"
)
n_sample_fan <- length(sample_names_fan)

data_list_fan <- list()
for(i in 1:n_sample_fan){
  data_list_fan[[i]] <- read10X(
    matFile = paste("./data/", sample_names_fan[i], "_matrix.mtx.gz", sep = ""), 
    cellFile = paste("./data/", sample_names_fan[i], "_barcodes.tsv.gz", sep = ""),
    geneFile = "./data/GSE102086_genes.tsv.gz", 
    suffix = sample_names_fan[i], sep = "-"
  )
}
names(data_list_fan) <- sample_names_fan


# Seurat object for each sample
data_S_list_fan <- lapply(data_list_fan, function(x){
  x_S <- CreateSeuratObject(counts =  x, min.features = 0, names.delim = "-", names.field = 3)
  x_S <- NormalizeData(x_S)
  x_S <- ScaleData(x_S)
  x_S <- FindVariableFeatures(x_S, selection.method = "mvp", do.plot = F)
  x_S <- RunPCA(x_S, npcs = 10, verbose = F)
  x_S <- RunTSNE(x_S, dims = 1:10)
  x_S <- FindNeighbors(x_S, dims = 1:10, verbose = F)
  x_S <- FindClusters(x_S, resolution = 0.1, verbose = F)
  
  return(x_S)
})
names(data_S_list_fan) <- sample_names_fan


# select only dermal cells based on Col1a1
dermal_cluster_fan <- list(
  c("0"), c("0","4")
)
data_derm_S_list_fan <- list()
for(i in 1:n_sample_fan){
  data_derm_S_list_fan[[i]] <- subset(
    data_S_list_fan[[i]], cells = which(data_S_list_fan[[i]]@active.ident %in% dermal_cluster_fan[[i]])
  )
}
names(data_derm_S_list_fan) <- sample_names_fan


# merge data
data_S_fan <- merge(x = data_derm_S_list_fan[[1]], y = data_derm_S_list_fan[[2]])
data_S_fan <- ScaleData(data_S_fan)

data_S_fan@meta.data$time <- gsub("_WT","",gsub("GSM[[:digit:]]+_","",data_S_fan@meta.data$orig.ident))
data_S_fan@meta.data$time <- factor(data_S_fan@meta.data$time, levels = c("E15","E13"))



## Gene modules

da_gene_modules <- lapply(Seurat_markers, FUN = function(x){
  rownames(x)[1:min(100,nrow(x))]
})
names(da_gene_modules) <- paste0("DA", names(da_gene_modules))

data_S_fan <- AddModuleScore(
  data_S_fan, features = da_gene_modules, assay = "RNA", name = names(da_gene_modules)
)
for(i in 1:n_da){
  colnames(data_S_fan@meta.data)[grep(names(da_gene_modules)[i],colnames(data_S_fan@meta.data))] <- 
    names(da_gene_modules)[i]
}





##=============================================##
## Generate plots

library(scales)
da_cols <- hue_pal()(n_da)
da_order <- order(da_regions$da.region.label)

tsne_embedding <- data_S@reductions$tsne@cell.embeddings

## TSNE plots

gg1 <- plotCellLabel(tsne_embedding, label = data_S@meta.data$time, size = 0.1, do.label = F) + theme_tsne
ggsave(gg1, filename = "figs/mouseSkin_a.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg1, legend.position = "right"), 
       filename = "figs/mouseSkin_a_legend.pdf", width = 0.5, height = 0.3, dpi = 1200)


gg2 <- da_cells$pred.plot + theme_tsne
ggsave(gg2, filename = "figs/mouseSkin_b.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg2, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm")), 
       filename = "figs/mouseSkin_b_legend.pdf", height = 30, width = 15, units = "mm", dpi = 1200)


gg3 <- plotCellLabel(
  tsne_embedding[da_order,], label = as.character(da_regions$da.region.label[da_order]), 
  size = 0.1, label.size = 2, label.plot = as.character(c(1:n_da))
) + scale_color_manual(
  values = c("gray", da_cols), breaks = c(1:n_da), labels = paste0("DA",c(1:n_da))
) + theme_tsne
ggsave(gg3, filename = "figs/mouseSkin_c.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg3), filename = "figs/mouseSkin_c_legend.pdf", width = 0.5, height = 0.7, dpi = 1200)



## Feature plot

gg4 <- plotCellScore(
  tsne_embedding, data_S@assays$RNA@data["Sox2",], cell.col = c("gray","blue"), size = 0.1
) + theme_tsne
ggsave(gg4, filename = "figs/mouseSkin_d.png", width = 50, height = 50, units = "mm", dpi = 1200)



## Dot plot

marker_genes <- list(
  "1" = c("Dkk1","Wif1"),
  "2" = c("Sox2","Cdkn1a","Bmp4","Ptch1"),
  "3" = c("Mgp","Six1"),
  "4" = c("Pitx2","Pitx1"),
  "5" = c("Emx2","Osr1")
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
DefaultAssay(data_S) <- "RNA"
gg5 <- DotPlot(
  data_S, features = unlist(marker_genes), cols = c("gray","blue"), group.by = "da"
) + guides(
    color = guide_colorbar(title = "Average Expression", order = 2), 
    size = guide_legend(title = "Percent Expressed", order = 1)
) + theme_dot + RotatedAxis()
ggsave(gg5, filename = "figs/mouseSkin_e.pdf", width = 80, height = 40, units = "mm", dpi = 1200)
ggsave(
  g_legend(gg5, legend.key.height = unit(0.15,"cm"), legend.key.width = unit(0.2,"cm"), legend.spacing = unit(0.5, 'cm')), 
  filename = "figs/mouseSkin_e_legend.pdf", height = 40, width = 30, units = "mm", dpi = 1200
)



## Violin plots
gg6 <- lapply(c(1:n_da), FUN = function(x){
  VlnPlot(data_S_fan, features = paste0("DA",x), group.by = "time", pt.size = 0) + 
    theme_tsne + ggtitle("")
})
ggsave(
  plot_grid(plotlist = gg6, nrow = 1), filename = "figs/mouseSkin_f.pdf", 
  width = 150, height = 35, units = "mm", dpi = 1200
)

# p-value
sapply(c(1:5), function(x){
  wilcox.test(x = data_S_fan@meta.data[which(data_S_fan@meta.data$time == "E15"),paste0("DA",x)], 
              y = data_S_fan@meta.data[which(data_S_fan@meta.data$time == "E13"),paste0("DA",x)])$p.value
})



