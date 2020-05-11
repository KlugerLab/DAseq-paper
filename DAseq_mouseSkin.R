### DA-seq on mouse skin data
### Original paper: https://www.sciencedirect.com/science/article/pii/S1534580718309882
### This script reproduces analysis presented in Figure 3

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
  x_S <- CreateSeuratObject(raw.data = x, min.genes = 1000, names.delim = "-", names.field = 3)
  x_S <- NormalizeData(x_S)
  x_S <- ScaleData(x_S)
  x_S <- FindVariableGenes(x_S, do.plot = F)
  x_S <- RunPCA(
    x_S, pc.genes = rownames(subset(x_S@hvg.info, gene.dispersion > 0.8)),
    pcs.compute = 10, do.print = F
  )
  x_S <- RunTSNE(x_S, dims.use = 1:10)
  x_S <- FindClusters(x_S, dims.use = 1:10, resolution = 0.1, print.output = F)
  
  return(x_S)
})
names(data_S_list) <- sample_names

lapply(data_S_list, function(x){
  table(x@meta.data$orig.ident)
})


# select only dermal cells based on Col1a1
dermal_cluster <- lapply(data_S_list, function(x){
  gene.ratio.cluster <- by(x@data["Col1a1",], INDICES = x@ident, FUN = function(xx) mean(xx > 0))
  return(names(gene.ratio.cluster)[gene.ratio.cluster > 0.85])
})
data_derm_S_list <- list()
for(i in 1:n_sample){
  data_derm_S_list[[i]] <- SubsetData(
    data_S_list[[i]], 
    cells.use = data_S_list[[i]]@cell.names[data_S_list[[i]]@ident %in% dermal_cluster[[i]]],
    subset.raw = T
  )
}



## Merge data

# merge replicates
E13_derm_S <- RunCCA(
  object = data_derm_S_list[[1]], object2 = data_derm_S_list[[2]], num.cc = 15
)
E13_derm_S@meta.data$time <- "E13.5"
E14_derm_S <- RunCCA(
  object = data_derm_S_list[[3]], object2 = data_derm_S_list[[4]], num.cc = 15
)
E14_derm_S@meta.data$time <- "E14.5"

data_derm_S <- RunCCA(
  object = E13_derm_S, object2 = E14_derm_S, num.cc = 15, group.by = "orig.ident"
)


# align CCA space
data_derm_S <- AlignSubspace(data_derm_S, grouping.var = "orig.ident", dims.align = 1:15)


# t-SNE on CCA.aligned
data_derm_S <- RunTSNE(
  data_derm_S, reduction.use = "cca.aligned", dims.use = 1:15, seed.use = 3
)

# Fig. 2A
TSNEPlot(data_derm_S, group.by = "time", pt.size = 0.1)

cca_embedding <- data_derm_S@dr$cca.aligned@cell.embeddings
tsne_embedding <- data_derm_S@dr$tsne@cell.embeddings





##=============================================##
## DA-Seq

## Get DA cells

da.cells <- getDAcells(
  X = cca_embedding, 
  cell.labels = data_derm_S@meta.data$orig.ident,
  labels.1 = c("GSM3453535_e13Control","GSM3453536_e13Control_replicate"), 
  labels.2 = c("GSM3453537_e14Control","GSM3453538_e14Control_replicate"), 
  k.vector = seq(50,500,50), 
  do.plot = T, plot.embedding = tsne_embedding,
  python.use = python2use, GPU = GPU
)

# Fig. 2B
da.cells$pred.plot
da.cells$da.cells.plot



## Get DA regions

da.regions <- getDAregion(
  X = cca_embedding, cell.idx = da.cells$da.cell.idx, 
  k = 2, alpha = 0.7, iter.max = 30, restr.fact = 1,
  cell.labels = data_derm_S@meta.data$orig.ident,
  labels.1 = c("GSM3453535_e13Control","GSM3453536_e13Control_replicate"), 
  labels.2 = c("GSM3453537_e14Control","GSM3453538_e14Control_replicate"), 
  plot.embedding = tsne_embedding
)

# Fig. 2C
da.regions$da.region.plot
da.regions$DA.stat

n.da <- length(unique(da.regions$cluster.res)) - 1



## DA markers

STG.markers <- STGmarkerFinder(
  X = as.matrix(data_derm_S@data), 
  cell.idx = da.cells$da.cell.idx, 
  da.region.label = da.regions$cluster.res, 
  lambda = 1.5, n.runs = 10, 
  python.use = python2use, GPU = GPU
)

# Fig. 2D
FeaturePlot(
  data_derm_S, features.plot = "Sox2", cols.use = c("gray","red"), pt.size = 0.1
)



## Dot plot

marker_genes <- c("Scx","Meox1","Mpped2","Pitx2","Lef1","Bmp4","Sox2","Ptch1")

# add "da" slot
data_derm_S@meta.data$da <- 0
data_derm_S@meta.data$da[da.cells$da.cell.idx] <- da.regions$cluster.res

# add STG information
STG.marker.info <- do.call(rbind, lapply(STG.markers$da.markers, function(x,inputgenes){
  as.numeric(inputgenes %in% x$gene)
}, inputgenes = rev(marker_genes)))
STG.marker.info <- rbind(0, STG.marker.info)
colnames(STG.marker.info) <- rev(marker_genes)
rownames(STG.marker.info) <- c(1:(n.da+1))
STG.marker.info[STG.marker.info == 0] <- NA

STG.marker.info.m <- melt(STG.marker.info)
STG.marker.info.m <- STG.marker.info.m[-which(is.na(STG.marker.info.m$value)),]
STG.marker.info.m$value <- 8 * STG.marker.info.m$value


# dot plot
ggdot <- DotPlot(
  data_derm_S, genes.plot = marker_genes, cols = c("gray","red"), group.by = "da", 
  x.lab.rot = T, do.return = T
)

# Fig. 2E
ggdot + geom_point(data = STG.marker.info.m, aes(x = Var2, y = Var1, size = value))


