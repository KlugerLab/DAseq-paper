### Run FIt-SNE with Seurat v3 object
library(rsvd)
runFItSNE <- function(
  object, reduction = "pca", dims.use = 1:10, seed.use = 0, 
  reduction.name = "tsne", reduction.key = "tSNE_", 
  fast.R.path = "~/git/FIt-SNE/fast_tsne.R", ...
){
  current.dir <- getwd()
  source(fast.R.path, chdir = T)
  X.in <- object@reductions[[reduction]]@cell.embeddings[,dims.use]
  X.out <- fftRtsne(
    X = X.in, rand_seed = seed.use, ...
  )
  X.out.ncol <- ncol(X.out)
  colnames(X.out) <- paste(reduction.key, c(1:X.out.ncol), sep = "")
  rownames(X.out) <- colnames(object)
  object@reductions[[reduction.name]] <- new(
    "DimReduc", cell.embeddings = X.out, 
    assay.used = object@reductions[[reduction]]@assay.used, key = reduction.key
  )
  setwd(current.dir)
  return(object)
}



### Plotting themes

library(grid)
library(ggplot2)
library(cowplot)

theme_tsne <- theme_cowplot() + theme(
  text = element_text(size = 7), 
  axis.title = element_blank(), axis.text = element_blank(),
  axis.ticks = element_blank(), plot.title =
    element_text(size = 10, face = "plain"),
  legend.position="none",
  plot.margin=unit(c(0.05,0.05,0.05,0.05 ), "cm"),
  axis.line.y.left=element_line(colour="black",size=0.2),
  axis.line.x.bottom=element_line(colour="black",size=0.2)
)

theme_dot <- theme_cowplot() + theme(
  axis.title = element_blank(), axis.text = element_text(size = 7), 
  legend.position = "none", plot.margin=unit(c(0.05,0.05,0.05,0.05 ), "cm"),
)

g_legend <- function(input.g, ...){
  input.g <- input.g + theme(
    legend.position = "right", legend.text = element_text(size = 7), 
    legend.title = element_text(size = 7), 
    legend.key.width = unit(0.05, 'cm'), legend.key.height = unit(0.05,"cm"),
    legend.spacing = unit(0.25, 'cm')
  ) + theme(...)
  get_legend(input.g)
}