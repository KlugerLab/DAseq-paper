### DA-seq on simulation data 2

library(DAseq)
library(ggplot2)
library(cowplot)

source("convenience.R")

## Set Python and GPU
python2use <- "/data/henry/henry_env/venv/bin/python"
GPU <- 2



## Load data
data.folder <- "./data/sim2/"
data.file <- grep("data", list.files(data.folder, full.names = T), value = T)
label.file <- grep("da",grep("labels", list.files(data.folder, full.names = T), value = T), 
                   invert = T, value = T)
da.label.file <- grep("_da_", list.files(data.folder, full.names = T), value = T)
tsne.file <- grep("tsne", list.files(data.folder, full.names = T), value = T)

data.1 <- as.matrix(read.table(data.file, header = F))
labels.1 <- read.table(label.file, header = F)[,1]
da.labels.1 <- read.table(da.label.file, header = F)[,1]
tsne.1 <- read.table(tsne.file, header = F)

plotCellLabel(
  X = tsne.1, label = as.character(labels.1), do.label = F
)
table(labels.1)

plotCellLabel(
  X = tsne.1, label = as.character(da.labels.1), cell.col = c("gray","red","blue"), do.label = F, size = 1
)
table(da.labels.1)



## DA-seq

da_cells <- getDAcells(
  X = data.1, cell.labels = as.character(labels.1), labels.1 = "0", labels.2 = "1",
  k.vector = round(seq(50,500,50)), plot.embedding = tsne.1
)
da_cells$pred.plot
da_cells$da.cells.plot


da_regions <- getDAregion(
  X = as.matrix(data.1), da.cells = da_cells, 
  cell.labels = as.character(labels.1), labels.1 = "0", labels.2 = "1", 
  plot.embedding = tsne.1
)
da_regions$da.region.plot
table(da_regions$da.region.label)
n_da <- length(unique(da_regions$da.region.label)) - 1


STG_markers <- STGmarkerFinder(
  X = t(data.1), da.regions = da_regions, 
  lambda = 0.5, n.runs = 5, return.model = T, 
  python.use = "/data/henry/henry_env/venv/bin/python", GPU = 2
)
plotCellScore(
  tsne.1, score = data.1[,4], cell.col = c("gray","blue")
)



## Generate plots

library(scales)
da_cols <- hue_pal()(n_da)
da_order <- order(da_regions$da.region.label)

cell_labels <- gsub("1","A",gsub("0","B",labels.1))
gg1 <- plotCellLabel(
  X = tsne.1, label = cell_labels, do.label = F, size = 0.1
) + theme_tsne
ggsave(gg1, filename = "figs/gaussian_a.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg1, legend.position = "right"), 
       filename = "figs/gaussian_a_legend.pdf", width = 0.3, height = 0.3, dpi = 1200)

gg2 <- plotCellLabel(
  X = tsne.1, label = as.character(da.labels.1), do.label = F, size = 0.1
) + scale_color_manual(values = c("gray",rev(hue_pal()(2))), breaks = c("1","2")) + theme_tsne
ggsave(gg2, filename = "figs/gaussian_b.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg2, legend.position = "right"), 
       filename = "figs/gaussian_b_legend.pdf", width = 0.3, height = 0.3, dpi = 1200)

gg3 <- da_cells$pred.plot + theme_tsne
ggsave(gg3, filename = "figs/gaussian_c.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg3, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm"), legend.title = element_blank()), 
       filename = "figs/gaussian_c_legend.pdf", height = 30, width = 15, units = "mm", dpi = 1200)

gg4 <- plotCellLabel(
  tsne.1[da_order,], label = as.character(da_regions$da.region.label[da_order]), 
  size = 0.1, label.size = 2, label.plot = as.character(c(1:n_da))
) + scale_color_manual(
  values = c("gray", da_cols), breaks = c(1:n_da), labels = paste0("DA",c(1:n_da))
) + theme_tsne
ggsave(gg4, filename = "figs/gaussian_d.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg4), filename = "figs/gaussian_d_legend.pdf", width = 0.5, height = 0.7, dpi = 1200)

gg5 <- lapply(c("V5","V10","V4","V9"), function(x){
  plotCellScore(tsne.1, score = data.1[,x], cell.col = c("gray","blue"), size = 0.05) + theme_tsne
})
ggsave(
  plot_grid(plotlist = gg5, nrow = 1), filename = "figs/gaussian_e.png", 
  width = 100, height = 25, units = "mm", dpi = 1200
)
