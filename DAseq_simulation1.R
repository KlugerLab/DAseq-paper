### DA-seq on simulation data 1

library(DAseq)
library(ggplot2)
library(cowplot)
library(cydar)

source("DAseq_simulation1_funs.R")
source("convenience.R")

## NOT run ##
# # select sites
# n.da <- 4
# set.seed(3)
# da.site <- selectSiteRandSize(
#   nn.list = cell.nn.list, n.site = n.da
# )
# 
# # create labels
# n.sample <- 48
# da.label <- generateLabel(
#   da.site = da.site, n.label = n.cells, n.sample = n.sample, prob = 0.9
# )



## Load data

n.da <- 4
da.site <- list()
for(i in 1:n.da){
  da.site[[i]] <- read.table(
    file = paste0("data/sim1/site",i), stringsAsFactors = F, header = F
  )[,1]
}

plotDAsite(
  X = X.2d.melanoma, 
  site.list = da.site, 
  size = .5, cols = c("red","pink","lightblue","blue")
)


da.label <- read.table("data/sim1/da.label", stringsAsFactors = F, header = F)[,1]
da.binary.label <- da.label
da.binary.label[da.label %in% c(1:24)] <- 0
da.binary.label[da.label %in% c(25:48)] <- 1
da.binary.label <- factor(da.binary.label, levels = c(1,0))



## DA-seq

da_cells <- getDAcells(
  X = pca_embedding, 
  cell.labels = as.character(da.label), 
  labels.1 = as.character(c(1:24)), labels.2 = as.character(c(25:48)), 
  k.vector = seq(50,500,50), 
  plot.embedding = tsne_embedding
)

da_cells <- updateDAcells(
  X = da_cells, 
  plot.embedding = tsne_embedding, size = 0.1
)
da_cells$pred.plot

da_regions <- getDAregion(
  X = pca_embedding, da.cells = da_cells, 
  cell.labels = as.character(da.label), 
  labels.1 = as.character(c(1:24)), labels.2 = as.character(c(25:48)), 
  plot.embedding = tsne_embedding
)
da_regions$da.region.plot

n_da <- length(unique(da_regions$da.region.label)) - 1



## Cydar

# get Cydar result
cydar.res <- runCydar(data.use = pca_embedding, label = as.character(da.label), n.sample = n.sample, tol = 2)
cydar.centers <- cydar.res[["obj"]]@cellData@rownames[cydar.res[["obj"]]@elementMetadata$center.cell]
cydar.sig.idx <- which(cydar.res[["res"]]$qval < 0.05)
cydar.numcells <- unlist(lapply(unpackIndices(cydar.res$obj@cellAssignments), length))


# get significant cells
cydar.sig.cells <- cydar.res[["obj"]]@cellData@rownames[unique(unlist(unpackIndices(
  cydar.res[["obj"]]@cellAssignments[cydar.sig.idx]
)))]


# Cydar with different tol
cydar_tol <- c(1.5,2,2.5,3,3.5)
cydar_res_list <- lapply(cydar_tol, FUN = function(x){
  runCydar(data.use = pca_embedding, label = as.character(da.label), n.sample = n.sample, 
           tol = x, BPPARAM=SnowParam(workers = 1))
})
names(cydar_res_list) <- as.character(cydar_tol)



## ROC Curve

assessCydar <- function(input.cydar, da.sites, cell.names, fdr.range = seq(0,1,0.01)){
  n.cells <- length(input.cydar$obj@cellData$cell.id)
  pos <- unlist(da.sites)
  neg <- setdiff(c(1:n.cells),pos)
  
  n.fdr <- length(fdr.range)
  tpr <- rep(0,n.fdr)
  fpr <- rep(0,n.fdr)
  for(ii in 1:n.fdr){
    fdr <- fdr.range[ii]
    sig.idx <- which(input.cydar[["res"]]$qval <= fdr)
    input.cells <- match(input.cydar[["obj"]]@cellData@rownames[unique(unlist(unpackIndices(
      input.cydar[["obj"]]@cellAssignments[sig.idx]
    )))], cell.names)
    tp <- length(intersect(input.cells, pos))
    fp <- length(input.cells) - tp
    fn <- length(setdiff(pos, input.cells))
    tn <- length(neg) - fp
    tpr[ii] <- tp/(tp+fn)
    fpr[ii] <- fp/(fp+tn)
  }
  return(list(tpr = tpr, fpr = fpr))
}

assessDA <- function(input.da, da.sites, th.range = seq(0,1,0.01)){
  n.cells <- length(input.da$da.pred)
  pos <- unlist(da.sites)
  neg <- setdiff(c(1:n.cells),pos)
  
  n.th <- length(th.range)
  tpr <- rep(0,n.th)
  fpr <- rep(0,n.th)
  for(ii in 1:n.th){
    th <- th.range[ii]
    input.cells <- updateDAcells(
      X = input.da, pred.thres = c(-th, th), force.thres = T, do.plot = F
    )
    input.cells <- c(input.cells$da.up, input.cells$da.down)
    tp <- length(intersect(input.cells, pos))
    fp <- length(input.cells) - tp
    fn <- length(setdiff(pos, input.cells))
    tn <- length(neg) - fp
    tpr[ii] <- tp/(tp+fn)
    fpr[ii] <- fp/(fp+tn)
  }
  return(list(tpr = tpr, fpr = fpr))
}

da_metric <- assessDA(input.da = da_cells, da.sites = da.site)

cydar_metric <- lapply(cydar_tol, function(x){
  x.m <- assessCydar(input.cydar = cydar_res_list[[as.character(x)]], 
                     da.sites = da.site, cell.names = rownames(pca_embedding))
  data.frame(TPR = x.m$tpr, FPR = x.m$fpr, method = paste0("Cydar tol=",x))
})

plot_roc <- rbind(
  data.frame(TPR = da_metric$tpr, FPR = da_metric$fpr, method = "DA-seq"),
  do.call("rbind", cydar_metric)
)

ggplot(data = plot_roc) + theme_cowplot() + 
  geom_line(aes(FPR, TPR, col = method))





##=============================================##
## Generate plots

library(scales)
da_cols <- hue_pal()(n_da)
da_order <- order(da_regions$da.region.label)

gg1 <- plotDAsite(
  X = tsne_embedding, site.list = da.site, size = .1, cols = c("red","pink","lightblue","blue")
) + theme_tsne
ggsave(gg1, filename = "figs/simulation_s_a.png", width = 50, height = 50, units = "mm", dpi = 1200)

gg2 <- plotCellLabel(tsne_embedding, label = da.binary.label, do.label = F, size = .1) + theme_tsne
ggsave(gg2, filename = "figs/simulation_s_b.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg2), filename = "figs/simulation_s_b_legend.pdf", width = 0.25, height = 0.3, dpi = 1200)

gg3 <- da_cells$pred.plot + theme_tsne
ggsave(gg3, filename = "figs/simulation_s_c.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg3, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm"), legend.title = element_blank()), 
       filename = "figs/simulation_s_c_legend.pdf", height = 30, width = 15, units = "mm", dpi = 1200)

gg4 <- plotCellLabel(
  tsne_embedding[da_order,], label = as.character(da_regions$da.region.label[da_order]), 
  size = 0.1, label.size = 2, label.plot = as.character(c(1:n_da))
) + scale_color_manual(
  values = c("gray", da_cols), breaks = c(1:n_da), labels = paste0("DA",c(1:n_da))
) + theme_tsne
ggsave(gg4, filename = "figs/simulation_s_d.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg4), filename = "figs/simulation_s_d_legend.pdf", width = 0.5, height = 0.7, dpi = 1200)



## Cydar plots

# plot Cydar hypersphere result
gg5 <- ggplot() + geom_point(data = data.frame(
  Dim1 = tsne_embedding[cydar.centers,1], Dim2 = tsne_embedding[cydar.centers,2],
  logFC = cydar.res[["res"]]$logFC, numCells = cydar.numcells
)[cydar.sig.idx, ], aes(Dim1, Dim2, col = logFC, size = numCells), alpha = 0.5, stroke = 0) + 
  scale_color_gradientn(colours = c("blue","white","red")) + 
  guides(size = guide_legend(title = "# Cells", override.aes = list(alpha=1))) + theme_tsne
ggsave(gg5, filename = "figs/simulation_s_e.png", width = 50, height = 50, units = "mm", dpi = 1200)
ggsave(g_legend(gg5, legend.key.height = unit(0.4,"cm"), legend.key.width = unit(0.4,"cm")), 
       filename = "figs/simulation_s_e_legend.pdf", height = 50, width = 15, units = "mm", dpi = 1200)

# plot Cydar cells in sig spheres
gg6 <- plotDAsite(tsne_embedding, site.list = list(cydar.sig.cells), size = 0.1) + theme_tsne
ggsave(gg6, filename = "figs/simulation_s_f.png", width = 50, height = 50, units = "mm", dpi = 1200)




