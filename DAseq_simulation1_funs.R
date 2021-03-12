# function to select DA sites with random size
selectSiteRandSize <- function(
  nn.list, n.site = 4, min.cell = 100, 
  r = NULL, cell.noise = NULL
){
  n.cell <- nrow(nn.list$nn.idx)
  max.cell <- ncol(nn.list$nn.idx)
  if(is.null(r)){
    r <- median(nn.list$nn.dist[,max.cell])
  }
  if(is.null(cell.noise)){
    cell.noise <- ceiling(max.cell / 5)
  }
  
  n.output <- 0
  site.output <- rep(0, n.site)
  site.used <- NULL
  cell.output <- list()
  
  for(i in 1:n.cell){
    if(n.output >= n.site){
      break
    }
    
    # randomly select a site
    site.idx.tmp <- sample(
      x = setdiff(1:n.cell, site.used), 
      size = 1
    )
    site.used <- c(site.used, site.idx.tmp)
    
    # randomly choose site size
    site.cell <- sum(nn.list$nn.dist[site.idx.tmp,] <= r)
    if(site.cell < min.cell | site.cell > max.cell - cell.noise){
      next
    }
    
    site.size <- sample(
      max(min.cell, site.cell - cell.noise):min(max.cell, site.cell + cell.noise), 
      size = 1
    )
    
    # check whether it overlaps with existing sites
    cell.tmp <- nn.list$nn.idx[site.idx.tmp,1:site.size]
    if(length(intersect(cell.tmp, unlist(cell.output))) > 0){
      next
    }
    
    # if not, select this site
    site.output[n.output + 1] <- site.idx.tmp
    cell.output[[n.output + 1]] <- cell.tmp
    n.output <- n.output + 1
  }
  
  return(cell.output)
}


# function to generate labels based on DA sites
generateLabel <- function(da.site, n.label, n.sample, prob = 1, seed = 0){
  set.seed(seed)
  
  n.da <- length(da.site)
  random.label <- sample(c(1:n.sample), size = n.label, replace = T)
  da.label <- random.label
  for(j in 1:n.da){
    if(j <= n.da/2){
      # create labels that belong to group 2
      da.label[da.site[[j]]] <- sample(
        x = 1:n.sample, 
        prob = c(rep((1 - prob)/n.sample*2, n.sample/2), rep(prob/n.sample*2, n.sample/2)), 
        size = length(da.site[[j]]), replace = T
      )
    } else {
      # create labels that belong to group 1
      da.label[da.site[[j]]] <- sample(
        x = 1:n.sample, 
        prob = c(rep(prob/n.sample*2, n.sample/2), rep((1 - prob)/n.sample*2, n.sample/2)), 
        size = length(da.site[[j]]), replace = T
      )
    }
  }
  
  return(da.label)
}


# function to test Cydar
runCydar <- function(data.use, label, n.sample, ...){
  # create data for Cydar
  data.use.list <- list()
  for(ii in 1:n.sample){
    data.use.list[[as.character(ii)]] <- data.use[label == ii,]
  }
  
  label2cond <- c(rep("1",n.sample/2), rep("2",n.sample/2))
  names(label2cond) <- as.character(c(1:n.sample))
  
  
  # Cydar obj
  cydar.obj <- prepareCellData(data.use.list)
  cydar.obj <- countCells(cydar.obj, ...)
  
  # add cell names in cydar object
  for(i in 1:ncol(cydar.obj@cellIntensities)){
    cydar.obj@cellData@rownames[i] <- rownames(
      data.use.list[[cydar.obj@cellData$sample.id[i]]]
    )[cydar.obj@cellData$cell.id[i]]
  }
  
  # DA analysis
  cydar.dge <- DGEList(assay(cydar.obj), lib.size = cydar.obj$totals)
  
  # QL framework (not usable here because of no replicates, use likelihood ratio test instead)
  cydar.dge.design <- model.matrix(~factor(label2cond[colnames(cydar.obj)]))
  cydar.dge <- estimateDisp(cydar.dge, cydar.dge.design)
  cydar.dge.fit <- glmQLFit(y = cydar.dge, design = cydar.dge.design, robust=TRUE)
  cydar.res <- glmQLFTest(cydar.dge.fit, coef = 2)
  cydar.res.table <- cydar.res$table
  
  # spatial-FDR
  cydar.qval <- spatialFDR(intensities(cydar.obj), cydar.res.table$PValue)
  cydar.res.table$qval <- cydar.qval
  
  return(list(
    obj = cydar.obj, 
    res = cydar.res.table
  ))
}
