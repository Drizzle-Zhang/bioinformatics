library(splatter)
library(scater)
library(Seurat)
library(kBET)
library(cluster)
library(FNN)
library(mclust)

# entropy
entropy <- function(vec){
    entropy <- 0
    len.vec <- length(vec)
    for (batch in unique(vec)) {
        len.batch <- length(vec[vec == batch])
        p.batch <- len.batch/len.vec
        entropy <- entropy + p.batch*log(p.batch)
    }
    return(entropy)
}

# mutual information
mutinfo <- function(x, y){
    mutinfo <- 0
    len.x <- length(x)
    len.y <- length(y)
    xy <- data.frame(x = x, y = y)
    for (element.x in unique(x)) {
        p.x <- length(x[x == element.x])/len.x
        for (element.y in unique(y)) {
            p.y <- length(y[y == element.y])/len.y
            p.xy <- dim(xy[(xy$x == element.x) & (xy$y == element.y),])[1]/len.y
            if (p.xy != 0) {
                mutinfo <- mutinfo + p.xy*log(p.xy/(p.x*p.y))
            }
        }
    }
    return(mutinfo)
}

# main function
generate_data_evaluate <- function(
    facLoc = 0.04, facScale = 0, num.batch = 3, num.cell = 500,
    num.gene = 15000, use.pc = 2, use.sig.pc = F){
    
    # simulate
    params <- newSplatParams(
        nGenes = num.gene, 
        batchCells = rep(num.cell, num.batch), 
        batch.facLoc = rep(facLoc, num.batch),
        batch.facScale = rep(facScale, num.batch))
    if (length(num.cell) != 1) {
        if (length(num.cell) != num.batch) {
            print('Error: batch number mismatch')
        } else {
            setParam(params, "batchCells", num.cell)
        }
    }
    sim.data <- splatSimulate(params, verbose = F)
    mtx.count <- counts(sim.data)
    # mtx.count[is.na(mtx.count)] <- 0
    # gene.name <- rownames(sim.data)
    # cell.name <- colnames(sim.data)
    batches <- sim.data$Batch
    
    # creat Seurat object and normlize data
    object <- CreateSeuratObject(mtx.count)
    object@meta.data$batch <- batches
    object <- NormalizeData(object, display.progress = F)
    object <- FindVariableGenes(
        object = object, mean.function = ExpMean, 
        dispersion.function = LogVMR, do.plot = F, display.progress = F,
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    length(object@var.genes)
    object <- ScaleData(
        object, genes.use = object@var.genes, vars.to.regress = c("nUMI"),
        display.progress = F)
    
    # PCA
    object <- RunPCA(
        object = object, pcs.compute = 100, 
        pc.genes = object@var.genes, do.print = F, 
        pcs.print = 1:5, genes.print = 5)
    pca.data <- object@dr$pca@cell.embeddings
    # title.pca <- paste0('Variance parameter: ', facLoc)
    # ggplot.pca <- PCAPlot(
    #     object, group.by = "batch", png.arguments = c(15, 15, 225), 
    #     do.return = T, plot.title = title.pca)
    # ggsave(
    #     plot = ggplot.pca, path = out.path,
    #     filename = paste0('PCA_plot_var_', facLoc, '.png'))

    # TSNE
    # object <- RunTSNE(object = object, perplexity = 10, do.fast = TRUE)
    # title.tsne.lab <- paste0('t-SNE plot')
    # ggo.tsne.lab <-
    #     TSNEPlot(object = object, group.by = 'batch', pt.size = 1,
    #              png.arguments = c(15, 15, 225), do.return = T,
    #              plot.title = title.tsne.lab)
    # ggsave(plot = ggo.tsne.lab, filename = 'tSNE_plot_lab.png',
    #        path = output)
    
    # pcReg
    batch.factor <- as.factor(object@meta.data$batch)
    tol <- 0.05
    p.value <- c()
    r.squared <- c()
    for (i in 1:dim(pca.data)[2]) {
        fit <- summary(lm(pca.data[, i] ~ batch.factor))
        r.squared <- c(r.squared, fit$r.squared)
        p.value <- c(p.value, 
                     1 - pf(fit$fstatistic[1], 
                            fit$fstatistic[2], 
                            fit$fstatistic[3]))
    }
    p.adj.value <- p.adjust(p.value)
    # select number of PC
    if (use.sig.pc) {
        pc.idx <- which(p.adj.value < tol)
        use.pc <- max(pc.idx, use.pc)
    }
    use.pca.data <- pca.data[,1:use.pc]
    # calculate pcReg
    sdev <- object@dr$pca@sdev[1:use.pc]
    var <- sdev^2 / sum(sdev^2)
    pcReg.var <- sum((var * r.squared[1:use.pc])[p.adj.value[1:use.pc] < tol])
    
    # k-means 
    fit.km <- kmeans(use.pca.data, centers = length(unique(batch.factor)), 
                     iter.max = 20)
    batch.km <- as.numeric(fit.km$cluster)
    batch.init <- as.numeric(batch.factor)
    
    # sil
    dissimilar.dist <- dist(use.pca.data)
    sil.out <- silhouette(batch.init, dissimilar.dist)
    sil <- abs(summary(sil.out)$avg.width)
    
    # kBET
    k0 <- floor(mean(table(batch.init))) #neighbourhood size: mean batch size 
    knn.kBET <- get.knn(use.pca.data, k = k0, algorithm = 'cover_tree')
    batch.estimate <- 
        kBET(use.pca.data, batch.init, k0 = k0, knn = knn.kBET, plot = F)
    RR.kBET <- batch.estimate$summary$kBET.observed[1] -
        batch.estimate$summary$kBET.expected[1]
    # AR.kBET <- 1 - (batch.estimate$summary$kBET.observed[1] -
    #                     batch.estimate$summary$kBET.expected[1])
    
    # Entropy of batch mixing
    select.cells <- sample(dimnames(mtx.count)[2][[1]], 100)
    query.data <- use.pca.data[select.cells,]
    knn.mixent <- get.knnx(
        use.pca.data, query = query.data, k = 100, algorithm = 'cover_tree')
    knn.idx <- knn.mixent$nn.index
    list.knn <- list()
    for (i in 1:dim(knn.idx)[1]) {
        list.knn[[i]] <- batch.init[knn.idx[i,]]
    }
    vec.entropy <- lapply(list.knn, entropy)
    mix.entropy <- mean(unlist(vec.entropy))
    # entropy of data mixed completely
    p.batch <- params@batchCells / params@nCells
    max.entropy <- sum(p.batch*log(p.batch))
    # normalize mixing entropy
    norm.mixent <- 1 - mix.entropy/max.entropy
    
    # ARI
    ARI <- adjustedRandIndex(batch.init, batch.km)
    
    # NMI
    NMI <- mutinfo(batch.init, batch.km) / 
        sqrt(entropy(batch.init) * entropy(batch.km))
    
    # output
    access.res <- 
        data.frame(pcReg.var = pcReg.var, sil = sil, RR.kBET = RR.kBET,
                   norm.mixent = norm.mixent, ARI = ARI, NMI = NMI)
    output <- list(access.res = access.res, object.pca = object)
    
    return(output)
    
}
