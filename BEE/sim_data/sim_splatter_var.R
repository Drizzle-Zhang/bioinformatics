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

# evaluate multiple methods for different variances
num.batch <- 3
num.cell <- 500
num.gene <- 15000
vector.facLoc <- c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.1)
vector.facScale <- 0

out.path <- '/home/drizzle_zhang/Desktop/single_cell/BEE/sim_data/variance'


for (i in 1:length(vector.facScale)) {
    # simulate
    params <- newSplatParams(
        nGenes = num.gene, 
        batchCells = rep(num.cell, num.batch), 
        batch.facLoc = rep(vector.facLoc[i], num.batch),
        batch.facScale = rep(vector.facScale, num.batch))
    sim.data <- splatSimulate(params, verbose = F)
    mtx.count <- counts(sim.data)
    # mtx.count[is.na(mtx.count)] <- 0
    # gene.name <- rownames(sim.data)
    # cell.name <- colnames(sim.data)
    batches <- sim.data$Batch
    
    # creat Seurat object and normlize data
    object <- CreateSeuratObject(mtx.count)
    object@meta.data$batch <- batches
    object <- NormalizeData(object)
    object <- FindVariableGenes(
        object = object, mean.function = ExpMean, 
        dispersion.function = LogVMR, do.plot = F,
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5)
    length(object@var.genes)
    object <- ScaleData(
        object, genes.use = object@var.genes, 
        vars.to.regress = c("nUMI"), num.cores = 2, do.par = TRUE)
    
    # PCA
    object <- RunPCA(
        object = object, pcs.compute = 50, 
        pc.genes = object@var.genes, do.print = F, 
        pcs.print = 1:5, genes.print = 5)
    pca.data <- object@dr$pca@cell.embeddings
    title.pca <- paste0('Variance parameter: ', vector.facLoc[i])
    ggplot.pca <- PCAPlot(
        object, group.by = "batch", png.arguments = c(15, 15, 225), 
        do.return = T, plot.title = title.pca)
    ggsave(
        plot = ggplot.pca, path = out.path,
           filename = paste0('PCA_plot_var_', vector.facLoc[i], '.png'))
    
    # TSNE
    # object <- RunTSNE(object = object, perplexity = 10, do.fast = TRUE)
    # title.tsne.lab <- paste0('t-SNE plot ', title.suffix, ' grouped by lab')
    # ggo.tsne.lab <- 
    #     TSNEPlot(object = object.scale,group.by = 'lab',pt.size = 1.5, 
    #              png.arguments = c(15, 15, 225), do.return = T,
    #              plot.title = title.tsne.lab)
    # ggsave(plot = ggo.tsne.lab, filename = 'tSNE_plot_lab.png', 
    #        path = output)
    
    # pcReg
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
    pc.idx <- which(p.adj.value < tol)
    max.pc <- max(pc.idx)
    use.pca.data <- pca.data[,1:max.pc]
    # calculate pcReg
    sdev <- object@dr$pca@sdev[1:max.pc]
    var <- sdev^2 / sum(sdev^2)
    pcReg.var <- sum((var * r.squared[1:max.pc])[p.adj.value[1:max.pc] < tol])
    
    # k-means 
    batch.factor <- as.factor(object@meta.data$batch)
    fit.km <- kmeans(use.pca.data, centers = length(unique(batch.factor)))
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
        kBET(use.pca.data, batch.init, k0 = k0, knn = knn, plot = F)
    RR.kBET <- batch.estimate$summary$kBET.observed[1] -
        batch.estimate$summary$kBET.expected[1]
    AR.kBET <- 1 - (batch.estimate$summary$kBET.observed[1] -
        batch.estimate$summary$kBET.expected[1])
    
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
    norm.mixent <- mix.entropy/max.entropy
    
    # ARI
    ARI <- adjustedRandIndex(batch.init, batch.km)
    
    # NMI
    NMI <- mutinfo(batch.init, batch.km) / 
        sqrt(entropy(batch.init) * entropy(batch.km))
    
    # 
    
}



