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
    facLoc = 0.04, facScale = 0, num.batch = 3, num.cell = 200,
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
    object <- NormalizeData(object, verbose = F)
    # object <- FindVariableFeatures(
    #     object = object, nfeatures = 5000, verbose = F)
    
    object <- FindVariableFeatures(
        object = object, selection.method = 'mvp',
        mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf),
        verbose = F)
    # VariableFeaturePlot(object = object)
    length(object@assays$RNA@var.features)
    # object <- ScaleData(
    #     object, feature = object@assays$RNA@var.features, 
    #     use.umi = T, verbose = F)
    object <- ScaleData(
        object, feature = object@assays$RNA@var.features, 
        vars.to.regress = 'nCount_RNA', verbose = F)
    
    # PCA
    object <- RunPCA(
        object = object, npcs = 50, 
        feature = object@assays$RNA@var.features, verbose = F)
    pca.data <- object@reductions$pca@cell.embeddings
    # title.pca <- paste0('Variance parameter: ', facLoc)
    # ggplot.pca <- DimPlot(
    #     object, reductions = 'pca', group.by = "batch", 
    #     pt.size = 1)
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
    
    # UMAP
    object <- RunUMAP(
        object, reduction = 'pca', dims = 1:50, n.components = 2, verbose = F)
    umap.data <- object@reductions$umap@cell.embeddings
    # DimPlot(
    #     object, reductions = 'umap', group.by = "batch", pt.size = 1)
    
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
    use.pc <- which(p.adj.value < tol)
    # calculate pcReg
    sdev <- object@reductions$pca@stdev[use.pc]
    var <- sdev^2 / sum(sdev^2)
    pcReg.var <- sum((var * r.squared[use.pc])[p.adj.value[use.pc] < tol])
    
    # k-means 
    fit.km <- kmeans(umap.data, centers = length(unique(batch.factor)), 
                     iter.max = 100, nstart = length(unique(batch.factor)))
    batch.km <- as.numeric(fit.km$cluster)
    batch.init <- as.numeric(batch.factor)
    batch.km.new <- rep(0, length(batch.km))
    for (batch in unique(batch.init)) {
        sub.batch.km <- batch.km[batch.init == batch]
        sub.mode <- 
            as.numeric(names(table(sub.batch.km)))[
                table(sub.batch.km) == max(table(sub.batch.km))]
        batch.km.new[batch.km == sub.mode] <- batch
    }
    batch.km <- batch.km.new
    
    # sil
    dissimilar.dist <- dist(umap.data)
    sil.out <- silhouette(batch.init, dissimilar.dist)
    sil <- abs(summary(sil.out)$avg.width)
    
    # kBET
    k0 <- floor(mean(table(batch.init))) #neighbourhood size: mean batch size 
    knn.kBET <- get.knn(umap.data, k = k0, algorithm = 'cover_tree')
    batch.estimate <- 
        kBET(umap.data, batch.init, k0 = k0, knn = knn.kBET, plot = F)
    RR.kBET <- batch.estimate$summary$kBET.observed[1]
    # AR.kBET <- 1 - (batch.estimate$summary$kBET.observed[1] -
    #                     batch.estimate$summary$kBET.expected[1])
    
    # Entropy of batch mixing
    select.cells <- sample(dimnames(mtx.count)[2][[1]], 100)
    query.data <- umap.data[select.cells,]
    knn.mixent <- get.knnx(
        umap.data, query = query.data, k = 100, algorithm = 'cover_tree')
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
