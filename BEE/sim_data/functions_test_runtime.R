library(splatter)
library(scater)
library(Seurat)
library(kBET)
library(cluster)
library(FNN)
library(mclust)
library(MASS)

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
calculate_evaluate_runtime <- function(
    facLoc = 0.06, facScale = 0, num.batch = 2, 
    num.cell = 200, num.gene = 15000, num.pc = 30,
    seed.splatter = 1234) {
    
    # simulate
    if (length(num.cell) != 1) {
        if (length(num.cell) != num.batch) {
            print('Error: batch number mismatch')
        } else {
            batchCells <- num.cell
        }
    } else {
        batchCells = rep(num.cell, num.batch)
    }
    params <- newSplatParams(
        nGenes = num.gene, 
        batchCells = batchCells, 
        batch.facLoc = rep(facLoc, num.batch),
        batch.facScale = rep(facScale, num.batch),
        seed = seed.splatter)
    sim.data <- splatSimulate(params, verbose = F)
    mtx.count <- counts(sim.data)
    batches <- sim.data$Batch
    
    # creat Seurat object and normlize data
    object <- CreateSeuratObject(mtx.count)
    object@meta.data$batch <- batches
    object <- NormalizeData(object, verbose = F)

    object <- FindVariableFeatures(
        object = object, selection.method = 'mvp',
        mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf),
        verbose = F)
    # VariableFeaturePlot(object = object)
    length(object@assays$RNA@var.features)
    object <- ScaleData(
        object, feature = object@assays$RNA@var.features, 
        vars.to.regress = 'nCount_RNA', verbose = F)
    batch.factor <- as.factor(object@meta.data$batch)
    batch.init <- as.numeric(batch.factor)
    # print('scale')
    
    # PCA
    time1 <- Sys.time()
    object <- RunPCA(
        object = object, npcs = num.pc, 
        feature = object@assays$RNA@var.features, verbose = F)
    pca.data <- object@reductions$pca@cell.embeddings
    time2 <- Sys.time()
    time.pca <- difftime(time2, time1, units = 'secs') 
    
    # UMAP
    time3 <- Sys.time()
    object <- RunUMAP(
        object, reduction = 'pca', dims = 1:num.pc, n.components = 2, 
        verbose = F)
    umap.data <- object@reductions$umap@cell.embeddings
    time4 <- Sys.time()
    time.umap <- difftime(time4, time3, units = 'secs') 
    # print('pca&umap')
    
    # pcReg
    time5 <- Sys.time()
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
    pcReg <- sum((var * r.squared[use.pc])[p.adj.value[use.pc] < tol])
    time6 <- Sys.time()
    time.pcReg <- difftime(time6, time5, units = 'secs') + time.pca
    
    # k-means 
    time7 <- Sys.time()
    fit.km <- kmeans(umap.data, centers = length(unique(batch.factor)), 
                     iter.max = 200, nstart = 10*length(unique(batch.factor)))
    batch.km <- as.numeric(fit.km$cluster)
    batch.km.new <- rep(0, length(batch.km))
    for (batch in unique(batch.init)) {
        sub.batch.km <- batch.km[batch.init == batch]
        sub.mode <- 
            as.numeric(names(table(sub.batch.km)))[
                table(sub.batch.km) == max(table(sub.batch.km))]
        batch.km.new[batch.km == sub.mode] <- batch
    }
    batch.km <- batch.km.new
    time8 <- Sys.time()
    time.kmeans <- difftime(time8, time7, units = 'secs')
    
    # sil
    time9 <- Sys.time()
    dissimilar.dist <- dist(umap.data)
    sil.out <- silhouette(batch.init, dissimilar.dist)
    sil <- abs(summary(sil.out)$avg.width)
    time10 <- Sys.time()
    time.sil <- difftime(time10, time9, units = 'secs') + time.umap
    
    # kBET
    time11 <- Sys.time()
    k0 <- floor(mean(table(batch.init))) #neighbourhood size: mean batch size 
    knn.kBET <- get.knn(umap.data, k = k0, algorithm = 'cover_tree')
    batch.estimate <- 
        kBET(umap.data, batch.init, k0 = k0, knn = knn.kBET, plot = F)
    RR.kBET <- batch.estimate$summary$kBET.observed[1]
    time12 <- Sys.time()
    time.kBET <- difftime(time12, time11, units = 'secs') + time.umap
    # print('kBET')

    # Entropy of batch mixing
    time13 <- Sys.time()
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
    time14 <- Sys.time()
    time.mixent <- difftime(time14, time13, units = 'secs') + time.umap
    
    # ARI
    time15 <- Sys.time()
    ARI <- adjustedRandIndex(batch.init, batch.km)
    time16 <- Sys.time()
    time.ARI <- difftime(time16, time15, units = 'secs') + 
        time.kmeans + time.umap
    
    # NMI
    time17 <- Sys.time()
    NMI <- mutinfo(batch.init, batch.km) / 
        sqrt(entropy(batch.init) * entropy(batch.km))
    time18 <- Sys.time()
    time.NMI <- difftime(time18, time17, units = 'secs') + 
        time.kmeans + time.umap
    # print('ARI')
    
    # LDA
    time19 <- Sys.time()
    lda.input <- as.data.frame(pca.data)
    lda.input$batch <- batch.init
    fit.lda <- lda(formula = batch ~ .,data = lda.input, method = 'mve')
    lda.data <- pca.data %*% fit.lda$scaling
    p.value <- c()
    r.squared <- c()
    for (i in 1:dim(lda.data)[2]) {
        fit <- summary(lm(lda.data[, i] ~ batch.factor))
        r.squared <- c(r.squared, fit$r.squared)
        p.value <- c(p.value, 
                     1 - pf(fit$fstatistic[1], 
                            fit$fstatistic[2], 
                            fit$fstatistic[3]))
    }
    p.adj.value <- p.adjust(p.value)
    # calculate ldaReg
    svd <- fit.lda$svd
    var <- svd^2 / sum(svd^2)
    ldaReg <- sum(var * r.squared)
    time20 <- Sys.time()
    time.ldaReg <- difftime(time20, time19, units = 'secs') + time.pca
    # print('LDA')
    
    # output
    access.res <- 
        data.frame(pcReg = pcReg, sil = sil, RR.kBET = RR.kBET,
                   norm.mixent = norm.mixent, ARI = ARI, NMI = NMI, 
                   ldaReg = ldaReg)
    output <- data.frame(
        pcReg = time.pcReg, sil = time.sil, kBET = time.kBET, 
        mixent = time.mixent, ARI = time.ARI, NMI = time.NMI, 
        ldaReg = time.ldaReg
    )

    return(output)
    
}
