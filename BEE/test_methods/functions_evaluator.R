library(Seurat)
library(cluster)
library(MASS)
library(reticulate)
use_python('python3')

# Pre-process data
preprocess.data <- function(input, batches, groups, num.pc = 30) {
    
    if (class(input) %in% c("data.frame", "matrix")) {
        mtx.count  <- input
        object <- CreateSeuratObject(mtx.count)
        object <- NormalizeData(object, verbose = F)
        object <- FindVariableFeatures(
            object = object, selection.method = 'mvp',
            mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf),
            verbose = F)
        # VariableFeaturePlot(object = object)
        length(object@assays$RNA@var.features)
    }
    if (class(input) == "Seurat") {
        object <- input
        if (is.null(object@commands$NormalizeData.RNA)) {
            object <- NormalizeData(object, verbose = F)
        }
        if (is.null(object@commands$FindVariableFeatures.RNA)) {
            object <- FindVariableFeatures(
                object = object, selection.method = 'mvp',
                mean.cutoff = c(0.0125, 3), dispersion.cutoff = c(0.5, Inf),
                verbose = F)
            # VariableFeaturePlot(object = object)
            length(object@assays$RNA@var.features)
        }
    }
    # creat Seurat object and normlize data
    object@meta.data$batch <- batches
    object@meta.data$group <- groups
    
    object <- ScaleData(
        object, feature = object@assays$RNA@var.features, verbose = F)
    
    # PCA
    object <- RunPCA(
        object = object, npcs = num.pc, 
        feature = object@assays$RNA@var.features, verbose = F)
    # ggplot.pca <- DimPlot(
    #     object, reductions = 'pca', group.by = "batch", 
    #     pt.size = 1)
    
    # UMAP
    object <- RunUMAP(
        object, reduction = 'pca', dims = 1:num.pc, n.components = 2, 
        verbose = F)
    # DimPlot(
    #     object, reductions = 'umap', group.by = "batch", pt.size = 1)
    
    return(object)
    
}

# Evaluate batch effect from two dimensions
evaluate.two.dims <- function(object) {
    
    # data
    pca.data <- object@reductions$pca@cell.embeddings
    umap.data <- object@reductions$umap@cell.embeddings
    
    # group
    group.factor <- as.factor(object@meta.data$group)
    group.init <- as.numeric(group.factor)
    
    # Silhouettes
    dissimilar.dist <- dist(umap.data)
    sil.out <- silhouette(group.init, dissimilar.dist)
    sil.group <- abs(summary(sil.out)$avg.width)
    
    # batch
    batch.factor <- as.factor(object@meta.data$batch)
    batch.init <- as.numeric(batch.factor)
    groups <- unique(group.init)
    df.group <- data.frame()
    for (group in groups) {
        sub.pca.data <- pca.data[group.init == group,]
        sub.batch.init <- batch.init[group.init == group]
        sub.batch.factor <- batch.factor[group.init == group]
        weight <- length(sub.batch.init)
        if (length(unique(sub.batch.init)) == 1) {
            df.group <- rbind(df.group,
                              data.frame(group = group, weight = 0,
                                         ldaReg = 0))
            next()
        }
        
        # LDA
        # remove collinear
        if (kappa(cor(sub.pca.data), exact = T) > 10000) {
            # pre.lda.data <- sub.pca.data
            # PCA <- import("sklearn.decomposition")
            # pca <- PCA$PCA()
            # pre.lda.data <-
            #     pca$fit(pre.lda.data)$transform(pre.lda.data)
            pre.lda.data <- 
                prcomp(sub.pca.data)$x[, 1:(dim(sub.pca.data)[2] - 1)]
            while (kappa(cor(pre.lda.data), exact = T) > 10000) {
                pre.lda.data <- pre.lda.data[, 1:(dim(pre.lda.data)[2] - 1)]
            }
        } else {
            pre.lda.data <- sub.pca.data
        }
        # print(group)
        # print(kappa(cor(pre.lda.data), exact = T))
        lda.input <- as.data.frame(pre.lda.data)
        lda.input$batch <- sub.batch.init
        fit.lda <- lda(formula = batch ~ ., data = lda.input, method = 'mve')
        lda.data <- pre.lda.data %*% fit.lda$scaling
        p.value <- c()
        r.squared <- c()
        for (i in 1:dim(lda.data)[2]) {
            fit <- summary(lm(lda.data[, i] ~ sub.batch.factor))
            r.squared <- c(r.squared, fit$r.squared)
            p.value <- c(p.value, 
                         1 - pf(fit$fstatistic[1], 
                                fit$fstatistic[2], 
                                fit$fstatistic[3]))
        }
        # p.adj.value <- p.adjust(p.value)
        # calculate ldaReg
        svd <- fit.lda$svd
        var <- svd^2 / sum(svd^2)
        ldaReg <- sum(var * r.squared)
        df.group <- rbind(df.group,
                          data.frame(group = group, weight = weight,
                                     ldaReg = ldaReg))
    }
    if (sum(df.group$weight) == 0) {
        print('Error: All batches only have one kinds of cell type')
    }
    ldaReg.batch <- sum((df.group$ldaReg*df.group$weight)/sum(df.group$weight))
    ldaReg.batch <- 1 - ldaReg.batch
    harmonic.mean <- 2*sil.group*ldaReg.batch/(sil.group + ldaReg.batch)
    
    output <- data.frame(cell.distance = sil.group, 
                         batch.effect.factor = ldaReg.batch,
                         harmonic.mean = harmonic.mean)
    
    return(output)
    
}


# Evaluate batch effect from two dimensions
evaluate.sil <- function(object) {
    
    # data
    umap.data <- object@reductions$umap@cell.embeddings
    
    # group
    group.factor <- as.factor(object@meta.data$group)
    group.init <- as.numeric(group.factor)
    
    # Silhouettes
    dissimilar.dist <- dist(umap.data)
    sil.out <- silhouette(group.init, dissimilar.dist)
    sil.group <- abs(summary(sil.out)$avg.width)
    
    # batch
    batch.factor <- as.factor(object@meta.data$batch)
    batch.init <- as.numeric(batch.factor)
    groups <- unique(group.init)
    df.group <- data.frame()
    for (group in groups) {
        sub.umap.data <- umap.data[group.init == group,]
        sub.batch.init <- batch.init[group.init == group]
        weight <- length(sub.batch.init)
        if (length(unique(sub.batch.init)) == 1) {
            df.group <- rbind(df.group,
                              data.frame(group = group, weight = 0,
                                         sub.sil = 0))
            next()
        }
        
        # Silhouettes
        dissimilar.dist <- dist(sub.umap.data)
        sil.out <- silhouette(sub.batch.init, dissimilar.dist)
        sub.sil <- abs(summary(sil.out)$avg.width)
        df.group <- rbind(df.group,
                          data.frame(group = group, weight = weight,
                                     sub.sil = sub.sil))
    }
    if (sum(df.group$weight) == 0) {
        print('Error: All batches only have one kinds of cell type')
    }
    sil.batch <- sum((df.group$sub.sil*df.group$weight)/sum(df.group$weight))
    sil.batch <- 1 - sil.batch
    harmonic.mean <- 2*sil.group*sil.batch/(sil.group + sil.batch)
    
    output <- data.frame(cell.distance = sil.group, 
                         batch.effect.factor = sil.batch,
                         harmonic.mean = harmonic.mean)
    
    return(output)
    
}


