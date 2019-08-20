library(Seurat)
library(cluster)
library(MASS)

# Pre-process data
preprocess.data <- function(mtx.count, batches, groups, num.pc = 30) {
    
    # creat Seurat object and normlize data
    object <- CreateSeuratObject(mtx.count)
    object@meta.data$batch <- batches
    object@meta.data$group <- groups
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
        
        # LDA
        lda.input <- as.data.frame(sub.pca.data)
        lda.input$batch <- sub.batch.init
        fit.lda <- lda(formula = batch ~ .,data = lda.input, method = 'mve')
        lda.data <- sub.pca.data %*% fit.lda$scaling
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
    ldaReg.batch <- sum((df.group$ldaReg*df.group$weight)/sum(df.group$weight))
    ldaReg.batch <- 1 - ldaReg.batch
    harmonic.mean <- 2*sil.group*ldaReg.batch/(sil.group + ldaReg.batch)
    
    output <- data.frame(cell.distance = sil.group, 
                         batch.effect.factor = ldaReg.batch,
                         harmonic.mean = harmonic.mean)
    
    return(output)
    
}


