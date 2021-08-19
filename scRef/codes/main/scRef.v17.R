########################################
#Author: Yu Zhang
#Email: zhang_yu18@fudan.edu.cn
#######################################

.get_overlap_genes <- function(exp_sc_mat, exp_ref_mat) {
    exp_ref_mat <- as.data.frame(exp_ref_mat)
    exp_sc_mat <- as.data.frame(exp_sc_mat)
    # get overlap genes
    exp_sc_mat <- exp_sc_mat[order(rownames(exp_sc_mat)),]
    exp_ref_mat <- exp_ref_mat[order(rownames(exp_ref_mat)),]
    gene_sc <- rownames(exp_sc_mat)
    gene_ref <- rownames(exp_ref_mat)
    gene_over <- gene_sc[which(gene_sc %in% gene_ref)]
    exp_sc_mat <- exp_sc_mat[gene_over,]
    exp_ref_mat <- exp_ref_mat[gene_over,]
    
    out.overlap <- list()
    out.overlap$exp_sc_mat <- exp_sc_mat
    out.overlap$exp_ref_mat <- exp_ref_mat
    out.overlap$gene_over <- gene_over
    return(out.overlap)
    
}


.get_high_variance_genes <- function(exp_ref_mat, num.genes = 2000, type_ref = 'sum-counts') {
    if (type_ref == 'sum-counts') {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref@assays$RNA@data <- exp_ref_mat
    }
    seurat.Ref <- FindVariableFeatures(
        seurat.Ref,
        selection.method = "vst",
        nfeatures = num.genes,
        verbose = F
    )
    return(VariableFeatures(seurat.Ref))
    
}


.one_multinomial <- function(i, exp_sc_mat, exp_ref_mat, colname_ref, 
                             verbose, print_step) {
    delta <- 0.5
    Refprob <- function(exp_sc, exp_ref) {
        log_p_sc_given_ref <- dmultinom(x = exp_sc, log = T, prob = exp_ref)
        return(log_p_sc_given_ref)
    }
    #################
    exp_sc <- as.array(exp_sc_mat[, i])
    log_p_sc_given_ref_list <- numeric(length = length(colname_ref))
    j = 1
    while (j <= length(colname_ref)) {
        exp_ref <- as.array(exp_ref_mat[, j])
        #####
        exp_ref[which(exp_ref == 0)] <- delta * min(exp_ref[which(exp_ref > 0)])
        #####
        log_p_sc_given_ref <- Refprob(exp_sc, exp_ref)
        log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
        j = j + 1
    }
    ################################
    if (verbose) {
        if (i %% print_step == 1) {
            print(i)
        }
    }
    
    return(log_p_sc_given_ref_list)
    
}


.get_log_p_sc_given_ref <- function(exp_sc_mat, exp_ref_mat, CPU = 4, 
                                    print_step = 100, gene_overlap = FALSE, 
                                    verbose = FALSE) {
    library(parallel, verbose = F)
    
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name

    #Step 1. get overlapped genes
    if (gene_overlap) {
        if (verbose) {
            print('Gene number of exp_sc_mat:')
            print(nrow(exp_sc_mat))
            print('Gene number of exp_ref_mat:')
            print(nrow(exp_ref_mat))
        }
        out.overlap <- .get_overlap_genes(exp_sc_mat, exp_ref_mat)
        exp_sc_mat <- out.overlap$exp_sc_mat
        exp_ref_mat <- out.overlap$exp_ref_mat
        if (verbose) {
            print('Number of overlapped genes:')
            print(nrow(exp_sc_mat))
        }
    }
    
    ###############
    colname_sc <- colnames(exp_sc_mat)
    colname_ref <- colnames(exp_ref_mat)
    
    #Step 2. calculate prob
    cl = makeCluster(CPU, outfile = '')
    RUN <- parLapply(
        cl = cl,
        1:length(exp_sc_mat[1,]),
        .one_multinomial,
        exp_sc_mat = exp_sc_mat,
        exp_ref_mat = exp_ref_mat,
        colname_ref = colname_ref,
        verbose = verbose,
        print_step = print_step
    )
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for (log_p_sc_given_ref_list in RUN) {
        LOG_P_SC_GIVEN_REF <-
            cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
    }
    #######################################
    rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
    colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
    ######2019.02.16 start ######
    LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <- min(LOG_P_SC_GIVEN_REF)
    ######2019.02.16 end ######
    return(LOG_P_SC_GIVEN_REF)
    
}


.one_get_corr <- function(barcode, exp_sc_mat, exp_ref_mat, colname_ref, 
                          method, verbose, print_step) {
    # calculate a correlation
    exp_sc <- as.array(exp_sc_mat[, barcode])
    log_p_sc_given_ref_list <- numeric(length = length(colname_ref))
    j <- 1
    while (j <= length(colname_ref)) {
        exp_ref <- as.array(exp_ref_mat[, j])
        if (method == 'kendall') {
            log_p_sc_given_ref <- cor.fk(exp_sc, exp_ref)
        } else {
            if (method == 'cosine') {
                log_p_sc_given_ref <-
                    sum(exp_sc * exp_ref) / sqrt(sum(exp_sc ^ 2) * sum(exp_ref ^ 2))
            } else {
                log_p_sc_given_ref <- cor(exp_sc, exp_ref, method = method)
            }
        }
        log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
        j <- j + 1
    }
    ################################
    if (verbose) {
        if (i %% print_step == 1) {
            print(i)
        }
    }
    
    gc()
    return(log_p_sc_given_ref_list)
    
}


.get_cor <- function(exp_sc_mat, exp_ref_mat, method = 'kendall', CPU = 4, 
                      print_step = 100, gene_overlap = FALSE, verbose = FALSE){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    
    #################
    library(parallel, verbose = F)
    #Step 1. get overlapped genes
    if (gene_overlap) {
        if (verbose) {
            print('Gene number of exp_sc_mat:')
            print(nrow(exp_sc_mat))
            print('Gene number of exp_ref_mat:')
            print(nrow(exp_ref_mat))
        }
        out.overlap <- .get_overlap_genes(exp_sc_mat, exp_ref_mat)
        exp_sc_mat <- out.overlap$exp_sc_mat
        exp_ref_mat <- out.overlap$exp_ref_mat
        if (verbose) {
            print('Number of overlapped genes:')
            print(nrow(exp_sc_mat))
        }
    }
    
    ###############
    colname_sc <- colnames(exp_sc_mat)
    colname_ref <- colnames(exp_ref_mat)
    exp_sc_mat <- as.matrix(exp_sc_mat)
    exp_ref_mat <- as.matrix(exp_ref_mat)
    #######################################
    
    #Step 2. calculate corr
    cl <- makeCluster(CPU, outfile = '')
    clusterEvalQ(cl, library(pcaPP))
    RUN <- parLapply(
        cl = cl,
        dimnames(exp_sc_mat)[[2]],
        .one_get_corr,
        exp_sc_mat = exp_sc_mat,
        exp_ref_mat = exp_ref_mat,
        colname_ref = colname_ref,
        method = method, 
        verbose = verbose, 
        print_step = print_step
    )
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF <- c()
    for (log_p_sc_given_ref_list in RUN) {
        LOG_P_SC_GIVEN_REF <-
            cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
    }
    #######################################
    rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
    colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
    ######2019.02.16 start ######
    LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <- -1
    ######2019.02.16 end ######
    
    return(LOG_P_SC_GIVEN_REF)
}


.get_tag_max <- function(P_REF_GIVEN_SC) {
    RN <- rownames(P_REF_GIVEN_SC)
    CN <- colnames(P_REF_GIVEN_SC)
    TAG <- cbind(CN, rep('NA', length(CN)))
    i <- 1
    while (i <= length(CN)) {
        this_rn_index <- which(P_REF_GIVEN_SC[, i] == max(P_REF_GIVEN_SC[, i]))[1]
        TAG[i, 2] <- RN[this_rn_index]
        i <- i + 1
    }
    colnames(TAG) <- c('cell_id', 'tag')
    return(TAG)
}


.get_tag_second <- function(P_REF_GIVEN_SC) {
    RN <- rownames(P_REF_GIVEN_SC)
    CN <- colnames(P_REF_GIVEN_SC)
    TAG <- cbind(CN, rep('NA', length(CN)))
    i <- 1
    while (i <= length(CN)) {
        TAG[i, 2] <- RN[order(P_REF_GIVEN_SC[,i], decreasing = T)[2]]
        i <- i + 1
    }
    colnames(TAG) <- c('cell_id', 'tag')
    return(TAG)
}


.get_tag_min <- function(P_REF_GIVEN_SC){
  RN=rownames(P_REF_GIVEN_SC)
  CN=colnames(P_REF_GIVEN_SC)
  TAG=cbind(CN,rep('NA',length(CN)))
  i=1
  while(i<=length(CN)){
    this_rn_index=which(P_REF_GIVEN_SC[,i] == min(P_REF_GIVEN_SC[,i]))[1]
    TAG[i,2]=RN[this_rn_index]
    i=i+1
  }
  colnames(TAG)=c('cell_id','tag')
  return(TAG)
}


.generate_ref <- function(exp_sc_mat, TAG, min_cell = 1, M = 'SUM', 
                          refnames = FALSE ){
    M <- M
    # print(M)
    min_cell <- min_cell
    refnames <- refnames
    exp_sc_mat <- exp_sc_mat
    TAG <- TAG
    NewRef <- c()
    TAG[, 2] <- as.character(TAG[, 2])
    if (refnames == FALSE) {
        refnames <- names(table(TAG[, 2]))
    }
    else{
        refnames <- refnames
    }
    outnames <- c()
    for (one in refnames) {
        this_col <- which(TAG[, 2] == one)
        if (length(this_col) >= min_cell) {
            outnames <- c(outnames, one)
            if (length(this_col) > 1) {
                if (M == 'SUM') {
                    this_new_ref <- apply(exp_sc_mat[, this_col], 1, sum)
                } else{
                    this_new_ref <- apply(exp_sc_mat[, this_col], 1, mean)
                }
            }
            else{
                this_new_ref <- exp_sc_mat[, this_col]
            }
            NewRef <- cbind(NewRef, this_new_ref)
        }
    }
    if (is.null(dim(NewRef))) {
        return(NULL)
    }
    rownames(NewRef) <- rownames(exp_sc_mat)
    colnames(NewRef) <- outnames
    # if (length(NewRef[1, ]) == 1) {
    #     NewRef <- cbind(NewRef[, 1], NewRef[, 1])
    #     rownames(NewRef) <- rownames(exp_sc_mat)
    #     colnames(NewRef) <- c(outnames, outnames)
    # }
    return(NewRef)
}


.getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE, 
                        block = NULL, p_adj = "fdr", fpkm = T) {
    # limma function
    if (is.null(esetm)) {
        cat(
            "esetm: gene expression matrix",
            "group: factor: \"c\"/\"d\"",
            "pair: TRUE/FALSE*",
            "block: e.g.1 2 2 1 if paired; blank if not",
            "p_adj: p.adjust, fdr* ",
            "fpkm: TRUE/FALSE*",
            sep = "\n"
        )
    } else{
        library(limma, verbose = F)
        if (pair) {
            design <- model.matrix( ~ block + group)
        } else{
            design <- model.matrix( ~ group)
        }
        suppressWarnings(fit <- lmFit(esetm, design))
        if (fpkm) {
            suppressWarnings(fit <- eBayes(fit, trend = T, robust = T))
        } else{
            suppressWarnings(fit <- eBayes(fit))
        }
        x <- topTable(fit, number = nrow(esetm), adjust.method = p_adj, 
                     coef = "group2")
        x <- x[!is.na(row.names(x)), ]
        x <- x[!duplicated(row.names(x)), ]
        return(x)
    }
}


.imoprt_outgroup <- function(out.group = 'MCA', use.RUVseq = T) {
    library(Seurat, verbose = F)
    if (class(out.group)[1] %in% c("data.frame", "matrix")) {
        df.out.group <- out.group
    } else {
        if (class(out.group)[1] == "character") {
            if (out.group %in% c("MCA", "HCA")) {
                if (out.group == "MCA") {
                    file.out.group <- './CellAtlas/MCA.txt'
                } else {
                    file.out.group <- './CellAtlas/HCA.txt'
                }
            } else {
                if (file.exists(out.group)) {
                    file.out.group <- out.group
                }
            }
            df.out.group <- 
                read.table(file.out.group, header = T, row.names = 1, 
                           sep = '\t', check.names = F)
        } else {
            stop('Error: incorrect input of outgroup')
        }
    }
    seurat.out.group <- 
        CreateSeuratObject(counts = df.out.group, project = "out.group", 
                           min.cells = 1, min.features = 5000)
    seurat.out.group <- 
        NormalizeData(seurat.out.group, normalization.method = "LogNormalize", 
                      scale.factor = 1e6, verbose = F)
    if (use.RUVseq) {
        # get stably expression genes
        seurat.out.group <- FindVariableFeatures(
            seurat.out.group, selection.method = "mvp",
            mean.cutoff = c(3, Inf), dispersion.cutoff = c(-0.05, 0.05), verbose = F)
    }

    return(seurat.out.group)
    
}


.find_markers_manual <- function(exp_ref_mat, ref_MCA_names, type_ref = 'sum-counts', 
                                 out.group = 'MCA', topN = NULL) {
    # check parameters
    if (!is.null(topN)) {
        topN <- topN
    } else {
        stop('Error in finding markers: provide incorrect parameters')
    }
    
    # regard MCA as reference of DEG
    seurat.MCA <- .imoprt_outgroup(out.group)
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    
    # transform count to fpm
    if (type_ref == 'sum-counts') {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
        exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
    }
    
    # overlap genes
    out.overlap <- .get_overlap_genes(fpm.MCA, exp_ref_mat)
    fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
    exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
    # print('Number of overlapped genes:')
    # print(nrow(exp_ref_mat))
    
    cell.MCA <- dimnames(fpm.MCA)[[2]]
    cell.ref <- dimnames(exp_ref_mat)[[2]]
    
    # match MCA names with ref names
    row.names(ref_MCA_names) <- ref_MCA_names$ref.name
    MCA.ref.names <- ref_MCA_names[cell.ref, 'MCA.name']
    
    # combat
    library(sva, verbose = F)
    mtx.in <- cbind(fpm.MCA, exp_ref_mat)
    names(mtx.in) <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
    batch <- c(rep(1, dim(fpm.MCA)[2]), rep(2, dim(exp_ref_mat)[2]))
    cov.cell <- c(cell.MCA, MCA.ref.names)
    mod <- model.matrix( ~ as.factor(cov.cell))
    mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T, ref.batch = 1)
    mtx.combat <- scale(mtx.combat)
    
    cells <- names(exp_ref_mat)
    
    list.cell.genes <- list()
    for (i in 1:length(cell.ref)) {
        cell <- cell.ref[i]
        vec.cell <- mtx.combat[, paste0('Ref.', cell)]
        vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, 0.85)]
        vec.ref <- exp_ref_mat[, cell]
        vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, 0.85)]
        # high expression genes
        genes.high <- intersect(names(vec.cell.high), names(vec.ref.high))
        # diff in reference
        cells <- c(setdiff(cell.ref, cell), cell)
        bool.cell <- cells
        bool.cell[bool.cell != cell] <- '1'
        bool.cell[bool.cell == cell] <- '2'
        res.limma.ref <- .getDEgeneF(exp_ref_mat[, cells], bool.cell)
        bool.cell <- c()
        genes.ref <-
            row.names(res.limma.ref[((res.limma.ref$P.Value < 0.1) &
                                         (res.limma.ref$logFC > 0.5)),])
        use.genes <- intersect(genes.high, genes.ref)
        
        mtx.combat.use <- mtx.combat[, cov.cell != MCA.ref.names[i]]
        mtx.limma <- cbind(mtx.combat.use, vec.cell)
        bool.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
        res.limma <- .getDEgeneF(mtx.limma, bool.cell)
        res.limma.high <- res.limma[use.genes,]
        res.limma.high <- res.limma.high[res.limma.high$logFC > 0,]
        res.limma.high <- res.limma.high[order(res.limma.high$P.Value),]
        
        genes.diff <- row.names(res.limma.high)[1:topN]
        list.cell.genes[[cell]] <- genes.diff
        
    }
    
    out <- list()
    out[['list.cell.genes']] <- list.cell.genes
    out[['exp_ref_mat']] <- exp_ref_mat
    return(out)
    
}


.find_markers_auto <- function(exp_ref_mat, out.group, list.localNeg = NULL,
                               type_ref = 'sum-counts', use.RUVseq = T, topN = 100) {
    library(parallel, verbose = F)
    library(Seurat, verbose = F)
    # check parameters
    if (!is.null(topN)) {
        topN <- topN
    } else {
        stop('Error in finding markers: provide incorrect parameters')
    }
    
    # transform count to fpm
    if (type_ref == 'sum-counts') {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
        exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm', 'bulk')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
    }
    
    ###### regard a outgroup (e.g. MCA/HCA) as reference of DEG
    seurat.MCA <- out.group
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    
    # overlap genes
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    out.overlap <- .get_overlap_genes(fpm.MCA, exp_ref_mat)
    fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
    exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
    if (use.RUVseq) {
        gene_overlap <- out.overlap$gene_over
        SEG.MCA <- VariableFeatures(seurat.MCA)
        gene.constant <- intersect(gene_overlap, SEG.MCA)
    }
    # print('Number of overlapped genes:')
    # print(nrow(exp_ref_mat))
    
    cell.MCA <- dimnames(fpm.MCA)[[2]]
    cell.ref <- dimnames(exp_ref_mat)[[2]]
    # use RUVseq to remove batch effect
    mtx.in <- cbind(fpm.MCA, exp_ref_mat)
    names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
    dimnames(mtx.in)[[2]] <- names.mix
    if (use.RUVseq) {
        library(RUVSeq, verbose = F)
        seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=1, isLog = T)
        mtx.combat <- seqRUVg$normalizedCounts
    } else {
        mtx.combat <- mtx.in
    }
    
    
    # topN <- 100
    list.cell.genes <- list()
    for (i in 1:length(cell.ref)) {
        cell <- cell.ref[i]
        # print(cell)
        vec.cell <- mtx.combat[, paste0('Ref.', cell)]
        vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, 0.7)]
        vec.ref <- exp_ref_mat[, cell]
        vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, 0.7)]
        # high expression genes
        genes.high <- intersect(names(vec.cell.high), names(vec.ref.high))
        # diff in reference
        cells <- c(setdiff(cell.ref, cell), cell)
        bool.cell <- cells
        bool.cell[bool.cell != cell] <- '1'
        bool.cell[bool.cell == cell] <- '2'
        res.limma.ref <- .getDEgeneF(exp_ref_mat[genes.high, cells], bool.cell)
        res.limma.ref <- res.limma.ref[res.limma.ref$logFC > 0,]
        genes.ref <- row.names(res.limma.ref[(res.limma.ref$P.Value < 0.1),])
        if (length(genes.ref) < (2*topN)) {
            res.limma.ref <- res.limma.ref[order(res.limma.ref$P.Value),]
            genes.ref <- row.names(res.limma.ref)[1:(2*topN)]
        }
        use.genes <- intersect(genes.high, genes.ref)
        
        # diff in negative reference
        if (!is.null(list.localNeg)) {
            count.neg <- list.localNeg[[cell]]
            if (!is.null(count.neg)) {
                seurat.neg <- CreateSeuratObject(counts = list.localNeg[[cell]], project = "Neg")
                seurat.neg <- NormalizeData(seurat.neg,normalization.method = "LogNormalize",
                                            scale.factor = 1e6, verbose = F)
                mat.ref.neg <- as.matrix(seurat.neg@assays$RNA@data)
                mat.neg <- cbind(mat.ref.neg[use.genes,], exp_ref_mat[use.genes, cell])
                bool.neg.cell <- c(rep('1', dim(mat.ref.neg)[2]), '2')
                res.limma.neg <- .getDEgeneF(mat.neg, bool.neg.cell)
                res.limma.neg <- res.limma.neg[res.limma.neg$logFC > 0,]
                genes.neg <- row.names(res.limma.neg[(res.limma.neg$P.Value < 0.05),])
                if (length(genes.neg) < (topN)) {
                    res.limma.neg <- res.limma.neg[order(res.limma.neg$P.Value),]
                    genes.neg <- row.names(res.limma.neg)[1:(topN)]
                }
                use.genes <- intersect(use.genes, genes.neg)
            }
        }
        
        # diff in Atlas
        if (length(use.genes) > topN) {
            mtx.combat.use <- mtx.combat[, names.mix != paste0('Ref.', cell)]
            mtx.limma <- cbind(mtx.combat.use, vec.cell)
            bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
            res.limma.MCA <- .getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
            res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
            res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
            genes.diff <- row.names(res.limma.MCA)[1:topN]
        } else {
            genes.diff <- use.genes
        }

        list.cell.genes[[cell]] <- genes.diff
        
    }
    
    out <- list()
    out[['list.cell.genes']] <- list.cell.genes
    out[['exp_ref_mat']] <- exp_ref_mat
    return(out)
    
}


.find_markers_auto_add_neg <- function(exp_ref_mat, seurat.out.group, list.localNeg = NULL,
                                       type_ref = 'sum-counts', use.RUVseq = T, 
                                       base.topN = 50, percent.high.exp = 0.8) {
    library(parallel, verbose = F)
    library(Seurat, verbose = F)
    # check parameters
    if (!is.null(base.topN)) {
        base.topN <- base.topN
    } else {
        stop('Error in finding markers: provide incorrect parameters')
    }
    
    # transform count to fpm
    if (type_ref == 'sum-counts') {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
        exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
    }
    if (type_ref %in% c('fpkm', 'tpm', 'rpkm', 'bulk')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
    }
    
    ###### regard a outgroup (e.g. MCA/HCA) as reference of DEG
    seurat.MCA <- seurat.out.group
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    
    # overlap genes
    fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data)
    out.overlap <- .get_overlap_genes(fpm.MCA, exp_ref_mat)
    fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
    exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
    if (use.RUVseq) {
        gene_overlap <- out.overlap$gene_over
        SEG.MCA <- VariableFeatures(seurat.MCA)
        gene.constant <- intersect(gene_overlap, SEG.MCA)
    }
    # print('Number of overlapped genes:')
    # print(nrow(exp_ref_mat))
    
    cell.MCA <- dimnames(fpm.MCA)[[2]]
    cell.ref <- dimnames(exp_ref_mat)[[2]]
    # use RUVseq to remove batch effect
    mtx.in <- cbind(fpm.MCA, exp_ref_mat)
    names.mix <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
    dimnames(mtx.in)[[2]] <- names.mix
    if (use.RUVseq) {
        library(RUVSeq, verbose = F)
        seqRUVg <- RUVg(as.matrix(mtx.in), gene.constant, k=1, isLog = T)
        mtx.combat <- seqRUVg$normalizedCounts
    } else {
        mtx.combat <- mtx.in
    }
    mtx.combat.use <- mtx.combat[, paste0('MCA.', cell.MCA)]
    
    
    topN <- base.topN
    list.cell.genes <- list()
    for (i in 1:length(cell.ref)) {
        cell <- cell.ref[i]
        # print(cell)
        vec.cell <- mtx.combat[, paste0('Ref.', cell)]
        vec.cell.high <- vec.cell[vec.cell > quantile(vec.cell, percent.high.exp)]
        vec.ref <- exp_ref_mat[, cell]
        vec.ref.high <- vec.ref[vec.ref > quantile(vec.ref, percent.high.exp)]
        # high expression genes
        genes.high <- intersect(names(vec.cell.high), names(vec.ref.high))
        
        # diff in reference
        cells <- c(setdiff(cell.ref, cell), cell)
        bool.cell <- cells
        bool.cell[bool.cell != cell] <- '1'
        bool.cell[bool.cell == cell] <- '2'
        res.limma.ref <- .getDEgeneF(exp_ref_mat[genes.high, cells], bool.cell)
        res.limma.ref <- res.limma.ref[res.limma.ref$logFC > 0,]
        genes.ref <- row.names(res.limma.ref[(res.limma.ref$P.Value < 0.1),])
        if (length(genes.ref) < (3*topN)) {
            res.limma.ref <- res.limma.ref[order(res.limma.ref$P.Value),]
            genes.ref <- row.names(res.limma.ref)[1:(3*topN)]
        }
        use.genes <- intersect(genes.high, genes.ref)
        
        # diff in negative reference
        if (!is.null(list.localNeg)) {
            count.neg <- list.localNeg[[cell]]
            if (!is.null(count.neg)) {
                seurat.neg <- CreateSeuratObject(counts = list.localNeg[[cell]], project = "Neg")
                seurat.neg <- NormalizeData(seurat.neg,normalization.method = "LogNormalize",
                                            scale.factor = 1e6, verbose = F)
                mat.ref.neg <- as.matrix(seurat.neg@assays$RNA@data)
                mat.neg <- cbind(mat.ref.neg[use.genes,], exp_ref_mat[use.genes, cell])
                bool.neg.cell <- c(rep('1', dim(mat.ref.neg)[2]), '2')
                res.limma.neg <- .getDEgeneF(mat.neg, bool.neg.cell)
                res.limma.neg <- res.limma.neg[res.limma.neg$logFC > 0,]
                genes.neg <- row.names(res.limma.neg[(res.limma.neg$P.Value < 0.05),])
                if (length(genes.neg) < (topN)) {
                    res.limma.neg <- res.limma.neg[order(res.limma.neg$P.Value),]
                    genes.neg <- row.names(res.limma.neg)[1:(topN)]
                }
                if (length(genes.neg) > (3*topN)) {
                    res.limma.neg <- res.limma.neg[order(res.limma.neg$P.Value),]
                    genes.neg <- row.names(res.limma.neg)[1:(3*topN)]
                }
                use.genes <- intersect(use.genes, genes.neg)
            }
        }
        
        # diff in Atlas
        if (length(use.genes) > topN) {
            mtx.limma <- cbind(mtx.combat.use, vec.cell)
            bool.atlas.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
            res.limma.MCA <- .getDEgeneF(mtx.limma[use.genes, ], bool.atlas.cell)
            res.limma.MCA <- res.limma.MCA[res.limma.MCA$logFC > 0,]
            genes.diff <- row.names(res.limma.MCA[(res.limma.MCA$P.Value < 0.001),])
            if (length(genes.diff) < (topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(topN)]
            }
            if (length(genes.diff) > (3*topN)) {
                res.limma.MCA <- res.limma.MCA[order(res.limma.MCA$P.Value),]
                genes.diff <- row.names(res.limma.MCA)[1:(3*topN)]
            }
        } else {
            genes.diff <- use.genes
        }
        
        list.cell.genes[[cell]] <- genes.diff
        # 
        # j = 1
        # while (j < i) {
        #     cell1 <- cell.ref[j]
        #     genes1 <- list.cell.genes[[cell1]]
        #     cell2 <- cell.ref[i]
        #     genes2 <- list.cell.genes[[cell2]]
        #     gene.intersect <- intersect(genes1, genes2)
        #     genes1.mod <- setdiff(genes1, gene.intersect)
        #     genes2.mod <- setdiff(genes2, gene.intersect)
        #     if (length(genes1.mod) > 30) {
        #         list.cell.genes[[cell1]] <- genes1.mod
        #     } else {
        #         warning(cell1, ' and ',  cell2, ": might be confused!")
        #     }
        #     if (length(genes2.mod) > 30) {
        #         list.cell.genes[[cell2]] <- genes2.mod
        #     } else {
        #         warning(cell1, ' and ',  cell2, ": might be confused!")
        #     }
        #     j = j + 1
        # }
        
    }
    
    out <- list()
    out[['list.cell.genes']] <- list.cell.genes
    out[['exp_ref_mat']] <- exp_ref_mat
    return(out)
    
}


.cluster_sc <- function(exp_sc_mat, cluster.num.pc = 75, cluster.resolution = 3) {
    # cluster
    library(Seurat, verbose = F)
    # data preparing
    seurat.unlabeled <- CreateSeuratObject(counts = exp_sc_mat)
    seurat.unlabeled <-
        NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize",
                      scale.factor = 10000, verbose = F)
    # print(seurat.unlabeled@assays$RNA@data[1:6,1:6])
    seurat.unlabeled <-
        FindVariableFeatures(seurat.unlabeled, selection.method = "vst",
                             nfeatures = 2000, verbose = F)
    seurat.unlabeled <- ScaleData(seurat.unlabeled, verbose = F)
    
    # PCA
    seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = cluster.num.pc, verbose = F)
    
    # cluster
    seurat.unlabeled <-
        FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:cluster.num.pc,
                      nn.eps = 0.5, verbose = F)
    seurat.unlabeled <-
        FindClusters(seurat.unlabeled, resolution = cluster.resolution, 
                     n.start = 10, n.iter = 100, verbose = F)
    
    out.cluster <-
        data.frame(
            cluster.id = as.character(seurat.unlabeled@meta.data$seurat_clusters),
            # original.label = seurat.unlabeled@meta.data$original.label,
            row.names = dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
        )
    return(out.cluster)
    
}


.cluster_increase_speed <- function(exp_sc_mat, df.cluster, cluster.cell = 5, CPU = 5) {
    library(parallel, verbose = F)
    sc.genes <- row.names(exp_sc_mat)
    df.cluster <- as.matrix(df.cluster)
    cluster.ids <- as.character(unique(df.cluster))
    cluster.cell <- cluster.cell
    
    .merge.one.cluster <- function(cluster.id) {
        # merge cells in one cluster
        cell.ids <- names(df.cluster[df.cluster[, 'cluster.id'] == cluster.id,])
        sub.exp <- exp_sc_mat[, cell.ids]
        # print(dim(sub.exp))
        
        sub.seurat <- CreateSeuratObject(counts = sub.exp)
        sub.seurat <-
            NormalizeData(sub.seurat, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = F)
        # print(sub.seurat@assays$RNA@data[1:6,1:6])
        sub.seurat <-
            FindVariableFeatures(sub.seurat, selection.method = "disp",
                                 nfeatures = 1000, verbose = F)
        # print('1')
        sub.seurat <- ScaleData(sub.seurat, 
                                features = VariableFeatures(sub.seurat), 
                                verbose = F)
        # print(head(sub.seurat@assays$RNA@var.features))
        # print('2')
        
        # PCA
        num.cell <- dim(sub.exp)[2]
        if (num.cell < 5) {
            stop('Error: You should provide a smaller resolution!')
        }
        sub.seurat <- RunPCA(sub.seurat, npcs = min(10, (num.cell-2)), verbose = F)
        sub.pca <- sub.seurat@reductions$pca@cell.embeddings
        
        num.clusters <- ceiling(num.cell / cluster.cell)
        if (num.clusters < 2) {
            sub.dict <- data.frame(
                cell.id = dimnames(sub.pca)[[1]],
                cluster.level1 = rep(cluster.id, num.cell),
                cluster.level2 = rep('1', num.cell)
            )
        } else {
            res.cluster <-
                kmeans(sub.pca, centers = num.clusters, nstart = 20, iter.max = 50)
            # names
            sub.dict <- data.frame(
                cell.id = names(res.cluster$cluster),
                cluster.level1 = rep(cluster.id, num.cell),
                cluster.level2 = res.cluster$cluster
            )
        }
        sub.dict$cluster.merge.id <-
            paste(sub.dict$cluster.level1, sub.dict$cluster.level2, sep = '-')
        row.names(sub.dict) <- sub.dict$cell.id
        
        # merge expression profile
        tag.in <- sub.dict[, c('cell.id', 'cluster.merge.id')]
        sub.exp.merge <- .generate_ref(sub.exp, tag.in)
        sub.exp.merge <- sub.exp.merge[sc.genes, ]
        # print(class(sub.exp.merge))
        if (class(sub.exp.merge)[1] %in% c('numeric', 'integer')) {
            sub.exp.merge <- as.data.frame(sub.exp.merge)
            names(sub.exp.merge) <- unique(sub.dict$cluster.merge.id)
        }
        
        sub.out <- list()
        sub.out$sub.dict <- sub.dict
        sub.out$sub.exp.merge <- sub.exp.merge
        return(sub.out)
        
    }
    # num.cpu <- floor(length(cluster.ids) / 10)
    
    # split dataset
    cl.input <- list()
    cl <- makeCluster(CPU, outfile = '')
    clusterExport(cl, '.generate_ref')
    clusterEvalQ(cl, library(Seurat))
    out.par <- parLapply(
        cl = cl,
        cluster.ids,
        .merge.one.cluster
        # exp_sc_mat = exp_sc_mat,
        # df.cluster = df.cluster,
        # cluster.cell = cluster.cell
    )
    stopCluster(cl)
    
    df.dict <- data.frame(stringsAsFactors = F)
    df.exp.merge <- data.frame(stringsAsFactors = F)
    i = 1
    for (sub.out in out.par) {
        df.dict <- rbind(df.dict, sub.out$sub.dict)
        if (i == 1) {
            df.exp.merge <- sub.out$sub.exp.merge
        } else {
            df.exp.merge <- cbind(df.exp.merge, sub.out$sub.exp.merge)
        }
        i = i + 1
    }
    
    out.merge <- list()
    out.merge$df.dict <- df.dict
    out.merge$df.exp.merge <- df.exp.merge
    return(out.merge)
    
}


.one_confirm_label <- function(barcode, exp_sc_mat, df.tag, list.cell.genes, 
                               method.test = 'wilcox.test') {
    expression.barcode <- exp_sc_mat[, barcode]
    bool.mark.gene <- rep(1, dim(exp_sc_mat)[1])
    cell <- df.tag[barcode, 'scRef.tag']
    genes.marker <- list.cell.genes[[cell]]
    bool.mark.gene[dimnames(exp_sc_mat)[[1]] %in% genes.marker] <- 2
    test.in <- cbind(expression.barcode, bool.mark.gene)
    test.in <- as.data.frame(test.in)
    names(test.in) <- c('expression.level', 'factor.mark.gene')
    test.in$factor.mark.gene <- as.factor(test.in$factor.mark.gene)
    if (method.test == 'ks.test') {
        vec.other <- test.in[test.in$factor.mark.gene == 1, 'expression.level']
        vec.marker <- test.in[test.in$factor.mark.gene == 2, 'expression.level']
        out.test <-
            ks.test(
                x = vec.other, y = vec.marker,
                data = test.in,
                alternative = 'greater'
            )
    } else {
        if (method.test == 'wilcox.test') {
            out.test <-
                wilcox.test(
                    formula = expression.level ~ factor.mark.gene,
                    data = test.in,
                    digits.rank = 0,
                    alternative = 'less'
                )
        } else {
            stop('Error: Please provide correct test method!')
        }
    }
    pvalue <- max(out.test$p.value, 1e-200)
    # library(ggplot2)
    # # ggplot(test.in, aes(x = expression.level, color = factor.mark.gene)) +
    # #     geom_line(stat = 'density') +
    # #     # ylim(0,0.5) +
    # #     xlim(0, 25)
    # test.in$marker.gene <- factor(as.numeric(test.in$factor.mark.gene), levels = c(1, 2),
    #                               labels = c('Other genes', 'Marker genes'))
    # plot.ecdf <- ggplot(test.in, aes(x = expression.level, color = marker.gene)) +
    #     stat_ecdf(size = 1.5) +
    #     scale_color_manual(breaks = c('Other genes', 'Marker genes'),
    #                        values = c('black', '#F8766D')) +
    #     ylim(0, 1) +
    #     xlim(0, 25) +
    #     labs(x = 'Counts', y = 'Empirical cumulative distribution function',
    #          color = '') +
    #     theme(panel.grid = element_blank(),
    #           panel.background = element_rect(fill='transparent', color='gray'),
    #           legend.key=element_rect(fill='transparent', color='transparent'),
    #           axis.text.x = element_text(size = 12),
    #           axis.text.y = element_text(size = 12),
    #           axis.title.x = element_text(size = 15),
    #           axis.title.y = element_text(size = 12),
    #           legend.text = element_text(size = 15))
    # ggsave(plot = plot.ecdf, path = '/home/zy/scRef/figure', filename = 'ecdf.png',
    #        units = 'cm', height = 10, width = 15)
    gc()
    return(data.frame(pvalue = pvalue, row.names = barcode))
    
}


.confirm_label <- function(exp_sc_mat, list.cell.genes, scRef.tag, CPU = 10) {
    library(parallel, verbose = F)
    # confirm label 
    exp_sc_mat <- as.matrix(exp_sc_mat)
    df.tag <- as.data.frame(scRef.tag)
    row.names(df.tag) <- df.tag$cell_id
    df.tag$cell_id <- NULL
    names(df.tag) <- 'scRef.tag'

    cl = makeCluster(CPU, outfile = '')
    RUN <- parLapply(
        cl = cl,
        row.names(df.tag),
        .one_confirm_label,
        exp_sc_mat = exp_sc_mat,
        df.tag = df.tag,
        list.cell.genes = list.cell.genes
    )
    stopCluster(cl)
    df.pval = data.frame()
    for (single.pval in RUN) {
        df.pval = rbind(df.pval, single.pval)
    }
    
    meta.tag <- merge(df.tag, df.pval, by = 'row.names')
    row.names(meta.tag) <- meta.tag$Row.names
    meta.tag$Row.names <- NULL
    meta.tag$log10Pval <- -log10(meta.tag$pvalue)
    
    gc()
    
    return(meta.tag)
    
}


.cutoff_GMM <- function(df.tags.in, num_cluster = 6, floor.cutoff = 5, ceiling.cutoff = 20, 
                        opt.strict = T) {
    library(mclust, verbose = F)
    cells <- unique(df.tags.in$scRef.tag)
    vec.cutoff <- c()
    for (i in 1:length(cells)) {
        cell <- cells[i]
        print(cell)
        df.sub <- df.tags.in[df.tags.in$scRef.tag == cell, ]
        if (length(unique(df.sub$log10Pval)) < 4) {
            vec.cutoff <- c(vec.cutoff, (ceiling.cutoff + floor.cutoff)/2)
            next()
        }
        model <- densityMclust(df.sub$log10Pval, 
                               G = min(num_cluster, length(unique(df.sub$log10Pval))-2), verbose = F)
        cluster.mean <- model$parameters$mean
        names(cluster.mean) <- as.character(1:length(cluster.mean))
        cluster.sd <- sqrt(model$parameters$variance$sigmasq)
        if (length(cluster.sd) == 1) {
            cluster.sd <- rep(sqrt(model$parameters$variance$sigmasq), length(cluster.mean))
        }
        names(cluster.sd) <- names(cluster.mean)
        all.clusters <- names(cluster.mean)
        cluster.unknown <- names(cluster.mean[cluster.mean < floor.cutoff])
        if (length(cluster.unknown) > 0) {
            unknown.right <- cluster.unknown[length(cluster.unknown)]
            cutoff.unknown <- 
                cluster.mean[unknown.right] + max(3 * cluster.sd[unknown.right], 4)
            for (cluster in setdiff(all.clusters, cluster.unknown)) {
                sub.mean <- cluster.mean[cluster]
                if (sub.mean < cutoff.unknown) {
                    cluster.unknown <- c(cluster.unknown, cluster)
                } else {
                    break()
                }
            }
        }
        cluster.known <- setdiff(all.clusters, cluster.unknown)
        cluster.final <- c(cluster.known[length(cluster.known)])
        if (length(cluster.known) > 1) {
            for (j in rev(1:(length(cluster.known) - 1))) {
                cluster <- cluster.known[j]
                logp.j.right <- cluster.mean[cluster] + max(3 * cluster.sd[cluster] ,15)
                known.left <- cluster.mean[cluster.known[j+1]]
                # print(cluster)
                # print(logp.j.right)
                # print(known.left)
                if (logp.j.right >= known.left) {
                    cluster.final <- c(cluster.final, cluster)
                } else {
                    if (sum(model$classification %in% cluster.final) < 5) {
                        cluster.final <- c(cluster.final, cluster)
                    } else {
                        break()
                    }
                }
            }
        }
        cluster.known <- rev(cluster.final)
        if (length(cluster.known) == 0) {
            vec.cutoff <- c(vec.cutoff, ceiling.cutoff)
            next()
        }
        df.sub$classification <- model$classification
        if (opt.strict) {
            sub.cutoff <- cluster.mean[cluster.known[1]]
            # sub.cutoff <- median(df.sub[df.sub$classification == cluster.known[1], 'log10Pval'])
        } else {
            sub.cutoff <- 
                min(df.sub[df.sub$classification %in% cluster.known, 'log10Pval'])
        }
        sub.cutoff <- max(sub.cutoff, floor.cutoff)
        sub.cutoff <- min(ceiling.cutoff, sub.cutoff)
        vec.cutoff <- c(vec.cutoff, sub.cutoff)
    }
    names(vec.cutoff) <- cells
    return(vec.cutoff)
}


.cutoff_GMM_add_neg <- function(df.tags.in, num_cluster = 6, cutoff.neg = 3, 
                                cutoff.pos = 5, ceiling.cutoff = 30, opt.strict = T) {
    library(mclust, verbose = F)
    cells <- as.character(unique(df.tags.in$scRef.tag))
    vec.cutoff <- c()
    vec.neg.cutoff <- c()
    for (i in 1:length(cells)) {
        cell <- cells[i]
        print(cell)
        df.sub <- df.tags.in[df.tags.in$scRef.tag == cell, ]
        if (length(unique(df.sub$log10Pval)) < 4) {
            vec.cutoff <- c(vec.cutoff, (ceiling.cutoff + cutoff.pos)/2)
            vec.neg.cutoff <- c(vec.neg.cutoff, neg.cutoff)
            next()
        }
        if (is.null(num_cluster)) {
            model <- densityMclust(df.sub$log10Pval, verbose = F)
        } else {
            model <- densityMclust(df.sub$log10Pval, 
                                   G = min(num_cluster, length(unique(df.sub$log10Pval))-2), 
                                   verbose = F)
        }
        cluster.mean <- model$parameters$mean
        names(cluster.mean) <- as.character(1:length(cluster.mean))
        cluster.sd <- sqrt(model$parameters$variance$sigmasq)
        if (length(cluster.sd) == 1) {
            cluster.sd <- rep(sqrt(model$parameters$variance$sigmasq), length(cluster.mean))
        }
        names(cluster.sd) <- names(cluster.mean)
        df.sub$classification <- model$classification
        
        # negative cutoff
        cluster.neg <- names(cluster.mean[cluster.mean < cutoff.neg])
        if (length(cluster.neg) == 0) {
            neg.cutoff <- NA
        } else {
            sub.neg.cutoff <- 
                max(df.sub[df.sub$classification %in% cluster.neg, 'log10Pval'])
            percent.neg <- 
                sum(df.sub$log10Pval < sub.neg.cutoff) / sum(df.sub$log10Pval < sub.neg.cutoff + 5)
            if (percent.neg > 0.7) {
                neg.cutoff <- min(sub.neg.cutoff, cutoff.neg)
            } else {
                neg.cutoff <- NA
            }
        }
        vec.neg.cutoff <- c(vec.neg.cutoff, neg.cutoff)
        
        # positive cutoff
        all.clusters <- names(cluster.mean)
        cluster.unknown <- names(cluster.mean[cluster.mean < cutoff.pos])
        cluster.known <- setdiff(all.clusters, cluster.unknown)
        cluster.final <- c(cluster.known[length(cluster.known)])
        if (length(cluster.known) > 1) {
            for (j in rev(1:(length(cluster.known) - 1))) {
                cluster <- cluster.known[j]
                logp.j.right <- cluster.mean[cluster] + 20
                known.left <- cluster.mean[cluster.known[j + 1]]
                if (logp.j.right >= known.left) {
                    cluster.final <- c(cluster.final, cluster)
                } else {
                    if (sum(model$classification %in% cluster.final) < 5) {
                        cluster.final <- c(cluster.final, cluster)
                    } else {
                        break()
                    }
                }
            }
        }
        cluster.known <- rev(cluster.final)
        
        if (length(cluster.known) == 0) {
            vec.cutoff <- c(vec.cutoff, ceiling.cutoff)
            next()
        }
        if (opt.strict) {
            sub.cutoff <- cluster.mean[cluster.known[1]]
            # sub.cutoff <- median(df.sub[df.sub$classification == cluster.known[1], 'log10Pval'])
        } else {
            sub.cutoff <- 
                min(df.sub[df.sub$classification %in% cluster.known, 'log10Pval'])
        }
        sub.cutoff <- max(sub.cutoff, cutoff.pos)
        sub.cutoff <- min(ceiling.cutoff, sub.cutoff)
        vec.cutoff <- c(vec.cutoff, sub.cutoff)
    }
    names(vec.cutoff) <- cells
    names(vec.neg.cutoff) <- cells
    out.cutoff <- list()
    out.cutoff$vec.cutoff <- vec.cutoff
    out.cutoff$vec.neg.cutoff <- vec.neg.cutoff
    return(out.cutoff)
    
}


.combine_tags <- function(df.tags1, df.tags2) {
    # concat reference pval and local pval
    pvalue1 <- df.tags1[, c('scRef.tag', 'pvalue')]
    names(pvalue1) <- c('scRef.tag.1', 'pvalue.1')
    pvalue2 <- df.tags2[, c('scRef.tag', 'pvalue')]
    names(pvalue2) <- c('scRef.tag.2', 'pvalue.2')
    pvalue <- merge(pvalue1, pvalue2, by = 'row.names')
    row.names(pvalue) <- pvalue$Row.names
    pvalue$Row.names <- NULL
    
    # select more confident tag
    mtx.tag <- as.matrix(pvalue[, c('scRef.tag.1', 'scRef.tag.2')])
    mtx.pval <- as.matrix(pvalue[, c('pvalue.1', 'pvalue.2')])
    mtx.rank <- apply(mtx.pval, 1, rank, ties.method = "first")
    tag.final <-
        apply(as.array(1:dim(mtx.tag)[1]), 1, function(i) {
            mtx.tag[i, mtx.rank[1, i]]
        })
    pval.final <-
        apply(as.array(1:dim(mtx.pval)[1]), 1, function(i) {
            mtx.pval[i, mtx.rank[1, i]]
        })
    tag.final <-
        data.frame(
            scRef.tag = tag.final,
            pvalue = pval.final,
            row.names = dimnames(mtx.tag)[[1]],
            stringsAsFactors = F
        )
    
    OUT <- list()
    OUT$pvalue <- pvalue
    OUT$tag.final <- tag.final
    return(OUT)
    
}


SCREF <- function(exp_sc_mat, exp_ref_mat, exp_ref_label = NULL, 
                  identify_unassigned = T, single_round = F, corr_use_HVGene = T,
                  type_ref = 'sc-counts', out.group = 'MCA', use.RUVseq = T, 
                  topN = 50, percent.high.exp = 0.8, 
                  cluster.num.pc =50, cluster.resolution = 3, 
                  cluster.speed = T, cluster.cell = 10,
                  method1 = 'kendall', method2 = 'multinomial', 
                  # cutoff.1 = 'default', cutoff.2 = 'default', 
                  GMM.num_cluster = 6, GMM.neg_cutoff = 3, GMM.floor_cutoff = 5, GMM.ceiling_cutoff = 20,
                  threshold.recall = 0.2,
                  min_cell = 20, CPU = 4, mod = 'debug') {
    library(parallel, verbose = F)
    # check parameters
    if (!type_ref %in% c('sc-counts', 'sum-counts', 'fpkm', 'tpm', 'rpkm')) {
        stop('Error: inexistent input of reference data format')
    }
    cutoff.1 = 'default'
    cutoff.2 = 'default'
    
    time1 <- Sys.time()
    # get sum-counts format
    if (type_ref == 'sc-counts') {
        print('Sum single cell counts matrix:')
        label.in <- data.frame(cell_id = colnames(exp_ref_mat), 
                               tag = as.character(exp_ref_label))
        exp_ref_mat.sum <- .generate_ref(exp_ref_mat, label.in, M='SUM')
        exp_ref_mat <- exp_ref_mat.sum
        type_ref <- 'sum-counts'
    }

    # get overlap genes
    out.overlap <- .get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    gene_over <- out.overlap$gene_over
    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))
    
    seurat.out.group <- .imoprt_outgroup(out.group = out.group, use.RUVseq = use.RUVseq)
    
    if (identify_unassigned) {
        # find markers of cell types in reference
        print('Find marker genes of cell types in reference:')
        suppressMessages(
            out.markers <-
                .find_markers_auto_add_neg(
                    exp_ref_mat,
                    seurat.out.group,
                    type_ref = type_ref,
                    use.RUVseq = use.RUVseq,
                    base.topN = topN, 
                    percent.high.exp = percent.high.exp
                ))
        list.cell.genes <- out.markers[['list.cell.genes']]
        genes.ref <- dimnames(out.markers[['exp_ref_mat']])[[1]]
        
        # overlap genes
        gene.overlap <- intersect(gene_over, genes.ref)
        exp_sc_mat <- exp_sc_mat[gene.overlap, ]
        exp_ref_mat <- exp_ref_mat[gene.overlap, ]
        print('Number of overlapped genes:')
        print(nrow(exp_sc_mat))
        
        # cluster analysis
        print('Start clustering :')
        df.cluster <-
            .cluster_sc(exp_sc_mat,
                        cluster.num.pc = cluster.num.pc,
                        cluster.resolution = cluster.resolution)
        print('Clustering completed!')
        
        # speed calculation
        if (cluster.speed) {
            print('Speed calculation by clustering:')
            out.merge <-
                .cluster_increase_speed(exp_sc_mat, df.cluster, 
                                        cluster.cell = cluster.cell, CPU = CPU)
            df.dict <- out.merge$df.dict
            df.exp.merge <- out.merge$df.exp.merge
            min_cell <- ceiling(min_cell / cluster.cell)
        } else {
            df.exp.merge <- exp_sc_mat
        }
    } else {
        df.exp.merge <- exp_sc_mat
    }
    
    # rm(exp_sc_mat)
    gc()
    df.exp.merge <- as.matrix(df.exp.merge)

    print('First-round annotation:')
    print(method1)
    if (corr_use_HVGene) {
        HVG <- .get_high_variance_genes(exp_ref_mat, type_ref = type_ref)
        similarity.in <- df.exp.merge[HVG, ]
        ref.in <- exp_ref_mat[HVG, ]
    } else {
        similarity.in <- df.exp.merge
        ref.in <- exp_ref_mat
    }
    if (method1 != 'multinomial') {
        out1 <- .get_cor(similarity.in, ref.in, method = method1, CPU = CPU)
    } else {
        out1 <- .get_log_p_sc_given_ref(similarity.in, ref.in, CPU = CPU)
    }
    
    gc()
    
    tag1 <- .get_tag_max(out1)

    if (identify_unassigned) {
        print('Build local reference')
        df.tags1 <- .confirm_label(df.exp.merge, list.cell.genes, tag1, CPU = CPU)
        cell_ids <- colnames(df.exp.merge)
        df.tags1 <- df.tags1[cell_ids, ]
        df.tags1$cluster.merge.id <- row.names(df.tags1)
        df.tags.merge <- merge(df.tags1, df.dict, by = 'cluster.merge.id')
        row.names(df.tags.merge) <- df.tags.merge$cell.id
        df.tags.merge$cell.id <- NULL
        df.view1 <- merge(label_sc, df.tags.merge, by = 'row.names')
        View(df.view1)
        
        # select cutoff.1 automatitically
        if (cutoff.1 == 'default') {
            # df.cutoff.1 <- .cutoff_GMM(df.tags1, num_cluster = GMM.num_cluster, 
            #                            floor.cutoff = GMM.floor_cutoff, 
            #                            ceiling.cutoff = GMM.ceiling_cutoff)
            out.cutoff <- .cutoff_GMM_add_neg(df.tags1, num_cluster = GMM.num_cluster, 
                                              cutoff.neg = GMM.neg_cutoff, 
                                              cutoff.pos = GMM.floor_cutoff,
                                              ceiling.cutoff = GMM.ceiling_cutoff)
            df.cutoff.1 <- out.cutoff$vec.cutoff
            neg.cutoff.1 <- out.cutoff$vec.neg.cutoff
            print('Default cutoff: ')
            print(df.cutoff.1)
            
            if (single_round) {
                df.tags <- df.tags1
                df.tags$scRef.tag.1 <- as.character(df.tags$scRef.tag)
                df.tags$scRef.tag <- df.tags$scRef.tag.1
                for (cell in names(df.cutoff.1)) {
                    sub.cutoff <- df.cutoff.1[cell]
                    df.tags$scRef.tag[(df.tags$scRef.tag.1 == cell) & 
                                          (df.tags$log10Pval < sub.cutoff)] <- 'Unassigned'
                }
                
                if (cluster.speed) {
                    df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
                    df.cluster <- unique(df.cluster)
                    df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                                             row.names = df.cluster$cluster.merge.id,
                                             stringsAsFactors = F)
                }
                df.tags <- merge(df.tags, df.cluster, by = 'row.names')
                row.names(df.tags) <- df.tags$Row.names
                df.tags$Row.names <- NULL
                
                if (cluster.speed) {
                    df.tags$cluster.merge.id <- row.names(df.tags)
                    df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
                                           by = 'cluster.merge.id')
                    df.tags.merge$cluster.merge.id <- NULL
                    row.names(df.tags.merge) <- df.tags.merge$cell.id
                    df.tags.merge$cell.id <- NULL
                    df.tags <- df.tags.merge
                }
                
                # recall Unassigned
                df.tags$scRef.tag.pre.recall <- df.tags$scRef.tag
                cluster.ids <- unique(df.cluster$cluster.id)
                info.cluster <- data.frame(stringsAsFactors = F)
                for (cluster.id in cluster.ids) {
                    sub.tag.cluster <- df.tags[df.tags$cluster.id == cluster.id,]
                    table.tag <- table(sub.tag.cluster$scRef.tag.1)
                    table.tag <- table.tag[order(table.tag, decreasing = T)]
                    main.cell <- names(table.tag[1])
                    percent.main.cell <- table.tag[1] / dim(sub.tag.cluster)[1]
                    num.Unassigned <- sum(sub.tag.cluster$scRef.tag.pre.recall == 'Unassigned')
                    percent.Unassigned <- num.Unassigned / nrow(sub.tag.cluster)
                    if (percent.Unassigned < threshold.recall) {
                        if (percent.main.cell > (1- threshold.recall - 0.05)) {
                            df.tags[(df.tags$cluster.id == cluster.id) & 
                                        (df.tags$scRef.tag == 'Unassigned'), 'scRef.tag'] <-
                                rep(main.cell, num.Unassigned)
                        } else {
                            df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag'] <-
                                df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag.1']
                        }
                    }
                    info.cluster <- rbind(
                        info.cluster,
                        data.frame(
                            cluster.id = cluster.id,
                            percent.Unassigned = percent.Unassigned,
                            main.cell = main.cell,
                            percent.main.cell = percent.main.cell,
                            stringsAsFactors = F
                        )
                    )
                }
                
                df.combine <- df.tags[, c("scRef.tag", "log10Pval")]
                cell_ids <- colnames(exp_sc_mat)
                df.combine <- df.combine[cell_ids, ]
                
                gc()
                
                time2 <- Sys.time()
                time.scRef <- difftime(time2, time1, units = 'secs')
                output <- list()
                output$tag1 <- tag1
                output$out1 <- out1
                output$combine.out <- df.tags
                output$info.cluster <- info.cluster
                if (cluster.speed) {
                    output$dict.cluster <- df.dict
                }
                output$ref.markers <- list.cell.genes
                output$final.out <- df.combine
                output$run.time <- time.scRef
                
                print('Finish!')
                
                return(output)
                
            }
            
            select.barcode <- c()
            for (cell in names(df.cutoff.1)) {
                sub.cutoff <- df.cutoff.1[cell]
                sub.select <- df.tags1[df.tags1$scRef.tag == cell, ]
                sub.select <- sub.select[sub.select$log10Pval >= sub.cutoff, ]
                select.barcode <- c(select.barcode, row.names(sub.select))
            }
            list.localNeg <- list()
            for (cell in names(neg.cutoff.1)) {
                print(cell)
                sub.cutoff <- neg.cutoff.1[cell]
                if (is.na(sub.cutoff)) {
                    next()
                }
                sub.select <- df.tags1[df.tags1$scRef.tag == cell, ]
                sub.select <- sub.select[sub.select$log10Pval <= sub.cutoff, ]
                neg.barcode <- row.names(sub.select)
                if (length(neg.barcode) < 2) {
                    next()
                }
                neg.exp <- df.exp.merge[, cell_ids %in% neg.barcode]
                if (cluster.speed) {
                    neg.tag1 <- df.dict[df.dict[, 'cluster.merge.id'] %in% neg.barcode, 
                                        c('cluster.merge.id', 'cluster.level1')]
                    names(neg.tag1) <- c('cell_id', 'tag')
                    neg.tag1 <- unique(neg.tag1)
                } else {
                    neg.tag1 <- data.frame(cell_id = neg.barcode, 
                                           tag = df.cluster[neg.barcode,])
                }
                LocalNeg <- .generate_ref(neg.exp, neg.tag1,  min_cell = min_cell + 2)
                if (is.null(LocalNeg)) {next()}
                if (ncol(LocalNeg) < 2) {next()}
                list.localNeg[[cell]] <- LocalNeg
            }
        } else {
            # cutoff.1 <- 1e-20
            # select.df.tags <- df.tags1[df.tags1$log10Pval >= cutoff.1, ]
            # select.barcode <- row.names(select.df.tags)
        }
        select.exp <- df.exp.merge[, cell_ids %in% select.barcode]
        select.tag1 <- tag1[tag1[, 'cell_id'] %in% select.barcode, ]
        LocalRef <- .generate_ref(select.exp, select.tag1,  min_cell = min_cell)

    } else {
        if (single_round) {
            df.combine <- as.data.frame(tag1)
            row.names(df.combine) <- df.combine$cell_id
            df.combine$cell_id <- NULL
            names(df.combine) <- 'scRef.tag'
            cell_ids <- colnames(exp_sc_mat)
            df.combine <- data.frame(scRef.tag = df.combine[cell_ids, ], row.names = cell_ids)
            
            time2 <- Sys.time()
            time.scRef <- difftime(time2, time1, units = 'secs')
            output <- list()
            output$tag1 <- tag1
            output$out1 <- out1
            output$final.out <- df.combine
            output$run.time <- time.scRef
            
            print('Finish!')
            
            return(output)
            
        }
        LocalRef <- .generate_ref(df.exp.merge, tag1, min_cell = min_cell)
    }
    print('Cell types in local reference:')
    print(dimnames(LocalRef)[[2]])
    
    #####
    gc()
    #####
    if (identify_unassigned) {
        # find local marker genes
        print('find local marker genes')
        out.markers <-
            .find_markers_auto_add_neg(
                LocalRef,
                seurat.out.group,
                list.localNeg = list.localNeg,
                use.RUVseq = use.RUVseq,
                base.topN = topN,
                percent.high.exp = percent.high.exp
            )
        local.cell.genes <- out.markers[['list.cell.genes']]
    }
    
    print('Second-round annotation:')
    print(method2)
    if (corr_use_HVGene) {
        HVG <- .get_high_variance_genes(LocalRef)
        similarity.in <- df.exp.merge[HVG, ]
        ref.in <- LocalRef[HVG, ]
    } else {
        similarity.in <- df.exp.merge
        ref.in <- LocalRef
    }
    if (method2 != 'multinomial') {
        out2 <- .get_cor(similarity.in, ref.in, method = method2, CPU = CPU)
    } else {
        out2 <- .get_log_p_sc_given_ref(similarity.in, ref.in, CPU = CPU)
    }
    tag2 <- .get_tag_max(out2)
    gc()
    
    if (identify_unassigned) {
        df.tags2 <- .confirm_label(df.exp.merge, local.cell.genes, tag2, CPU = CPU)
        cell_ids <- colnames(df.exp.merge)
        df.tags2 <- df.tags2[cell_ids, ]
        gc()

        print('Combine reference and local result:')
        # combine reference pval and local pval
        del.cells <- setdiff(colnames(exp_ref_mat), colnames(LocalRef))
        df.tags1$pvalue[df.tags1$scRef.tag %in% del.cells] <- 1
        out.combine <- .combine_tags(df.tags1, df.tags2)
        tag.final <- out.combine$tag.final
        # tag.final <- df.tags2
        pvalue <- out.combine$pvalue
       
        # modify tags and combine pval
        df.tags <- merge(pvalue, tag.final, by = 'row.names')
        row.names(df.tags) <- df.tags$Row.names
        df.tags$Row.names <- NULL
        
        df.tags$log10Pval <- -log10(df.tags$pvalue)
        
        if (cluster.speed) {
            df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
            df.cluster <- unique(df.cluster)
            df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                                     row.names = df.cluster$cluster.merge.id,
                                     stringsAsFactors = F)
        }
        
        # select cutoff.2 automatitically
        if (cutoff.2 == 'default') {
            out.cutoff.2 <- .cutoff_GMM_add_neg(df.tags, num_cluster = GMM.num_cluster, 
                                              cutoff.neg = GMM.neg_cutoff, 
                                              cutoff.pos = GMM.floor_cutoff,
                                              ceiling.cutoff = GMM.ceiling_cutoff, 
                                              opt.strict = F)
            df.cutoff.2 <- out.cutoff.2$vec.cutoff
            print('Default cutoff: ')
            print(df.cutoff.2)
            df.tags$scRef.tag.12 <- as.character(df.tags$scRef.tag)
            df.tags$scRef.tag <- df.tags$scRef.tag.12
            for (cell in names(df.cutoff.2)) {
                sub.cutoff <- df.cutoff.2[cell]
                df.tags$scRef.tag[(df.tags$scRef.tag.12 == cell) & 
                                      (df.tags$log10Pval < sub.cutoff)] <- 'Unassigned'
            }
        } else {
            # print('Cutoff: ')
            # print(cutoff.2)
            # df.tags$scRef.tag.12 <- as.character(df.tags$scRef.tag.12)
            # df.tags$scRef.tag <- df.tags$scRef.tag.12
            # df.tags$scRef.tag[df.tags$log10Pval < cutoff.2] <- 'Unassigned'
        }
        
        df.tags <- merge(df.tags, df.cluster, by = 'row.names')
        row.names(df.tags) <- df.tags$Row.names
        df.tags$Row.names <- NULL
        
        if (cluster.speed) {
            df.tags$cluster.merge.id <- row.names(df.tags)
            df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
                                   by = 'cluster.merge.id')
            df.tags.merge$cluster.merge.id <- NULL
            row.names(df.tags.merge) <- df.tags.merge$cell.id
            df.tags.merge$cell.id <- NULL
            df.tags <- df.tags.merge
        }
        
        # recall Unassigned and other cell types
        df.tags$scRef.tag.pre.recall <- df.tags$scRef.tag
        cluster.ids <- unique(df.cluster$cluster.id)
        info.cluster <- data.frame(stringsAsFactors = F)
        for (cluster.id in cluster.ids) {
            sub.tag.cluster <- df.tags[df.tags$cluster.id == cluster.id,]
            num.cell <- nrow(sub.tag.cluster)
            table.tag <- table(sub.tag.cluster$scRef.tag.pre.recall)
            table.tag <- table.tag[order(table.tag, decreasing = T)]
            main.cell <- names(table.tag[1])
            percent.main.cell <- table.tag[1] / dim(sub.tag.cluster)[1]
            if ((main.cell == 'Unassigned') & (length(table.tag) > 1)) {
                main.cell <- names(table.tag[2])
                percent.main.cell <- table.tag[2] / dim(sub.tag.cluster)[1]
            }
            num.Unassigned <- sum(sub.tag.cluster$scRef.tag.pre.recall == 'Unassigned')
            percent.Unassigned <- num.Unassigned / num.cell
            num.other <- sum(sub.tag.cluster$scRef.tag.pre.recall != main.cell)
            percent.other <- num.other / num.cell
            cell.other <- paste(setdiff(names(table.tag), c(main.cell, 'Unassigned')), 
                                collapse = '|')
            if (percent.Unassigned < threshold.recall) {
                if (percent.main.cell > (1 - threshold.recall - 0.2)) {
                    df.tags[(df.tags$cluster.id == cluster.id) & 
                                (df.tags$scRef.tag.pre.recall == 'Unassigned'), 'scRef.tag'] <-
                        rep(main.cell, num.Unassigned)
                } else {
                    df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag'] <-
                        df.tags[df.tags$cluster.id == cluster.id, 'scRef.tag.pre.recall']
                }
            }
            if (percent.main.cell > (0.8)) {
                df.tags[(df.tags$cluster.id == cluster.id) & 
                            (df.tags$scRef.tag.pre.recall != main.cell), 'scRef.tag'] <-
                    rep(main.cell, num.other)
            }
            info.cluster <- rbind(
                info.cluster,
                data.frame(
                    cluster.id = cluster.id,
                    num.cell = num.cell,
                    main.cell = main.cell,
                    percent.main.cell = percent.main.cell,
                    percent.other = percent.other,
                    percent.Unassigned = percent.Unassigned,
                    cell.other = cell.other,
                    stringsAsFactors = F
                )
            )
        }

        df.combine <- df.tags[, c("scRef.tag", "log10Pval")]
        cell_ids <- colnames(exp_sc_mat)
        df.combine <- df.combine[cell_ids, ]
        df.tags <- df.tags[cell_ids, ]
        
    } else {
        df.combine <- as.data.frame(tag2)
        row.names(df.combine) <- df.combine$cell_id
        df.combine$cell_id <- NULL
        names(df.combine) <- 'scRef.tag'
        cell_ids <- colnames(exp_sc_mat)
        df.combine <- data.frame(scRef.tag = df.combine[cell_ids, ], row.names = cell_ids)
    }
    
    # df.view <- merge(label_sc, df.tags, by = 'row.names')
    
    # df.view <- merge(df.combine, df.tags, by = 'row.names')
    # row.names(df.view) <- df.view$Row.names
    # df.view$Row.names <- NULL
    # df.view <- merge(label.filter, df.view, by = 'row.names')
    
    #####
    gc()
    #####
    time2 <- Sys.time()
    time.scRef <- difftime(time2, time1, units = 'secs')
    output <- list()
    output$tag1 <- tag1
    output$out1 <- out1
    output$tag2 <- tag2
    output$out2 <- out2
    if (identify_unassigned) {
        output$pvalue1 <- df.tags1
        output$pvalue2 <- df.tags2
        output$combine.out <- df.tags
        output$info.cluster <- info.cluster
        if (cluster.speed) {
            output$dict.cluster <- df.dict
        }
        output$cutoff.1 <- df.cutoff.1
        output$cutoff.neg.1 <- neg.cutoff.1
        output$cutoff.2 <- df.cutoff.2
        output$ref.markers <- list.cell.genes
        output$local.markers <- local.cell.genes
    }
    output$final.out <- df.combine
    output$run.time <- time.scRef
    if (mod == 'debug') {
        output$df.exp.merge <- df.exp.merge
    }
    
    print('Finish!')
    
    return(output)
    
}


supervised.UMAP <- function(mtx.in, labels) {
    library(Seurat)
    # use python in R
    library(reticulate)
    use_python('/home/zy/tools/anaconda3/bin/python3', required = T)
    # py_config()
    # import python package: umap
    print(paste0('Whether umap package is imported: ', py_module_available('umap')))
    
    # data preparing
    seurat.unlabeled <- CreateSeuratObject(counts = mtx.in)
    seurat.unlabeled <- NormalizeData(seurat.unlabeled, normalization.method = "LogNormalize", 
                                      scale.factor = 10000)
    seurat.unlabeled <- FindVariableFeatures(seurat.unlabeled, selection.method = "vst", nfeatures = 2000)
    seurat.unlabeled <- ScaleData(seurat.unlabeled)
    seurat.unlabeled@meta.data$label <- labels
    
    # PCA
    seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 75, verbose = F)
    mat.pca <- seurat.unlabeled@reductions$pca@cell.embeddings

    # transform character label to index
    vec.labels <- c(setdiff(unique(labels), 'Unassigned'), 'Unassigned')
    df.label.idx <- data.frame(label = vec.labels, idx = c(1:(length(vec.labels) - 1), -1))
    for (j in 1:dim(df.label.idx)[1]) {
        labels[labels == df.label.idx[j, 'label']] <- df.label.idx[j, 'idx']
    }
    
    # supervised UMAP
    umap <- import('umap')
    class.umap <- umap$UMAP()
    embedding <- class.umap$fit_transform(X = mat.pca, y = as.numeric(labels))
    dimnames(embedding)[[1]] <- dimnames(mat.pca)[[1]]
    umap.label <- CreateDimReducObject(embeddings = embedding, key = 'UMAP_')
    seurat.unlabeled@reductions$umap.label <- umap.label
    
    plot.out <- DimPlot(seurat.unlabeled, reduction = "umap.label", label = T, group.by = 'label')
    
    return(plot.out)
    
}
