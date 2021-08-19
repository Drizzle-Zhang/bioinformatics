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


.get_log_p_sc_given_ref <-
  function(exp_sc_mat, exp_ref_mat, CPU = 4, print_step = 10, verbose = FALSE) {
    delta <- 0.5
    
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    ##################
    library(parallel)
    ##################
    Refprob <- function(exp_sc, exp_ref) {
      log_p_sc_given_ref <- dmultinom(x = exp_sc, log = T, prob = exp_ref)
      return(log_p_sc_given_ref)
    }
    #################
    #Step 1. get overlapped genes
    out.overlap <- .get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    colname_sc <- colnames(exp_sc_mat)
    colname_ref <- colnames(exp_ref_mat)
    #Step 2. calculate prob
    SINGLE <- function(i) {
      exp_sc <- as.array(exp_sc_mat[, i])
      #log_p_sc_given_ref_list=c()
      log_p_sc_given_ref_list <- rep(0, length(colname_ref))
      j = 1
      while (j <= length(colname_ref)) {
        exp_ref <- as.array(exp_ref_mat[, j])
        #####
        exp_ref[which(exp_ref == 0)] <- 
            delta * min(exp_ref[which(exp_ref > 0)])
        #####
        log_p_sc_given_ref <- Refprob(exp_sc, exp_ref)
        #log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
        log_p_sc_given_ref_list[j] <- log_p_sc_given_ref
        j = j + 1
      }
      ################################
      if (verbose) {
        if (i %% print_step == 1) {
          print(i)
        }
        return(log_p_sc_given_ref_list)
      }
    }
    #######################################
    cl = makeCluster(CPU, outfile = '')
    RUN = parLapply(cl = cl, 1:length(exp_sc_mat[1, ]), SINGLE)
    stopCluster(cl)
    #RUN = mclapply(1:length(colname_sc), SINGLE, mc.cores=CPU)
    LOG_P_SC_GIVEN_REF = c()
    for (log_p_sc_given_ref_list in RUN) {
      LOG_P_SC_GIVEN_REF <- cbind(LOG_P_SC_GIVEN_REF, log_p_sc_given_ref_list)
    }
    #######################################
    rownames(LOG_P_SC_GIVEN_REF) <- colname_ref
    colnames(LOG_P_SC_GIVEN_REF) <- colname_sc
    ######2019.02.16 start ######
    LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] <-
      min(LOG_P_SC_GIVEN_REF)
    ######2019.02.16 end ######
    return(LOG_P_SC_GIVEN_REF)
  }


.get_cor  <- function(exp_sc_mat, exp_ref_mat, method = 'kendall', CPU = 4, 
                      print_step = 100, verbose = FALSE){
    #method = "pearson", "kendall", "spearman"
    #exp_sc_mat: single-cell gene expression matrix; row is gene; col is sample; should have row name and col name
    #exp_ref_mat: reference gene expression matrix; row is gene; col is sample; should have row name and col name
    ##################
    if (verbose) {
        print('Gene number of exp_sc_mat:')
        print(nrow(exp_sc_mat))
        print('Gene number of exp_ref_mat:')
        print(nrow(exp_ref_mat))
    }
    
    #################
    library(parallel)
    #Step 1. get overlapped genes
    out.overlap <- .get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    colname_sc <- colnames(exp_sc_mat)
    colname_ref <- colnames(exp_ref_mat)
    ###############
    if (verbose) {
        print('Number of overlapped genes:')
        print(nrow(exp_sc_mat))
    }
    ###################
    #Step 2. calculate prob
    SINGLE <- function(i) {
        library('pcaPP')
        exp_sc <- as.array(exp_sc_mat[, i])
        #log_p_sc_given_ref_list=c()
        log_p_sc_given_ref_list <- rep(0, length(colname_ref))
        j <- 1
        while (j <= length(colname_ref)) {
            exp_ref <- as.array(exp_ref_mat[, j])
            #####
            #exp_ref[which(exp_ref==0)]=min(exp_ref[which(exp_ref>0)])
            #####
            #if(method=='rococo'){log_p_sc_given_ref=rococo(exp_sc,exp_ref)} else
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
            #log_p_sc_given_ref_list=c(log_p_sc_given_ref_list, log_p_sc_given_ref)
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
    #######################################
    cl <- makeCluster(CPU, outfile = '')
    RUN <- parLapply(cl = cl, 1:length(exp_sc_mat[1, ]), SINGLE)
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
    LOG_P_SC_GIVEN_REF[which(is.na(LOG_P_SC_GIVEN_REF))] = -1
    ######2019.02.16 end ######
    
    return(LOG_P_SC_GIVEN_REF)
}


.get_tag_max <- function(P_REF_GIVEN_SC){
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
    rownames(NewRef) <- rownames(exp_sc_mat)
    colnames(NewRef) <- outnames
    if (length(NewRef[1, ]) == 1) {
        NewRef <- cbind(NewRef[, 1], NewRef[, 1])
        rownames(NewRef) <- rownames(exp_sc_mat)
        colnames(NewRef) <- c(outnames, outnames)
    }
    return(NewRef)
}


# limma function
.getDEgeneF <- function(esetm = NULL, group = NULL, pair = FALSE, 
                        block = NULL, p_adj = "fdr", fpkm = T) {
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
        fit <- lmFit(esetm, design)
        if (fpkm) {
            fit <- eBayes(fit, trend = T, robust = T)
        } else{
            fit <- eBayes(fit)
        }
        x <-
            topTable(fit, number = nrow(esetm), adjust.method = p_adj, 
                     coef = "group2")
        x <- x[!is.na(row.names(x)), ]
        x <- x[!duplicated(row.names(x)), ]
        return(x)
    }
}


.imoprt_outgroup <- function(out.group = 'MCA', num.use.gene = NULL) {
    library(Seurat, verbose = F)
    if (class(out.group)[1] %in% c("data.frame", "matrix")) {
        df.out.group <- out.group
    } else {
        if (out.group %in% c("MCA", "HCA")) {
            if (out.group == "MCA") {
                file.out.group <- './CellAtlas/MCA.txt'
            } else {
                file.out.group <- './CellAtlas/HCA.txt'
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
                           min.cells = 1, min.features = 1000)
    seurat.out.group <- 
        NormalizeData(seurat.out.group, normalization.method = "LogNormalize", 
                      scale.factor = 1e6, verbose = F)
    if (!is.null(num.use.gene)) {
        seurat.out.group <- FindVariableFeatures(
            seurat.out.group, selection.method = "vst",
            nfeatures = num.use.gene, verbose = F)
        # VariableFeaturePlot(seurat.out.group)
        use.genes <- VariableFeatures(seurat.out.group)
        fpm.out.group <- as.matrix(seurat.out.group@assays$RNA@data[use.genes,])
    } else {
        fpm.out.group <- as.matrix(seurat.out.group@assays$RNA@data)
    }
    
    return(fpm.out.group)
    
}


.find_markers_manual <- function(exp_ref_mat, type_ref = 'count', out.group = 'MCA', 
                                 topN = NULL, cutoff.fc = NULL, cutoff.pval = NULL) {
  ###### regard MCA as reference of DEG
  fpm.MCA <- .imoprt_outgroup(out.group)
  
  # transform count to fpm
  if (type_ref == 'count') {
    exp_ref_mat <- as.data.frame(exp_ref_mat)
    coldata.ref <- DataFrame(row.names = names(exp_ref_mat))
    obj.DESeq.ref <- DESeqDataSetFromMatrix(countData = exp_ref_mat, colData = coldata.ref, 
                                            design = ~ 1)
    exp_ref_mat <- fpm(obj.DESeq.ref, robust = T)
    exp_ref_mat <- as.data.frame(exp_ref_mat)
  }
  
  # overlap genes
  fpm.MCA=fpm.MCA[order(rownames(fpm.MCA)),]
  exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
  gene_MCA=rownames(fpm.MCA)
  gene_ref=rownames(exp_ref_mat)
  gene_over= gene_MCA[which(gene_MCA %in% gene_ref)]
  fpm.MCA=fpm.MCA[which(gene_MCA %in% gene_over),]
  exp_ref_mat=exp_ref_mat[which(gene_ref %in% gene_over),]
  print('Number of overlapped genes:')
  print(nrow(exp_ref_mat))
  
  cell.MCA <- dimnames(fpm.MCA)[[2]]
  cell.ref <- names(exp_ref_mat)
  cell.overlap <- intersect(cell.MCA, cell.ref)
  # combat
  library(sva)
  mtx.in <- cbind(fpm.MCA, exp_ref_mat)
  names(mtx.in) <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
  batch <- c(rep(1, dim(fpm.MCA)[2]), rep(2, dim(exp_ref_mat)[2]))
  cov.cell <- c(cell.MCA, names(exp_ref_mat))
  mod <- model.matrix(~ as.factor(cov.cell))
  mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T, ref.batch = 1)
  mtx.combat <- scale(mtx.combat)
  
  cells <- names(exp_ref_mat)
  
  # cutoff.fc <- 1.5
  # cutoff.pval <- 0.05
  # topN <- 100
  # cutoff.fc <- NULL
  # cutoff.pval <- NULL
  list.cell.genes <- list()
  for (cell in cells) {
    vec.cell <- mtx.combat[, paste0('Ref.', cell)]
    exp.top10 <- quantile(vec.cell, 0.9)
    genes.high <- gene_over[vec.cell > exp.top10]
    mtx.combat.use <- mtx.combat[, !(dimnames(mtx.combat)[[2]] %in% c(paste0('MCA.', cell), paste0('Ref.', cell)))]
    mtx.limma <- cbind(mtx.combat.use, vec.cell)
    mtx.limma.in <- mtx.limma[genes.high,]
    bool.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
    res.limma <- .getDEgeneF(mtx.limma.in, bool.cell)
    if (is.null(topN) & (!is.null(cutoff.fc) & !is.null(cutoff.pval))) {
      df.diff <- res.limma[
        ((res.limma$logFC > cutoff.fc) & (res.limma$adj.P.Val < cutoff.pval)),]
      genes.diff <- row.names(df.diff)
    }
    if (!is.null(topN) & (is.null(cutoff.fc) & is.null(cutoff.pval))) {
      genes.diff <- row.names(res.limma)[1:topN]
    }
    if (!is.null(topN) & (!is.null(cutoff.fc) & !is.null(cutoff.pval))) {
      print('Error: provide too many parameters')
      return()
    }
    if (is.null(topN) & (is.null(cutoff.fc) & is.null(cutoff.pval))) {
      print('Error: provide too few parameters')
      return()
    }
    list.cell.genes[[cell]] <- genes.diff
  }
  
  out <- list()
  out[['list.cell.genes']] <- list.cell.genes
  out[['exp_ref_mat']] <- exp_ref_mat
  return(out)
  
}


.find_markers_auto1 <- function(exp_ref_mat, type = 'count', out.group = 'MCA', 
                                topN = NULL, cutoff.fc = NULL, cutoff.pval = NULL) {
  ###### regard MCA as reference of DEG
  file.MCA <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_combined_mouse_uniform.txt'
  # file.MCA.outer <- '/home/disk/scRef/MouseAtlas_SingleCell_Han2018/combinedMCA/MCA_concat_outer.txt'
  df.MCA <- read.table(file.MCA, header=T, row.names=1, sep='\t', check.name=F)
  library(Seurat, verbose = F)
  seurat.MCA <- CreateSeuratObject(counts = df.MCA, project = "MCA", 
                                   min.cells = 1, min.features = 1000)
  seurat.MCA <- NormalizeData(seurat.MCA, normalization.method = "LogNormalize", 
                              scale.factor = 1e6, verbose = F)
  seurat.MCA <- FindVariableFeatures(seurat.MCA, selection.method = "vst", 
                                     nfeatures = 10000)
  # VariableFeaturePlot(seurat.MCA)
  use.genes <- VariableFeatures(seurat.MCA)
  fpm.MCA <- as.matrix(seurat.MCA@assays$RNA@data[use.genes,])
  
  # transform count to fpm
  if (type == 'count') {
    seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref", 
                                     min.cells = 1, min.features = 1000)
    seurat.Ref <- NormalizeData(seurat.Ref, normalization.method = "LogNormalize", 
                                scale.factor = 1e6, verbose = F)
    exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
  }
  if (type == 'fpkm') {
    exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
  }
  
  # overlap genes
  fpm.MCA=fpm.MCA[order(rownames(fpm.MCA)),]
  exp_ref_mat=exp_ref_mat[order(rownames(exp_ref_mat)),]
  gene_MCA=rownames(fpm.MCA)
  gene_ref=rownames(exp_ref_mat)
  gene_over= gene_MCA[which(gene_MCA %in% gene_ref)]
  fpm.MCA=fpm.MCA[gene_over,]
  exp_ref_mat=exp_ref_mat[gene_over,]
  # print('Number of overlapped genes:')
  # print(nrow(exp_ref_mat))
  
  # generate cov
  out <- .get_cor(fpm.MCA, exp_ref_mat, method='pearson', CPU=4, print_step=10)
  ref.names <- dimnames(out)[[1]]
  MCA.names <- dimnames(out)[[2]]
  tag <- .get_tag_max(out)
  out.rank <- apply(out, 1, rank)
  vec.rank <- apply(tag, 1, function(x) {out.rank[x[1], x[2]]})
  vec.corr <- apply(tag, 1, function(x) {out[x[2], x[1]]})
  tag <- as.data.frame(tag)
  tag$rank <- vec.rank
  tag$corr <- vec.corr
  # rank cutoff
  rank.cutoff <- length(MCA.names) - length(ref.names)
  corr.cutoff <- min(tag[tag$rank > rank.cutoff, 'corr'])
  cor.max <- c()
  new.ref.names <- c()
  for (i in 1:dim(out)[1]) {
    sub.tag <- tag[tag$tag == ref.names[i],]
    max.rank <- max(sub.tag$rank)
    if (max.rank > rank.cutoff) {
      new.ref.names <- c(new.ref.names, sub.tag[sub.tag$rank == max.rank, 'cell_id'])
      cor.max <- c(cor.max, sub.tag[sub.tag$rank == max.rank, 'corr'])
    } else {
      sub.max <- max(out[ref.names[i],])
      sub.ref.name <- MCA.names[out[ref.names[i],] == sub.max]
      if (sub.max > corr.cutoff) {
        new.ref.names <- c(new.ref.names, sub.ref.name)
        cor.max <- c(cor.max, sub.max)
      } else {
        new.ref.names <- c(new.ref.names, 'NA')
        cor.max <- c(cor.max, sub.max)
      }
    }
  }
  print('Reference cells:')
  print(ref.names)
  print('Matched cells:')
  print(new.ref.names)
  print('Pearson coefficient:')
  print(cor.max)
  
  
  cell.MCA <- dimnames(fpm.MCA)[[2]]
  cell.ref <- dimnames(exp_ref_mat)[[2]]
  # cell.overlap <- intersect(cell.MCA, cell.ref)
  # combat
  library(sva, verbose = F)
  mtx.in <- cbind(fpm.MCA, exp_ref_mat)
  dimnames(mtx.in)[[2]] <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
  batch <- c(rep(1, dim(fpm.MCA)[2]), rep(2, dim(exp_ref_mat)[2]))
  cov.cell <- c(cell.MCA, new.ref.names)
  mod <- model.matrix(~ as.factor(cov.cell))
  mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T, ref.batch = 1)
  mtx.combat <- scale(mtx.combat)
  
  # cutoff.fc <- 1.5
  # cutoff.pval <- 0.05
  # topN <- 100
  # cutoff.fc <- NULL
  # cutoff.pval <- NULL
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
    res.limma.ref <- .getDEgeneF(exp_ref_mat[,cells], bool.cell)
    bool.cell <- c()
    genes.ref <- row.names(res.limma.ref[((res.limma.ref$P.Value < 0.1) & (res.limma.ref$logFC > 0.5)),])
    use.genes <- intersect(genes.high, genes.ref)
    
    mtx.combat.use <- mtx.combat[, cov.cell != new.ref.names[i]]
    mtx.limma <- cbind(mtx.combat.use, vec.cell)
    bool.cell <- as.factor(c(rep('1', dim(mtx.combat.use)[2]), '2'))
    res.limma <- .getDEgeneF(mtx.limma, bool.cell)
    res.limma.high <- res.limma[use.genes,]
    res.limma.high <- res.limma.high[res.limma.high$logFC > 0,]
    res.limma.high <- res.limma.high[order(res.limma.high$P.Value),]
    
    if (is.null(topN) & (!is.null(cutoff.fc) & !is.null(cutoff.pval))) {
      df.diff <- res.limma.high[
        ((res.limma.high$logFC > cutoff.fc) & (res.limma.high$adj.P.Val < cutoff.pval)),]
      genes.diff <- row.names(df.diff)
    }
    if (!is.null(topN) & (is.null(cutoff.fc) & is.null(cutoff.pval))) {
      genes.diff <- row.names(res.limma.high)[1:topN]
    }
    if (!is.null(topN) & (!is.null(cutoff.fc) & !is.null(cutoff.pval))) {
      print('Error: provide too many parameters')
      return()
    }
    if (is.null(topN) & (is.null(cutoff.fc) & is.null(cutoff.pval))) {
      print('Error: provide too few parameters')
      return()
    }
    list.cell.genes[[cell]] <- genes.diff
  }
  
  out <- list()
  out[['list.cell.genes']] <- list.cell.genes
  out[['exp_ref_mat']] <- exp_ref_mat
  return(out)
  
}


.find_markers_auto <- function(exp_ref_mat, type_ref = 'count', 
                               out.group = 'MCA', topN = NULL) {
    # check parameters
    if (!is.null(topN)) {
        topN <- topN
    } else {
        stop('Error in finding markers: provide incorrect parameters')
    }
    
    ###### regard a outgroup (e.g. MCA/HCA) as reference of DEG
    fpm.MCA <- .imoprt_outgroup(out.group)
    
    # transform count to fpm
    if (type_ref == 'count') {
        seurat.Ref <- CreateSeuratObject(counts = exp_ref_mat, project = "Ref")
        seurat.Ref <- NormalizeData(seurat.Ref,normalization.method = "LogNormalize",
                                    scale.factor = 1e6, verbose = F)
        exp_ref_mat <- as.matrix(seurat.Ref@assays$RNA@data)
    }
    if (type_ref %in% c('fpkm', 'tpm')) {
        exp_ref_mat <- as.matrix(log1p(exp_ref_mat))
    }
    
    # overlap genes
    out.overlap <- .get_overlap_genes(fpm.MCA, exp_ref_mat)
    fpm.MCA <- as.matrix(out.overlap$exp_sc_mat)
    exp_ref_mat <- as.matrix(out.overlap$exp_ref_mat)
    # print('Number of overlapped genes:')
    # print(nrow(exp_ref_mat))
    
    # generate cov
    out <- .get_cor(fpm.MCA, exp_ref_mat, method = 'pearson', CPU = 4)
    ref.names <- dimnames(out)[[1]]
    MCA.names <- dimnames(out)[[2]]
    new.ref.names <- c()
    for (i in 1:dim(out)[1]) {
        sub.max <- max(out[i,])
        # print(paste0('Reference cells: ', ref.names[i]))
        sel.names <-
            MCA.names[order(out[i,], decreasing = T)[1]]
        print(paste0(
            'Reference cells: ', ref.names[i],
            ' | ',
            'Matched cells: ', sel.names
        ))
        print(paste0('Pearson coefficient: ', round(sub.max, digits = 4)))
        new.ref.names <- c(new.ref.names, sel.names)
    }
    
    
    cell.MCA <- dimnames(fpm.MCA)[[2]]
    cell.ref <- dimnames(exp_ref_mat)[[2]]
    # cell.overlap <- intersect(cell.MCA, cell.ref)
    # combat
    library(sva, verbose = F)
    mtx.in <- cbind(fpm.MCA, exp_ref_mat)
    dimnames(mtx.in)[[2]] <- c(paste0('MCA.', cell.MCA), paste0('Ref.', cell.ref))
    batch <- c(rep(1, dim(fpm.MCA)[2]), rep(2, dim(exp_ref_mat)[2]))
    cov.cell <- c(cell.MCA, new.ref.names)
    mod <- model.matrix(~ as.factor(cov.cell))
    mtx.combat <- ComBat(mtx.in, batch, mod, par.prior = T, ref.batch = 1)
    mtx.combat <- scale(mtx.combat)
    
    # cutoff.fc <- 1.5
    # cutoff.pval <- 0.05
    # topN <- 100
    # cutoff.fc <- NULL
    # cutoff.pval <- NULL
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
        
        mtx.combat.use <- mtx.combat[, cov.cell != new.ref.names[i]]
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


.cluster_sc <- function(exp_sc_mat) {
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
    seurat.unlabeled <- RunPCA(seurat.unlabeled, npcs = 75, verbose = F)
    
    # cluster
    seurat.unlabeled <-
        FindNeighbors(seurat.unlabeled, reduction = "pca", dims = 1:75,
                      nn.eps = 0.5, verbose = F)
    seurat.unlabeled <-
        FindClusters(seurat.unlabeled, resolution = 3, n.start = 20,
                     n.iter = 100, verbose = F)
    
    out.cluster <-
        data.frame(
            cluster.id = as.character(seurat.unlabeled@meta.data$seurat_clusters),
            # original.label = seurat.unlabeled@meta.data$original.label,
            row.names = dimnames(seurat.unlabeled@assays$RNA@counts)[[2]]
        )
    return(out.cluster)
    
}


.cluster_increase_speed <- function(exp_sc_mat, df.cluster, cluster.cell = 5) {
    library(foreach, verbose = F)
    library(doParallel, verbose = F)
    # merge cells in one cluster
    merge.one.cluster <- function(cluster.id) {
        library(Seurat, verbose = F)
        
        cell.ids <- names(df.cluster[df.cluster[, 'cluster.id'] == cluster.id,])
        sub.exp <- exp_sc_mat[, cell.ids]
        # print(dim(sub.exp))
        
        sub.seurat <- CreateSeuratObject(counts = sub.exp)
        sub.seurat <-
            NormalizeData(sub.seurat, normalization.method = "LogNormalize",
                          scale.factor = 10000, verbose = F)
        # print(sub.seurat@assays$RNA@data[1:6,1:6])
        sub.seurat <-
            FindVariableFeatures(sub.seurat, selection.method = "vst",
                                 nfeatures = 1000, verbose = F)
        sub.seurat <- ScaleData(sub.seurat, verbose = F)
        # print(head(sub.seurat@assays$RNA@var.features))
        
        # PCA
        sub.seurat <- RunPCA(sub.seurat, npcs = 10, verbose = F)
        sub.pca <- sub.seurat@reductions$pca@cell.embeddings
        
        num.cell <- dim(sub.pca)[1]
        num.clusters <- floor(num.cell / cluster.cell)
        res.cluster <-
            kmeans(sub.pca, centers = num.clusters, nstart = 20, iter.max = 50)
        
        # names
        sub.dict <- data.frame(
            cell.id = names(res.cluster$cluster),
            cluster.level1 = rep(cluster.id, num.cell),
            cluster.level2 = res.cluster$cluster
        )
        sub.dict$cluster.merge.id <-
            paste(sub.dict$cluster.level1, sub.dict$cluster.level2, sep = '-')
        
        # merge expression profile
        tag.in <- sub.dict[, c('cell.id', 'cluster.merge.id')]
        sub.exp.merge <- .generate_ref(sub.exp, tag.in)
        
        sub.out <- list()
        sub.out$sub.dict <- sub.dict
        sub.out$sub.exp.merge <- sub.exp.merge[sc.genes, ]
        return(sub.out)
        
    }
    
    sc.genes <- row.names(exp_sc_mat)
    df.cluster <- as.matrix(df.cluster)
    cluster.ids <- as.character(unique(df.cluster))
    cluster.cell <- cluster.cell
    num.cpu <- floor(length(cluster.ids) / 10)
    cl <- makeCluster(num.cpu, outfile = '')
    clusterExport(cl, '.generate_ref')
    # registerDoParallel(cl)
    # out.par <- foreach(cluster.id = cluster.ids[1], .inorder = F, .packages = 'Seurat',
    #                    .export = c('.generate_ref')) %dopar%
    #   merge.one.cluster(cluster.id)
    out.par = parLapply(cl = cl, cluster.ids, merge.one.cluster)
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


.confirm_label <- function(exp_sc_mat, list.cell.genes, scRef.tag, CPU=10) {
    # confirm label
    exp_sc_mat <- as.matrix(exp_sc_mat)
    df.tag <- as.data.frame(scRef.tag)
    row.names(df.tag) <- df.tag$cell_id
    df.tag$cell_id <- NULL
    names(df.tag) <- 'scRef.tag'
    
    SINGLE <- function(barcode, exp_sc_mat, df.tag, list.cell.genes) {
        expression.barcode <- exp_sc_mat[, barcode]
        bool.mark.gene <- rep(1, dim(exp_sc_mat)[1])
        cell <- df.tag[barcode, 'scRef.tag']
        genes.marker <- list.cell.genes[[cell]]
        bool.mark.gene[dimnames(exp_sc_mat)[[1]] %in% genes.marker] <- 2
        test.in <- cbind(expression.barcode, bool.mark.gene)
        test.in <- as.data.frame(test.in)
        names(test.in) <- c('expression.level', 'factor.mark.gene')
        test.in$factor.mark.gene <- as.factor(test.in$factor.mark.gene)
        out.test <-
            wilcox.test(
                formula = expression.level ~ factor.mark.gene,
                data = test.in,
                alternative = 'less',
                digits.rank = 7
            )
        pvalue <- out.test$p.value
        gc()
        return(data.frame(pvalue = pvalue, row.names = barcode))
        
    }
    
    cl = makeCluster(CPU, outfile = '')
    RUN <- parLapply(
        cl = cl,
        dimnames(exp_sc_mat)[[2]],
        SINGLE,
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
    
    gc()
    
    return(meta.tag)
    
}


SCREF <- function(exp_sc_mat, exp_ref_mat, type_ref = 'count', out.group = 'MCA', 
                  cluster.speed = F, cluster.cell = 5,
                  method1 = 'kendall', method2 = 'kendall', 
                  cutoff.1 = 'default', cutoff.2 = 'default',
                  min_cell = 20, CPU = 4, print_step = 100) {
    # check parameters
    if (!type_ref %in% c('count', 'fpkm', 'tpm')) {
        stop('Error: inexistent input of reference data format')
    }
    
    time1 <- Sys.time()
    library(doParallel, verbose = F)
    
    # get overlap genes
    out.overlap <- .get_overlap_genes(exp_sc_mat, exp_ref_mat)
    exp_sc_mat <- out.overlap$exp_sc_mat
    exp_ref_mat <- out.overlap$exp_ref_mat
    gene_over <- out.overlap$gene_over
    print('Number of overlapped genes:')
    print(nrow(exp_sc_mat))
    
    # find markers of cell types in reference
    print('Find marker genes of cell types in reference:')
    out.markers <-
        .find_markers_auto(exp_ref_mat, type_ref = type_ref, topN = 100)
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
    df.cluster <- .cluster_sc(exp_sc_mat)
    
    # speed calculation
    if (cluster.speed) {
        print('Speed calculation by clustering:')
        # time1 <- Sys.time()
        out.merge <-
            .cluster_increase_speed(exp_sc_mat, df.cluster, cluster.cell = cluster.cell)
        # time2 <- Sys.time()
        # difftime(time2, time1, units = 'secs')
        df.dict <- out.merge$df.dict
        df.exp.merge <- out.merge$df.exp.merge
        exp_sc_mat <- df.exp.merge
        min_cell <- ceiling(min_cell / cluster.cell)
    }
    
    print('First-round annotation:')
    print(method1)
    if (method1 != 'multinomial') {
        out1 <- .get_cor(exp_sc_mat, exp_ref_mat, method = method1,
                         CPU = CPU, print_step = print_step)
    } else {
        out1 <- .get_log_p_sc_given_ref(exp_sc_mat, exp_ref_mat,
                                        CPU = CPU, print_step = print_step)
    }
    
    gc()
    
    tag1 <- .get_tag_max(out1)
    # df.tag.1 <- as.data.frame(tag1)
    # row.names(df.tag.1) <- df.tag.1$cell_id
    # df.tag.1$cell_id <- NULL
    # names(df.tag.1) <- 'scRef.tag'
    # df.view.1 <- merge(df.tag.1, label.filter, by = 'row.names')
    
    print('Build local reference')
    exp_sc_mat <- as.matrix(exp_sc_mat)
    df.tags1 <- .confirm_label(exp_sc_mat, list.cell.genes, tag1, CPU = CPU)
    gc()
    df.tags1 <- df.tags1[dimnames(exp_sc_mat)[[2]], ]
    df.tags1$log10Pval <- -log10(df.tags1$pvalue)
    # df.tags1.view <- merge(df.tags1, label.filter, by = 'row.names')
    # select cutoff.1 automatitically
    if (cutoff.1 == 'default') {
        tag.base <- data.frame(cell_id = names(list.cell.genes),
                               tag = names(list.cell.genes))
        df.base <- .confirm_label(exp_ref_mat, list.cell.genes, tag.base, CPU = 2)
        cutoff.1 <- -(log10(median(df.base$pvalue)) * 0.5)
        print('Default cutoff: ')
        print(cutoff.1)
    }
    # cutoff.1 <- 1e-20
    select.df.tags <- df.tags1[df.tags1$log10Pval > cutoff.1, ]
    select.barcode <- row.names(select.df.tags)
    LocalRef = .generate_ref(exp_sc_mat[, select.barcode],
                             tag1[tag1[, 'cell_id'] %in% select.barcode, ], 
                             min_cell = min_cell)
    print('Cell types in local reference:')
    print(dimnames(LocalRef)[[2]])
    
    #####
    gc()
    #####
    # find local marker genes
    print('find local marker genes')
    out.markers <- .find_markers_auto(LocalRef, topN = 100)
    local.cell.genes <- out.markers[['list.cell.genes']]
    
    print('Second-round annotation:')
    print(method2)
    if (method2 != 'multinomial') {
        out2 <- .get_cor(exp_sc_mat, LocalRef, method = method2,
                         CPU = CPU, print_step = print_step)
    } else {
        out2 <- .get_log_p_sc_given_ref(exp_sc_mat, LocalRef,
                                       CPU = CPU, print_step = print_step)
    }
    tag2 <- .get_tag_max(out2)
    gc()
    
    df.tags2 <- .confirm_label(exp_sc_mat, local.cell.genes, tag2, CPU = CPU)
    gc()
    df.tags2 <- df.tags2[dimnames(exp_sc_mat)[[2]], ]
    # df.tags2.view <- merge(df.tags2, label.filter, by = 'row.names')
    # row.names(df.tags2.view) <- df.tags2.view$Row.names
    # df.tags2.view$Row.names <- NULL
    # names(df.tags2.view) <- c('scRef2.tag', 'pvalue', 'qvalue', 'label')
    
    print('Combine reference and local result:')
    # concat reference pval and local pval
    pvalue1 <- df.tags1[, c('scRef.tag', 'pvalue')]
    names(pvalue1) <- c('scRef.tag.1', 'pvalue.1')
    pvalue2 <- df.tags2
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
            scRef.tag.12 = tag.final,
            pvalue = pval.final,
            row.names = dimnames(mtx.tag)[[1]]
        )
    
    # modify tags and combine pval
    df.tags <- merge(pvalue, tag.final, by = 'row.names')
    row.names(df.tags) <- df.tags$Row.names
    df.tags$Row.names <- NULL
    
    df.tags$log10Pval <- -log10(df.tags$pvalue)
    
    if (cluster.speed) {
        df.cluster <- df.dict[, c("cluster.merge.id", "cluster.level1")]
        df.cluster <- unique(df.cluster)
        df.cluster <- data.frame(cluster.id = df.cluster$cluster.level1,
                                 row.names = df.cluster$cluster.merge.id)
    }
    
    # select cutoff.2 automatitically
    if (cutoff.2 == 'default') {
        df.tag.cluster <- merge(df.tags, df.cluster, by = 'row.names')
        cluster.ids <- unique(df.cluster$cluster.id)
        meta.cluster <- data.frame(stringsAsFactors = F)
        for (cluster.id in cluster.ids) {
            sub.tag.cluster <- df.tag.cluster[df.tag.cluster$cluster.id == cluster.id,]
            mean.log10Pval <- mean(sub.tag.cluster$log10Pval)
            # sd.log10Pval <- sd(sub.tag.cluster$log10Pval)
            table.tag <- table(sub.tag.cluster$scRef.tag.12)
            table.tag <- table.tag[order(table.tag, decreasing = T)]
            main.cell <- names(table.tag[1])
            percent.main.cell <- table.tag[1] / dim(sub.tag.cluster)[1]
            min.log10Pval <-
                min(sub.tag.cluster[sub.tag.cluster$scRef.tag.12 == main.cell, 'log10Pval'])
            meta.cluster <- rbind(
                meta.cluster,
                data.frame(
                    cluster.id = cluster.id,
                    mean.log10Pval = mean.log10Pval,
                    min.log10Pval = min.log10Pval,
                    main.cell = main.cell,
                    percent.main.cell = percent.main.cell
                )
            )
        }
        meta.cluster <- meta.cluster[order(meta.cluster$mean.log10Pval, decreasing = T), ]
        cluster.ids <- meta.cluster$cluster.id
        vec.min <- c()
        vec.cell <- c()
        for (cluster.id in cluster.ids) {
            one.cell <- meta.cluster[meta.cluster$cluster.id == cluster.id, 'main.cell']
            if (!(one.cell %in% vec.cell)) {
                one.min <- meta.cluster[meta.cluster$cluster.id == cluster.id, 'min.log10Pval']
                vec.cell <- c(vec.cell, one.cell)
                vec.min <- c(vec.min, one.min)
            }
        }
    }
    cutoff.2 <- median(vec.min)
    print('Default cutoff: ')
    print(cutoff.2)
    df.tags$scRef.tag <- df.tags$scRef.tag.12
    df.tags$scRef.tag[df.tags$log10Pval < cutoff.2] <- 'Unassigned'
    
    if (cluster.speed) {
        df.tags$cluster.merge.id <- row.names(df.tags)
        df.tags.merge <- merge(df.tags, df.dict[, c('cluster.merge.id', 'cell.id')],
                               by = 'cluster.merge.id')
        df.tags.merge$cluster.merge.id <- NULL
        row.names(df.tags.merge) <- df.tags.merge$cell.id
        df.tags.merge$cell.id <- NULL
        df.tags <- df.tags.merge
    }
    
    
    df.combine <- df.tags[, c("scRef.tag", "log10Pval")]
    
    # df.view <- merge(label.filter, combine.out, by = 'row.names')
    
    # df.view <- merge(df.combine, df.tags, by = 'row.names')
    # row.names(df.view) <- df.view$Row.names
    # df.view$Row.names <- NULL
    # df.view <- merge(label.filter, df.view, by = 'row.names')
    
    #####
    gc()
    #####
    time2 <- Sys.time()
    time.scRef <- difftime(time2, time1, units = 'secs')
    output = list()
    output$tag1 = tag1
    output$out1 = out1
    output$pvalue1 = df.tags1
    output$tag2 = tag2
    output$out2 = out2
    output$pvalue2 = df.tags2
    output$final.out <- df.combine
    output$combine.out <- df.tags
    output$run.time <- time.scRef
    if (cluster.speed) {
        output$dict.cluster <- df.dict
    }
    
    print('Finish!')
    
    return(output)
}
