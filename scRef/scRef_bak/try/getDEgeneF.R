# esetm: gene expression matrix
# group: factor: "case"/"control"
# pair: TRUE or FALSE*
# block: e.g.1 2 2 1 if paired; blank if not
# p_adj: p.adjust (fdr*)
# fpkm: TRUE (RNAseq FPKM or RSEM value) or FALSE* (chipseq)
getDEgeneF<-function(esetm=NULL,group=NULL,pair=FALSE,block=NULL,p_adj="fdr",fpkm=T){
        if (is.null(esetm)) {
                cat("esetm: gene expression matrix",
                           "group: factor: \"c\"/\"d\"",
                           "pair: TRUE/FALSE*",
                           "block: e.g.1 2 2 1 if paired; blank if not",
                           "p_adj: p.adjust, fdr* ",
                           "fpkm: TRUE/FALSE*",        
                           sep="\n")
        }else{
                library(limma)
                if(pair){
                        design<-model.matrix(~block+group)
                }else{
                        design<-model.matrix(~group)
                }
                fit<-lmFit(esetm,design)
                if(fpkm){
                        fit<-eBayes(fit,trend=T,robust=T)
                }else{
                        fit<-eBayes(fit)
                }
                x<-topTable(fit,number=nrow(esetm),adjust.method=p_adj,coef="group1")
                x<-x[!is.na(row.names(x)),]
                x<-x[!duplicated(row.names(x)),]
                return(x)    
        }
}

