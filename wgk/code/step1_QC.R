###################
##  QC: 每个样本单独质控
###################
library(Seurat)
library(dplyr)
library(ggplot2)

batchID <- c(paste('C',1:6,sep=''),paste('M',1:3,sep=''))
Gene.min <- c(1400,1500,1500,1100,1300,1500,800,800,900)
Gene.max <- c(6000,6000,6000,6000,6000,5000,6000,6000,6200)
UMI.min <- c(4200,4000,5200,6200,6000,4500,2000,3000,3000)
UMI.max <- c(50000,50000,50000,50000,50000,25000,50000,50000,50000)
percent.mt <- c(15,15,15,10,15,15,15,15,10)

names(Gene.min)=names(Gene.max)=names(UMI.min)=names(UMI.max)=names(percent.mt)=batchID

ID='C1'
DATA <- Read10X(paste('/home/yzj/JingMA_ORI/data/',ID,'/',sep=''))
pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = ID)
pbmc[["percent.mt"]] <- PercentageFeatureSet(pbmc, pattern = "^MT-")
pbmc@meta.data$batch <- rep(ID,ncol(pbmc))

###############
# nfeature distribution
h_df <- data.frame(pbmc$batch, pbmc$nFeature_RNA)
colnames(h_df) <- c("Sample", "nFeature_RNA")
p.Gene <- ggplot(h_df, aes(x=nFeature_RNA, color=Sample)) + geom_histogram(binwidth = 100, fill = "white") +
  theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Cell Number")+
  geom_vline(aes(xintercept=Gene.min[ID]), colour="black", linetype="dashed")+
  geom_vline(aes(xintercept=Gene.max[ID]), colour="black", linetype="dashed")+
  annotate('text',x=Gene.min[ID],y=200,label=as.character(Gene.min[ID]),parse=T)+
  annotate("text",x=Gene.max[ID],y=200,label=as.character(Gene.max[ID]),parse=T)
###############

###############
# ncount distribution
h2_df <- data.frame(pbmc$batch, pbmc$nCount_RNA)
colnames(h2_df) <- c("Sample", "nCount_RNA")
p.UMI <- ggplot(h2_df, aes(x=nCount_RNA, color=Sample)) + geom_histogram(binwidth = 100, fill="white") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Cell Number") +
  geom_vline(aes(xintercept=UMI.min[ID]), colour="black", linetype="dashed")+
  geom_vline(aes(xintercept=UMI.max[ID]), colour="black", linetype="dashed")+
  annotate('text',x=UMI.min[ID],y=40,label=as.character(UMI.min[ID]),parse=T)+
  annotate("text",x=UMI.max[ID],y=40,label=as.character(UMI.max[ID]),parse=T)
###############

###############
# mt distribution
h3_df <- data.frame(pbmc$batch, pbmc$percent.mt)
colnames(h3_df) <- c("Sample", "percent.mt")
p.mt <- ggplot(h3_df, aes(x=percent.mt, color=Sample)) + geom_histogram(binwidth = 2, fill="white") + theme(panel.background = element_blank(), axis.line = element_line(colour = "black")) + ylab("Cell Number") +
  geom_vline(aes(xintercept=percent.mt[ID]), colour="black", linetype="dashed")+
  annotate("text",x=percent.mt[ID],y=2000,label=as.character(percent.mt[ID]),parse=T)
###############

############
subpbmc <- subset(pbmc, subset = nFeature_RNA > Gene.min[ID] & nFeature_RNA < Gene.max[ID] &
                    nCount_RNA > UMI.min[ID] & nCount_RNA < UMI.max[ID] & 
                    percent.mt < percent.mt[ID])
############

############
saveRDS(subpbmc,paste('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_',ID,'.RDS',sep=''))
############



############
### merge
############
pbmc_C1 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_C1.RDS')
pbmc_C2 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_C2.RDS')
pbmc_C3 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_C3.RDS')
pbmc_C4 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_C4.RDS')
pbmc_C5 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_C5.RDS')
pbmc_C6 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_C6.RDS')

pbmc_M1 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_M1.RDS')
pbmc_M2 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_M2.RDS')
pbmc_M3 <- readRDS('/home/yzj/JingMA_NEW/res/QC/ALL/RDS/pbmc_QC_M3.RDS')

pbmc.combined <- merge(pbmc_C1, y = list(pbmc_C2,pbmc_C3,pbmc_C4,pbmc_C5,pbmc_C6,pbmc_M1,pbmc_M2,pbmc_M3),
                       add.cell.ids = c("C1", "C2", "C3", "C4","C5" ,"C6","M1", "M2", "M3"), project = "Microtia")
print(table(pbmc.combined$orig.ident))
pbmc.combined@meta.data$batch <- c(rep('C1',ncol(pbmc_C1)),rep('C2',ncol(pbmc_C2)),rep('C3',ncol(pbmc_C3)),
                                   rep('C4',ncol(pbmc_C4)),rep('C5',ncol(pbmc_C5)),rep('C6',ncol(pbmc_C6)),
                                   rep('M1',ncol(pbmc_M1)),rep('M2',ncol(pbmc_M2)),rep('M3',ncol(pbmc_M3)))
pbmc.combined@meta.data$type <- c(rep('Control',ncol(pbmc_C1)+ncol(pbmc_C2)+ncol(pbmc_C3)+ncol(pbmc_C4)+ncol(pbmc_C5)+ncol(pbmc_C6)),
                                  rep('Microtia',ncol(pbmc_M1)+ncol(pbmc_M2)+ncol(pbmc_M3)))

system("mkdir -p JingMA_NEW/res/QC/ALL/RDS/")
saveRDS(pbmc.combined,'JingMA_NEW/res/QC/ALL/RDS/PBMC_QC.RDS')



############
## QC info of sample
############
library(Seurat)
library(dplyr)
# before
before.info <- c()
sampleID <- c(paste('C',1:6,sep=''),paste('M',1:3,sep=''))
for(i in 1:length(sampleID)){
  ID=sampleID[i]
  DATA <- Read10X(paste('/home/yzj/JingMA_ORI/data/',ID,'/',sep=''))
  pbmc=CreateSeuratObject(counts = DATA, min.cells = 0, min.features = 0, project = ID)
  info <- c(ncol(pbmc),median(pbmc$nFeature_RNA),median(pbmc$nCount_RNA))
  before.info <- rbind(info,before.info)
}
before.info <- as.data.frame(before.info)
colnames(before.info) <- c('cellNumber','median_nFeature','median_nCount')
rownames(before.info) <- sampleID

# after
pbmc <- readRDS('/home/yzj/JingMA_NEW/res/Harmony/ALL/RDS/seurat_celltype.Rdata')
meta.df <- pbmc@meta.data
meta.df <- group_by(meta.df, batch)

after.info <- summarise(meta.df,
                        cellNumber = n(), 
                        median_nFeature = median(nFeature_RNA),    #平均值
                        median_nCount = median(nCount_RNA)) 

QC.info <- list(before=before.info,afer=after.info)
saveRDS(QC.info,'/home/yzj/JingMA_NEW/res/QC/ALL/RDS/QC_info.RDS')



