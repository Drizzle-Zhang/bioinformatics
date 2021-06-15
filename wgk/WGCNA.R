setwd('/home/zy/my_git/bioinformatics/wgk')
.libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(dplyr)
library(WGCNA)
library(ggplot2)
options(stringsAsFactors = FALSE);
#out_wp='/home/yzj/JingMA/res/WGCNA/Control/'
out_wp='/home/disk/drizzle/wgk/microtia_chon_child_M1M2M3/WGCNA/'
if (!file.exists(out_wp)) {
  dir.create(out_wp)
}
setwd(out_wp)
# enableWGCNAThreads(nThreads=6)

#=====================================================================================
#1.Data_input_and_cleaning
#=====================================================================================
#EXP <- readRDS('/home/yzj/JingMA/res/VECTOR/ALL/Control_EXP_PS.RDS')
# EXP <- readRDS('/home/yzj/JingMA_ORI/res/VECTOR/ALL/Microtia_EXP_PS.RDS')
path.data <- '/home/disk/drizzle/wgk/data/AllSample_2_merge/'
path.lineage <- paste0(path.data, 'chon_lineage/')
file.chon <- paste0(path.lineage, 'seurat_celltype.Rdata')
seurat.chon <- readRDS(file.chon)
seurat.child <- subset(seurat.chon, subset = batch %in% c('C4', 'C6', 'M1', 'M2', 'M3'))
# seurat.child <- NormalizeData(seurat.child)
# seurat.child <- FindVariableFeatures(seurat.child, nfeatures = 3000)
highvar.genes <- VariableFeatures(seurat.child)

EXP <- seurat.child@assays$RNA@data
# sample.cells <- sample(colnames(EXP), 3000)
sample.cells <- colnames(EXP)
EXP <- EXP[, sample.cells]
dim(EXP)

# m.mad <- apply(EXP,1,mad)
# dataExprVar <- EXP[which(m.mad >
#                                 max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
# dataExprVar <- EXP[which(m.mad > 0.01),]
# dim(dataExprVar)
dataExprVar <- EXP[highvar.genes,]
dim(dataExprVar)

dataExpr <- as.data.frame(t(as.matrix(dataExprVar)))
#check the missing values
# gsg = goodSamplesGenes(dataExpr, verbose = 3);
# if (!gsg$allOK){
#   if (sum(!gsg$goodGenes)>0) 
#     printFlush(paste("Removing genes:", paste(names(datExpr0)[!gsg$goodGenes], collapse = ", ")));
#   if (sum(!gsg$goodSamples)>0) 
#     printFlush(paste("Removing samples:", paste(rownames(datExpr0)[!gsg$goodSamples], collapse = ", ")));
#   datExpr0 = dataExpr[gsg$goodSamples, gsg$goodGenes]
# }
datExpr0 <- dataExpr

#*** clust the samples and check the outlier . You will got fig.1 ***#
#------------------------------------------------------------------------------------#
sampleTree = hclust(dist(datExpr0), method = "average");
file.sampleTree <- "./sampleTree.Rdata"
save(sampleTree, file = file.sampleTree)
sampleTree <- readRDS(file.sampleTree)

pdf(file = paste(out_wp,"Plots/fig1.pdf",sep=''), width = 10, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5, cex.axis = 1.5, cex.main = 2)
# Plot a line to show the cut
abline(h = 250, col = "red");
dev.off()
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 250, minSize = 10)
table(clust)
keepSamples = (clust==1)
datExpr = datExpr0[keepSamples, ]
#datExpr = datExpr0 # not exclude any samples
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

#***  set your traits matrix , import the traits information ***#
#------------------------------------------------------------------------------------#
datTraits <- seurat.child@meta.data[sample.cells, c("batch", "type", "celltype")]
# datTraits <- data.frame(PS=as.numeric(rownames(datExpr)),Type=0,row.names = rownames(datExpr))
head(datTraits)
save(datTraits, file = "datTraits.RData")
collectGarbage();

#*** to do clust analisis. you wil get fig.2 ***#
#------------------------------------------------------------------------------------#
# Re-cluster samples
# sampleTree2 = hclust(dist(datExpr), method = "average")
# traitColors = numbers2colors(datTraits, signed = FALSE);
# 
# 
# # Plot the sample dendrogram and the colors underneath.that means put them together
# pdf(paste(out_wp,"Plots/fig2.pdf",sep=''), width = 10, height = 9)
# plotDendroAndColors(sampleTree2, traitColors,
#                     groupLabels = names(datTraits), 
#                     main = "Sample dendrogram and trait heatmap")
# dev.off()
save(datExpr, datTraits, file = "OUT-01-dataInput.RData")


#=====================================================================================
#  2. Network construction and module detection(here use one-step, and also have step by step and block-wise network)
#=====================================================================================
allowWGCNAThreads(nThreads=6)
# lnames = load(file = "OUT-01-dataInput.RData");
load(file = "OUT-01-dataInput.RData")
nGenes = ncol(datExpr)
nSamples = nrow(datExpr)

# Choose a set of soft-thresholding powers ***#
#------------------------------------------------------------------------------------#
powers = c(c(1:10), seq(from = 12, to=20, by=2))
# Call the network topology analysis function
sft = pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)
power = sft$powerEstimate
print(power)

pdf(paste(out_wp,'Plots/softThres.pdf',sep=''),12,9)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red");
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

#***  One-step network construction and module detection  ***#
#------------------------------------------------------------------------------------#
net = blockwiseModules(datExpr, power = sft$powerEstimate,
                       TOMType = "unsigned", minModuleSize = 20,
                       reassignThreshold = 0, mergeCutHeight = 0.15,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "OUT_TOM", 
                       verbose = 3)

#***  Module visualization = > cluster tree  ***#
#------------------------------------------------------------------------------------#
pdf(paste(out_wp,'Plots/fig3.pdf',sep=''),12,9)
mergedColors = labels2colors(net$colors)

plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)
dev.off()

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(net,MEs, moduleLabels, moduleColors, geneTree, file = "OUT-02-networkConstruction-auto.RData")

load("OUT-02-networkConstruction-auto.RData")
module = "turquoise";
module = "blue";
module = "brown";
module = "yellow";
# Select module probes
probes = colnames(datExpr) ## 我们例子里面的probe就是基因
inModule = (moduleColors==module);
modProbes = probes[inModule]; 
sort(modProbes)

#=====================================================================================
#3.Relating modules2external clinical traits
#=====================================================================================
# Recalculate MEs with color labels
MEs0 = moduleEigengenes(datExpr, moduleColors)$eigengenes
MEs = orderMEs(MEs0)

#*** Caculate the correlation of MEs(Module Eigengenes) and traits.  moduleTraitCor ***#
#------------------------------------------------------------------------------------#

datTraits.1 <- apply(datTraits[, c('type', 'celltype')], 1, 
                     function(x) {paste(x[1], x[2], sep = '_')})
moduleTraitCor = cor(MEs, datTraits.1, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples);

# Plot the heatmap for moduleTraitCor.  You will get fig. 4  ***#
#------------------------------------------------------------------------------------#
pdf(paste(out_wp,'Plots/fig4.pdf',sep=''),8,16)
textMatrix =  paste(signif(moduleTraitCor, 2), "\n(",
                    signif(moduleTraitPvalue, 1), ")", sep = "");
dim(textMatrix) = dim(moduleTraitCor)
par(mar = c(6, 8.5, 3, 3));
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = names(datTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = greenWhiteRed(50),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text = 1,
               zlim = c(-1,1),
               main = paste("Module-trait relationships"))
dev.off()


#***   Caculate the correlation of genes and Module (MM: Module Membership)  ***#
#------------------------------------------------------------------------------------#
# Define variable PS containing the PS column of datTrait
PS = as.data.frame(datTraits$PS);
names(PS) = "PS"
modNames = substring(names(MEs), 3)

geneModuleMembership = as.data.frame(cor(datExpr, MEs, use = "p"));
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples));

names(geneModuleMembership) = paste("MM", modNames, sep="");
names(MMPvalue) = paste("p.MM", modNames, sep="");

#***   Caculate the correlation of genes and traits (GS: Gene Significance)  ***#
#------------------------------------------------------------------------------------#
geneTraitSignificance = as.data.frame(cor(datExpr, PS, use = "p"));
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples));

names(geneTraitSignificance) = paste("GS.", names(PS), sep="");
names(GSPvalue) = paste("p.GS.", names(PS), sep="");


#***   Get geneInfo***#
#------------------------------------------------------------------------------------#
# Create the starting data frame
geneInfo0 = data.frame(ID = colnames(datExpr),
                       moduleColor = moduleColors,
                       geneTraitSignificance,
                       GSPvalue)
# Order modules by their significance for PS
modOrder = order(-abs(cor(MEs, PS, use = "p")));
# Add module membership information in the chosen order
for (mod in 1:ncol(geneModuleMembership)){
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[, modOrder[mod]], 
                         MMPvalue[, modOrder[mod]]);
  names(geneInfo0) = c(oldNames, paste("MM.", modNames[modOrder[mod]], sep=""),
                       paste("p.MM.", modNames[modOrder[mod]], sep=""))
}
# Order the genes in the geneInfo variable first by module color, then by geneTraitSignificance
geneOrder = order(geneInfo0$moduleColor, -abs(geneInfo0$GS.PS));
geneInfo = geneInfo0[geneOrder, ]
write.csv(geneInfo, file = "geneInfo.csv")


# Pick hubgenes which cor > 0.7 with p <0.05 in grey module
geneInfo <- read.csv('geneInfo.csv')
modules <- unique(moduleColors)
for(i in 1:length(modules)){
  targetModule <- modules[i]
  print(targetModule)
  if(targetModule != 'grey'){
    geneInfo_target <- geneInfo[which(geneInfo$moduleColor==targetModule),
                                c(1:4,
                                  which(colnames(geneInfo)==paste('MM.',targetModule,sep='')),
                                  which(colnames(geneInfo)==paste('p.MM.',targetModule,sep='')))]
    
    hubgeneInfo_target <- geneInfo_target[which(geneInfo_target[,ncol(geneInfo_target)-1] > 0.7 & 
                                                  geneInfo_target[,ncol(geneInfo_target)] < 0.05),]
    
    print(dim(hubgeneInfo_target))
    write.table(hubgeneInfo_target,file=paste('hubgene_info_',targetModule,'.txt',sep=''),row.names=F,col.names=T,sep='\t',quote=F)
  }
}


#=====================================================================================
#  4. Module pattern along PS
#=====================================================================================
rm(list=ls())
out_wp='/home/yzj/JingMA/res/WGCNA/Control/'
#out_wp='/home/yzj/JingMA/res/WGCNA/Microtia/'
setwd(out_wp)
lnames = load(file = "OUT-01-dataInput.RData");
geneInfo<- read.csv('geneInfo.csv')
dim(geneInfo)
modules <- unique(geneInfo$moduleColor)

EXP_mean <- c()
for(i in 1:length(modules)){
  targetModule <- modules[i]
  print(targetModule)
  gene_target <- geneInfo$ID[which(geneInfo$moduleColor==targetModule)]
  exp <- dplyr::select(as.data.frame(datExpr),gene_target)
  exp_mean <- apply(exp,1,mean)
  EXP_mean <-cbind(EXP_mean,exp_mean)
}
colnames(EXP_mean) <- modules

max_ps <- c()
for(i in 1:ncol(EXP_mean)){
  index<- which(EXP_mean[,i]==max(EXP_mean[,i]))
  max_ps <- c(max_ps,index)
}
names(max_ps) <- colnames(EXP_mean)
modules_sort <- names(sort(max_ps))

moduleID <- paste('M',1:length(modules_sort),sep='')
names(moduleID) <- modules_sort
save(max_ps,moduleID, file = "moduleID.RData")

geneInfo$moduleID <- moduleID[geneInfo$moduleColor]
dim(geneInfo)
write.csv(geneInfo, file = "geneInfo.csv",row.names = F)


pdf('Plots/Modules_PS.pdf',width = 20,height = 20)
RCN=trunc(sqrt(length(moduleID))+1)
par(mfrow=c(RCN,RCN))
for(i in 1:length(moduleID)){
  targetModuleID <- moduleID[i]
  targetModule <- names(targetModuleID)
  print(targetModule)
  if(targetModule != 'grey'){
    gene_target <- geneInfo$ID[which(geneInfo$moduleColor==targetModule)]
    exp <- dplyr::select(as.data.frame(datExpr),gene_target)
    exp_mean <- apply(exp,1,mean)
    plot(exp_mean,type='l',col='blue',main=targetModuleID,col.main=targetModule,
         lwd=2,cex.axis=1.5,cex.main=2.5,xlab = 'PS',ylab = 'Mean EXP')
  }#col.main=targetModule
}
dev.off()


#######################
GN <- geneInfo[,c(2:3,ncol(geneInfo))]
GN$moduleID <- factor(GN$moduleID,levels = moduleID)
ggplot() + geom_bar(data = GN, aes(x = moduleID), stat = "count",fill= names(moduleID))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = 12))+
  theme(axis.text.y = element_text(face = 'bold',size = 12))

rm(list=ls()) 



#######################
## Correlation heatmap of para samples
# Get upper triangle of the correlation matrix
library(reshape2)
get_upper_tri <- function(cormat){
  cormat[lower.tri(cormat)]<- NA
  return(cormat)
}

reorder_cormat <- function(cormat){
  # Use correlation between variables as distance
  dd <- as.dist((1-cormat)/2)
  hc <- hclust(dd)
  cormat <-cormat[hc$order, hc$order]
  return(cormat)
}


df <- as.data.frame(EXP_mean,stringsAsFactors = F)
df <- dplyr::select(df,names(moduleID))
colnames(df) <- moduleID
cormat <- round(cor(df),2)
save(cormat,file = '/home/yzj/JingMA/res/WGCNA/Control/cormat_M2M.RData')

#cormat <- reorder_cormat(cormat)
upper_tri <- get_upper_tri(cormat)
# Melt the correlation matrix
melted_cormat <- melt(upper_tri, na.rm = TRUE)

pdf('/home/yzj/JingMA/res/WGCNA/Control/Plots/COR.pdf')
ggheatmap <- ggplot(melted_cormat, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Pearson\nCorrelation") +
  theme_minimal()+ # minimal theme
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 9, hjust = 1))+
  coord_fixed()
ggheatmap + 
  geom_text(aes(Var2, Var1, label = value), color = "black", size = 2) +
  theme(axis.title.x = element_blank(),axis.title.y = element_blank(),
        panel.grid.major = element_blank(),panel.border = element_blank(),
        panel.background = element_blank(),axis.ticks = element_blank(),
        legend.justification = c(1, 0),legend.position = c(0.6, 0.7),legend.direction = "horizontal")+
  guides(fill = guide_colorbar(barwidth = 7, barheight = 1,title.position = "top", title.hjust = 0.5))
dev.off()




#=====================================================================================
#  5. Levels of modules acrrodiing TF-Genes relationship
#=====================================================================================
geneInfo<- read.csv('/home/yzj/JingMA/res/WGCNA/Control/geneInfo.csv')
lnames = load(file = "/home/yzj/JingMA/res/WGCNA/Control/moduleID.RData")
lnames = load(file = "/home/yzj/JingMA/res/WGCNA/Control/OUT-01-dataInput.RData")
lnames = load(file = '/home/yzj/JingMA/res/WGCNA/Control/cormat_M2M.RData')

####
lnames = load(file = "/home/yzj/JingMA/res/VECTOR/ALL/SIG_geneTF.RData");
lnames = load(file='/home/yzj/publicData/TRRUST/trrust_rawdata.human_Known.RData')

####
G2M_map <- geneInfo$moduleID
names(G2M_map) <- geneInfo$ID
####

sigGene <- intersect(geneInfo$ID,SIG_gene$Gene)
length(sigGene)

sigTF <- intersect(geneInfo$ID,SIG_TF$Gene)
length(sigTF)

## 0. BG: 所有Gene之间的调控关系(bgTF2TF)  >> bgGene所在的模块之间的调控关系(M2M)
bgGene <- colnames(datExpr)
bgG2G <- TF_Gene[TF_Gene$Gene %in% bgGene & TF_Gene$TF %in% bgGene,]
bgG2G <- bgG2G[-which(bgG2G$TF==bgG2G$Gene),]
dim(bgG2G)

# 0.1 M2M and M2M_SELECT
M2M_all <- bgG2G
M2M_all$TF <- G2M_map[M2M_all$TF]
M2M_all$Gene <- G2M_map[M2M_all$Gene]

M2M_Freq <- as.data.frame(table(M2M_all$TF,M2M_all$Gene,M2M_all$RelationShip))
colnames(M2M_Freq) <- c('TF','Gene','RelationShip','Freq')
M2M_nonDup <- M2M_all[!duplicated(M2M_all[,c(1,2)]),]
M2M <- merge(M2M_Freq,M2M_nonDup,by=c('TF','Gene','RelationShip'))
dim(M2M)

Combn <- combn(union(M2M$TF,M2M$Gene), 2, simplify=F)

M2M_pick <- c()
# 四选一
for(i in 1:length(Combn)){
  combn <- Combn[[i]]
  Freq <- M2M[(M2M$TF==combn[1] & M2M$Gene==combn[2]) | (M2M$TF==combn[2] & M2M$Gene==combn[1]),]
  if(nrow(Freq) > 0){
    if(nrow(Freq) >1 & length(unique(Freq$Freq)) >1){
      index <- which(Freq$Freq==max(Freq$Freq))
      M2M_pick <- rbind(M2M_pick,Freq[index,])
    }else{
      M2M_pick <- rbind(M2M_pick,Freq)
    }
  }
}

# ## 四选二
# for(i in 1:length(Combn)){
#   combn <- Combn[[i]]
#   Freq_F <- M2M[(M2M$TF==combn[1] & M2M$Gene==combn[2]),]
#   Freq_A <- M2M[(M2M$TF==combn[2] & M2M$Gene==combn[1]),]
#   if(nrow(Freq_F) > 0){
#     if(nrow(Freq_F) >1 & length(unique(Freq_F$Freq)) >1){
#       index <- which(Freq_F$Freq==max(Freq_F$Freq))
#       M2M_pick <- rbind(M2M_pick,Freq_F[index,])
#     }else{
#       M2M_pick <- rbind(M2M_pick,Freq_F)
#     }
#   }
#   if(nrow(Freq_A) > 0){
#     if(nrow(Freq_A) >1 & length(unique(Freq_A$Freq)) >1){
#       index <- which(Freq_A$Freq==max(Freq_A$Freq))
#       M2M_pick <- rbind(M2M_pick,Freq_A[index,])
#     }else{
#       M2M_pick <- rbind(M2M_pick,Freq_A)
#     }
#   }
# }


dim(M2M_pick)
tmp <- M2M_pick[,c(1,3,2)]
write.table(tmp,'/home/yzj/JingMA/res/WGCNA/Control/M2M_pick.txt',row.names = FALSE,col.names = TRUE,sep='\t',quote = F)



## 进一步CUT network
PS_Time <- c(rep('Up',5),rep('M',5),rep('Down',5))
names(PS_Time) <- 1:15

M2M_SELECT <- c()
for(i in 1:nrow(M2M_pick)){
  source_m <- M2M_pick$TF[i]
  source_col <- names(moduleID)[which(moduleID==source_m)] 
  source_maxPS <- max_ps[source_col]
  source_PSTime <- PS_Time[source_maxPS]
  
  target_m <- M2M_pick$Gene[i]
  target_col <- names(moduleID)[which(moduleID==target_m)] 
  target_maxPS <- max_ps[target_col]
  target_PSTime <- PS_Time[target_maxPS]
  
  RS <- M2M_pick$RelationShip[i]
  
  C1 <- source_PSTime=='Up' & target_PSTime=='M' & RS=='Activation'
  C2 <- source_PSTime=='Up' & target_PSTime=='Down' & RS=='Repression'
  C3 <- source_PSTime=='M' & target_PSTime=='Down' & RS=='Activation'
  C4 <- source_PSTime=='M' & target_PSTime=='Up' & RS=='Notsure'
  C5 <- source_PSTime=='Down' & target_PSTime=='M' & RS=='Notsure'
  C6 <- source_PSTime=='Down' & target_PSTime=='Up' & RS=='Repression'
  C7 <- ( (source_PSTime=='Up' & target_PSTime=='Up') | (source_PSTime=='M' & target_PSTime=='M') |
            (source_PSTime=='Down' & target_PSTime=='Down') ) &  
    (as.integer(names(source_PSTime)) <= as.integer(names(target_PSTime))) & RS=='Activation'
  C8 <- ( (source_PSTime=='Up' & target_PSTime=='Up') | (source_PSTime=='M' & target_PSTime=='M') |
            (source_PSTime=='Down' & target_PSTime=='Down') ) &  
    (as.integer(names(source_PSTime)) > as.integer(names(target_PSTime))) & RS=='Notsure'
  
  if(C1|C2|C3|C4|C5|C6|C7|C8){
    M2M_SELECT <- rbind(M2M_SELECT,M2M_pick[i,])
  }
}

tmp <- M2M_SELECT[,c(1,3,2)]
dim(tmp)
write.table(tmp,'/home/yzj/JingMA/res/WGCNA/Control/M2M_SELECT.txt',row.names = FALSE,col.names = TRUE,sep='\t',quote = F)


## CEBPD
bgG2G[bgG2G$Gene == 'CEBPD' & bgG2G$RelationShip == 'Activation',]


# 0.2 dotplot: in-out degree plot
OUT=as.data.frame(table(M2M_SELECT$TF))
IN=as.data.frame(table(M2M_SELECT$Gene))
D <- merge(IN,OUT,by='Var1',all = T)
D[is.na(D)] <- 0
colnames(D) <- c('moduleID','IN','OUT')
FillCol <- names(moduleID)[which(moduleID %in% D$moduleID)]

ggplot(D, aes(x=IN, y=OUT, colour=FillCol)) + 
  geom_point() + 
  geom_text(aes(label=moduleID), size=5) + 
  scale_fill_discrete(guide=FALSE)+
  geom_abline(a=1, b=1, col="red")+
  guides(fill=FALSE)+
  theme_bw(base_size = 15)



## 0.3 barplot: bgGene distribution in module
GN <- geneInfo[geneInfo$ID %in% bgGene,c(2:3,ncol(geneInfo))]
GN$moduleID <- factor(GN$moduleID,levels = moduleID[moduleID %in% unique(GN$moduleID)])

ggplot() + geom_bar(data = GN, aes(x = moduleID), stat = "count",
                    fill= names(moduleID[moduleID %in% unique(GN$moduleID)]))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = 12))+
  theme(axis.text.y = element_text(face = 'bold',size = 12))




## 1. 显著Gene(包括TF)之间的调控关系(G2G)  >> Gene所在的模块之间的调控关系(M2M)
G2G <- TF_Gene[TF_Gene$Gene %in% sigGene & TF_Gene$TF %in% sigGene,]
G2G <- G2G[-which(G2G$TF==G2G$Gene),]
dim(G2G)

G2G$TF <- G2M_map[G2G$TF]
G2G$Gene <- G2M_map[G2G$Gene]

# 1.1 barplot: sigGene distribution in module (G2G)
GN <- geneInfo[geneInfo$ID %in% sigGene,c(2:3,ncol(geneInfo))]
GN$moduleID <- factor(GN$moduleID,levels = moduleID)

ggplot() + geom_bar(data = GN, aes(x = moduleID), stat = "count",fill= names(moduleID))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = 12))+
  theme(axis.text.y = element_text(face = 'bold',size = 12))


# 1.2 dotplot: in-out degree plot
OUT=as.data.frame(table(G2G$TF))
IN=as.data.frame(table(G2G$Gene))
D <- merge(IN,OUT,by='Var1',all = T)
D[is.na(D)] <- 0
colnames(D) <- c('moduleID','IN','OUT')
FillCol <- names(moduleID)[which(moduleID %in% D$moduleID)]

ggplot(D, aes(x=IN, y=OUT, colour=FillCol)) + 
  geom_point() + 
  geom_text(aes(label=moduleID), size=6) + 
  scale_fill_discrete(guide=FALSE)+
  geom_abline(a=1, b=1, col="red" )+
  guides(fill=FALSE)+
  theme_bw(base_size = 15)



M2M <- G2G[-which(G2G$TF==G2G$Gene),]
M2M <- M2M[!duplicated(M2M[,c(1,2)]),]
dim(M2M)
write.table(M2M,'/home/yzj/JingMA/res/WGCNA/Control/M2M.txt',row.names = FALSE,col.names = TRUE,sep='\t',quote = F)



## 2. 显著TF之间的调控关系(TF2TF)  >> TF所在的模块之间的调控关系(M2M)
TF2TF <- TF_Gene[TF_Gene$Gene %in% sigTF & TF_Gene$TF %in% sigTF,]
dim(TF2TF)
TF2TF <- TF2TF[-which(TF2TF$TF==TF2TF$Gene),]
dim(TF2TF)

TF2TF$TF <- G2M_map[TF2TF$TF]
TF2TF$Gene <- G2M_map[TF2TF$Gene]
write.table(TF2TF,'/home/yzj/JingMA/res/WGCNA/Control/TF2TF.txt',row.names = FALSE,col.names = TRUE,sep='\t',quote = F)

# 2.1 barplot: sigTF distribution in module (TF2TF)
GN <- geneInfo[geneInfo$ID %in% sigTF,c(2:3,ncol(geneInfo))]
GN$moduleID <- factor(GN$moduleID,levels = moduleID[moduleID %in% unique(GN$moduleID)])

ggplot() + geom_bar(data = GN, aes(x = moduleID), stat = "count",
                    fill= names(moduleID[moduleID %in% unique(GN$moduleID)]))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = 12))+
  theme(axis.text.y = element_text(face = 'bold',size = 12))

tmp <- c()
for(i in 1:length(sigTF)){
  tf <-sigTF[i]
  mod <- geneInfo$moduleID[geneInfo$ID ==sigTF[i]]
  tmp[i] <- mod
  names(tmp)[i] <- tf
  #print(paste(tf,mod,sep=':'))
}

Gene='TNFAIP3'
upTFs <- TF_Gene$TF[TF_Gene$Gene %in% Gene,]
for(i in 1:lenght(upTFs)){
  print()
}





# 2.2 dotplot: in-out degree plot
OUT=as.data.frame(table(TF2TF$TF))
IN=as.data.frame(table(TF2TF$Gene))
D <- merge(IN,OUT,by='Var1',all = T)
D[is.na(D)] <- 0

colnames(D) <- c('moduleID','IN','OUT')
FillCol <- names(moduleID)[which(moduleID %in% D$moduleID)]

ggplot(D, aes(x=IN, y=OUT, colour=FillCol)) + 
  geom_point() + 
  geom_text(aes(label=moduleID), size=5) + 
  scale_fill_discrete(guide=FALSE)+
  theme_bw(base_size = 15)+
  guides(fill=FALSE)+
  geom_abline(a=1, b=1, col="red" )


## 3. 显著Gene的UP TF之间的调控关系(upTF2TF)  >> upTF所在的模块之间的调控关系(M2M)
upTF <- unique(TF_Gene$TF[TF_Gene$Gene %in% sigGene])
print(length(upTF))

upTF2TF <- TF_Gene[TF_Gene$Gene %in% upTF & TF_Gene$TF %in% upTF,]
upTF2TF <- upTF2TF[-which(upTF2TF$TF==upTF2TF$Gene),]
dim(upTF2TF)

upTF2TF$TF <- G2M_map[upTF2TF$TF]
upTF2TF$Gene <- G2M_map[upTF2TF$Gene]

# 3.1 barplot: upTF distribution in module (TF2TF)
GN <- geneInfo[geneInfo$ID %in% upTF,c(2:3,ncol(geneInfo))]
GN$moduleID <- factor(GN$moduleID,levels = moduleID[moduleID %in% unique(GN$moduleID)])

ggplot() + geom_bar(data = GN, aes(x = moduleID), stat = "count",
                    fill= names(moduleID[moduleID %in% unique(GN$moduleID)]))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = 12))+
  theme(axis.text.y = element_text(face = 'bold',size = 12))

# 3.2 dotplot: in-out degree plot
OUT=as.data.frame(table(upTF2TF$TF))
IN=as.data.frame(table(upTF2TF$Gene))
D <- merge(IN,OUT,by='Var1',all = T)
D[is.na(D)] <- 0

colnames(D) <- c('moduleID','IN','OUT')
FillCol <- names(moduleID)[which(moduleID %in% D$moduleID)]

ggplot(D, aes(x=IN, y=OUT, colour=FillCol)) + 
  geom_point() + 
  geom_text(aes(label=moduleID), size=6) + 
  scale_fill_discrete(guide=FALSE)+
  theme_bw(base_size = 15)+
  guides(fill=FALSE)+
  geom_abline(a=1, b=1, col="red" )



## 4. 所有TF之间的调控关系(bgTF2TF)  >> bgTF所在的模块之间的调控关系(M2M)
bgTF <- TF
bgTF2TF <- TF_Gene[TF_Gene$Gene %in% bgTF & TF_Gene$TF %in% bgTF,]
bgTF2TF <- bgTF2TF[-which(bgTF2TF$TF==bgTF2TF$Gene),]
dim(bgTF2TF)

bgTF2TF$TF <- G2M_map[bgTF2TF$TF]
bgTF2TF$Gene <- G2M_map[bgTF2TF$Gene]

# 4.1 barplot: bgTF distribution in module (TF2TF)
GN <- geneInfo[geneInfo$ID %in% bgTF,c(2:3,ncol(geneInfo))]
GN$moduleID <- factor(GN$moduleID,levels = moduleID[moduleID %in% unique(GN$moduleID)])

ggplot() + geom_bar(data = GN, aes(x = moduleID), stat = "count",
                    fill= names(moduleID[moduleID %in% unique(GN$moduleID)]))+
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1,face = 'bold',size = 12))+
  theme(axis.text.y = element_text(face = 'bold',size = 12))

## 4.2 M2M
M2M <- bgTF2TF[-which(bgTF2TF$TF==bgTF2TF$Gene),]
M2M <- M2M[!duplicated(M2M[,c(1,2)]),]
dim(M2M)
print(sort(table(M2M$TF)))
write.table(M2M,'/home/yzj/JingMA/res/WGCNA/Control/M2M.txt',row.names = FALSE,col.names = TRUE,sep='\t',quote = F)







#=====================================================================================
#  6. UP TF targeted with CEBPD in C
#=====================================================================================
rm(list=ls())
geneInfo<- read.csv('/home/yzj/JingMA/res/WGCNA/Control/geneInfo.csv')
lnames = load(file = "/home/yzj/JingMA/res/WGCNA/Control/moduleID.RData")

####
TF_Gene <- read.table('/home/yzj/publicData/Harmonizome/merge_Transcription_Factor_Targets_copy.txt',header = TRUE,sep='\t',stringsAsFactors = F)
TF_Gene <- TF_Gene[,-3]

gene='CEBPD'
UP_TF <- TF_Gene$TF[TF_Gene$Gene==gene]
print(length(UP_TF))

#### 
C_UP_Module <- c()
for(i in 1:length(UP_TF)){
  up_tf <- UP_TF[i]
  print(up_tf)
  if(up_tf %in% geneInfo$ID){
    C_UP_Module[i] <- geneInfo$moduleColor[geneInfo$ID==up_tf]
  }else{
    C_UP_Module[i] <- NA
  }
  
}
names(C_UP_Module) <- UP_TF
C_UP_Module <- na.omit(C_UP_Module)
print(length(C_UP_Module))


MN <- as.data.frame(table(C_UP_Module))
MN$moduleID <- moduleID[as.character(MN$C_UP_Module)]
colnames(MN) <- c('moduleColor','MN','moduleID')

ggplot(MN,aes(x=moduleID,y=MN,fill=moduleColor))+
  geom_bar(stat="identity")+
  coord_polar(theta="x")+
  scale_fill_manual(values = levels(MN$moduleColor))+
  theme_bw()+
  theme(axis.text = element_text(face = 'bold',size = 12))+
  labs(title="Freq of Module of UP TF targeted CEBPD")+
  theme(plot.title = element_text(hjust = 0.5,face = 'bold'))+
  guides(fill=FALSE)




####################### 
T <- as.data.frame(table(geneInfo$moduleColor[geneInfo$ID %in% UP_TF]))
ALL <- as.data.frame(table(geneInfo$moduleColor))

all(T$Var1==ALL$Var1)
ratio <- T$Freq/ALL$Freq
names(ratio) <- T$Var1

for(i in 1:length(UP_TF)){
  tf <-UP_TF[i]
  mod <- geneInfo$moduleColor[geneInfo$ID ==UP_TF[i]]
  print(paste(tf,mod,sep=':'))
}

