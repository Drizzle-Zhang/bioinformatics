setwd('/home/zy/my_git/bioinformatics/wgk')
library(Seurat)

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
seurat.first <- readRDS(file.seurat)

# feature plot 
# bone morphogenesis
FeaturePlot(seurat.first, features = c('SFRP2', 'COL1A1', 'COCH', 'SERPINH1'))
# 
FeaturePlot(seurat.first, features = c('HSPA1A', 'HSPA1B', 'CRYAB', 'HSP90AA1'))
# response to temperature stimulus
FeaturePlot(seurat.first, features = c('HSPA6', 'HSPA1A', 'HSPA1B', 'DNAJB1'))
# skeletal system development
FeaturePlot(seurat.first, features = c('SFRP2', 'COL1A1', 'COCH', 'CYTL1'))

FeaturePlot(seurat.first, features = c('GATA6'))

# regulation of stem cell differentiation
FeaturePlot(seurat.first, features = c('ANK3', 'HSP90AA1', 'CTNNA1', 'MYOT'))

FeaturePlot(seurat.first, 
            features = c('CFD', 'MMP10', 'COCH', 'BASP1', 'SFRP4', 'FGL2', 'FBLN1', 'DPT'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('HSPA6', 'COL1A1', 'IER5L', 'GADD45G', 'IGFBP5', 'OGN', 
                         'ASPN', 'PI16', 'CRABP2', 'FBLN2', 'CLIC2', 'GPNMB'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('FST', 'CTGF', 'FRZB', 'FGFBP2', 'PLA2G2A', 'S100A1', 
                         'MT1M', 'S100B', 'PRELP'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('LAMB3', 'NOS2', 'LOX', 'ELN', 'COL2A1', 'CX3CL1', 
                         'COL11A2', 'SCUBE3', 'COLGALT2'), 
            ncol = 3)
FeaturePlot(seurat.first, 
            features = c('CYTL1', 'SNHG12', 'HSPB1', 'DLGAP1-AS2', 'NDRG2', 
                         'RASL11B', 'DUSP2', 'HSPA6', 'MDFI', 'GPRC5C'), 
            ncol = 3)

FeaturePlot(seurat.first, 
            features = c('ACTA2', 'ACTG2'), 
            ncol = 2)

##### extracellular matrix
# Fibulin
FeaturePlot(seurat.first, features = c('FBLN1', 'FBLN2', 'FBLN3', 'FBLN4', 'FBLN5', 'FBLN6', 'FBLN7'))
# Fibronectin

# Laminin
FeaturePlot(seurat.first, features = c('LAMA1', 'LAMA2', 'LAMA3', 'LAMA4', 'LAMA5', 
                                       'LAMB1', 'LAMB2', 'LAMB3', 'LAMB4',
                                       'LAMC1', 'LAMC2', 'LAMC3'))


# Thrombospondin
FeaturePlot(seurat.first, features = c('THBS1', 'THBS2', 'THBS3', 'THBS4'))


##### growth factor binding
# IGF
FeaturePlot(seurat.first, features = c('IGFBP1', 'IGFBP2', 'IGFBP3', 'IGFBP4', 'IGFBP5', 'IGFBP6', 'IGFBP7'))


##### COL
FeaturePlot(seurat.first, features = c('COL4A1', 'COL4A2', 'COL18A1', 'COL5A2'))


##### genes in cluster
# 9
unlist(strsplit((list.kegg$`9`)['hsa05022', 'geneID'], '/'))
list.marker$`9`[df.geneid2symbol[unlist(strsplit((list.kegg$`9`)['hsa04350', 'geneID'], '/')),],]
FeaturePlot(seurat.first, features = c('ID1', 'ID2', 'ID3', 'ID4'))
FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL14A1', 'COL3A1'))
FeaturePlot(seurat.first, features = c('CXCL14', 'CXCL12', 'FGF18', 'OGN'))
FeaturePlot(seurat.first, features = c('FGFR1', 'FGFR2', 'FGFR3', 'FGFR4'))
FeaturePlot(seurat.first, features = c('FBLN1', 'FBLN2', 'FBLN3', 'FBLN4', 'FBLN5', 'FBLN6', 'FBLN7'))
FeaturePlot(seurat.first, features = c('SKIL', 'KLF4'))
FeaturePlot(seurat.first, features = c('FOS', 'FOSB'))

FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL3A1'), ncol = 3)
FeaturePlot(seurat.first, features = c('ID3', 'SKIL', 'MYC'), ncol = 3)
FeaturePlot(seurat.first, features = c('CLK1', "MAP3K2", "MAP3K8"), ncol = 3)
FeaturePlot(seurat.first, features = c('S100A4', 'GPNMB', "CBX4", "ATF3"))
FeaturePlot(seurat.first, features = c('FOS', 'FOSB', "JUN", "JUNB"))
FeaturePlot(seurat.first, features = c('FBLN1', 'FBLN2'), ncol = 3)
FeaturePlot(seurat.first, features = c('ASPN', 'OGN'))

FeaturePlot(seurat.first, features = c('CDKN1A', 'CDKN1B', 'CDKN1C', 
                                       'CDKN2A', 'CDKN2B', 'CDKN2C', 'CDKN2D',
                                       'CDKN3',
                                       'CDK1', 'CDK2', 'CDK3', 'CDK4', 'CDK5', 'CDK6'))
FeaturePlot(seurat.first, features = c('CCNA1', 'CCNA2', 'CCNB1', 'CCNB2', 'CCNB3', 
                                       'CCND1', 'CCND2', 'CCND3', 'CCNE1', 'CCNE2'))
FeaturePlot(seurat.first, features = c('PLK1', 'E2F1', 'FOXM1', 'MKI67', 'BUB1', 'TOP2A'), ncol = 3)

VlnPlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL3A1', 
                                   'ID3', 'SKIL', 'MYC',
                                   'CLK1', "MAP3K2", "MAP3K8", 'S100A4', 'GPNMB', "CBX4", "ATF3",
                                   'FOS', 'FOSB', "JUN", "JUNB",
                                   'FBLN1', 'FBLN2', 'ASPN', 'OGN', 'DCN'), 
        group.by = 'RNA_snn_res.0.8', pt.size = 0, ncol = 6)



# 12
FeaturePlot(seurat.first, features = c('ACTA2', 'ACTG2', 'MCAM'), ncol = 3)
FeaturePlot(seurat.first, features = c('MYH11', 'MYL9', "CNN1", "SYNPO2", 'LMOD1', 'CAP2'), ncol = 3)
VlnPlot(seurat.first, features = c('ACTA2', 'ACTG2', 'MCAM', 
                                   'MYH11', 'MYL9', "CNN1", "SYNPO2", 'LMOD1', 'CAP2', 'CDH6', 'PLN'), 
        group.by = 'RNA_snn_res.0.8', pt.size = 0, ncol = 3)


# 11
FeaturePlot(seurat.first, features = c('MCAM', 'SYNM', 'ADRA2A', 'FILIP1', 'SRL', 'SORBS2', 
                                       'FSTL3', 'PLN', 'SLC25A19'))
FeaturePlot(seurat.first, features = c('APOLD1', 'SYNM', 'SYNPO2', 'OLFM3', 'SRL', 'SORBS2', 
                                       'FSTL3', 'PLN', 'SLC25A19'))
VlnPlot(seurat.first, features = c('APOLD1', 'SYNM', 'SYNPO2', 'OLFM3', 'SRL', 'SORBS2', 
                                   'FSTL3', 'PLN', 'SLC25A19'), 
        group.by = 'RNA_snn_res.0.8', pt.size = 0, ncol = 3)


# 7
FeaturePlot(seurat.first, features = c('EDNRB', 'APOLD1', 'PDGFRB', 'SYNPO2', 'OLFM2', 'PDE3A'), ncol = 3)
VlnPlot(seurat.first, features = c('EDNRB', 'APOLD1', 'PDGFRB', 'SYNPO2', 'OLFM2', 'PDE3A'), 
        group.by = 'RNA_snn_res.0.8', pt.size = 0, ncol = 3)


# 10
FeaturePlot(seurat.first, features = c('APOE', 'CFD', 'CTSC', 'MEDAG', 'SFRP4', 'VCAN', 
                                       'TYMP', 'ABI3BP', 'CHSY1'))
FeaturePlot(seurat.first, features = row.names(list.marker$`10`)[1:20], ncol = 5)


# 4
FeaturePlot(seurat.first, features = c('APOE', 'CFD', 'CTSC', 'MEDAG', 'SFRP4', 'VCAN', 
                                       'TYMP', 'ABI3BP', 'CHSY1'))
FeaturePlot(seurat.first, features = row.names(list.marker$`4`)[1:20], ncol = 5)
FeaturePlot(seurat.first, features = row.names(list.marker$`4`)[21:40], ncol = 5)



# 16
FeaturePlot(seurat.first, features = c('ANK3', 'NTM', "SLITRK6", 'CIT', 'NRXN1', "NGFR",
                                       'SORCS1', 'PLP1', 'PMP2'), ncol = 3)
VlnPlot(seurat.first, features = c('ANK3', 'NTM', "SLITRK6", 'CIT', 'NRXN1', "NGFR",
                                   'SORCS1', 'PLP1', 'PMP2'), 
        group.by = 'RNA_snn_res.0.8', pt.size = 0, ncol = 3)

# 2/3
FeaturePlot(seurat.first, features = c('COL2A1', 'COL9A3', 'COL11A2', 'COL9A2', 'COL11A1'))


# 0/1
FeaturePlot(seurat.first, features = row.names(list.marker$`0`)[1:20], ncol = 5)

# 7


