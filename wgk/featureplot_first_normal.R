setwd('/home/zy/my_git/bioinformatics/wgk')
library(Seurat)

file.seurat <- '/home/disk/drizzle/wgk/data/ear.seurat.first.Rdata'
seurat.first <- readRDS(file.seurat)

# feature plot 
# 1
# df.marker.1 <- list.marker$`1`[list.marker$`1`]
FeaturePlot(seurat.first, features = row.names(list.marker$`1`)[1:20], ncol = 5)
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

##### FOS
FeaturePlot(seurat.first, features = c('FOS', 'FOSB'))


##### genes in cluster
# 9
unlist(strsplit((list.kegg$`9`)['hsa05022', 'geneID'], '/'))
list.marker$`9`[df.geneid2symbol[unlist(strsplit((list.kegg$`9`)['hsa04350', 'geneID'], '/')),],]
FeaturePlot(seurat.first, features = c('ID1', 'ID2', 'ID3', 'ID4'))
FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL14A1', 'COL3A1'))
FeaturePlot(seurat.first, features = c('CXCL14', 'CXCL12', 'FGF18', 'OGN'))
FeaturePlot(seurat.first, features = c('FGFR1', 'FGFR2', 'FGFR3', 'FGFR4'))
FeaturePlot(seurat.first, features = c('FBLN1', 'FBLN2', 'FBLN3', 'FBLN4', 'FBLN5', 'FBLN6', 'FBLN7'))
FeaturePlot(seurat.first, features = c('IGF1', 'IGF2'))

FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL14A1', 'COL3A1'))

FeaturePlot(seurat.first, features = c('COL1A1', 'COL1A2', 'COL14A1', 'COL3A1',
                                       'ID3', 'CRABP2',
                                       'FBLN1', 'FBLN2', 'ASPN', 'OGN',
                                       'CXCL14', 'CXCL12'))
FeaturePlot(seurat.first, features = c('CLK1', 'COL1A2', 'COL14A1', 'COL3A1',
                                       'ID3', 'CRABP2',
                                       'FBLN1', 'FBLN2', 'ASPN', 'OGN',
                                       'CXCL14', 'CXCL12'))


# 12
FeaturePlot(seurat.first, features = c('CDKN1B', 'CDKN2B', 'CDKN1C', 'CDKN2D', 'CDKN2C'))
FeaturePlot(seurat.first, features = c('ACTA2', 'MYH11', 'CNN1', 'SYNPO2', 'ACTG2', 'IGFBP5'))

# 11
FeaturePlot(seurat.first, features = c('MCAM', 'HES4', 'PLAU', 'KCNE4', 'ADRA2A', 'PI15'))


