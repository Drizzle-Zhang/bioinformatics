setwd('/home/zy/my_git/bioinformatics/wgk')
# .libPaths('/home/yzj/R/x86_64-pc-linux-gnu-library/4.0')
library(Seurat)
library(ggplot2)
require("RColorBrewer")

# FeaturePlot(seurat.all.filter, c('COL9A3', 'COL9A2', 'COL8A2', 'COL8A1', 
#                                  'COL6A2', 'COL6A1','COL11A1', 'COL11A2', 'COL24A1'))

sel.pair <- 
    c('TGFB1_TGFbeta receptor1', 
      # 'TGFB1_TGFbeta receptor2', 'TGFB1_TGFBR3',
      'TGFB2_TGFbeta receptor1', 
      # 'TGFB2_TGFbeta receptor2', 'TGFB2_TGFBR3',
      # 'TGFB3_TGFbeta receptor1', 'TGFB3_TGFbeta receptor2', 'TGFB3_TGFBR3',
      # 'FGFR1_FGF7', 'FGF7_FGFR3', 'FGF7_FGFR2', 'FGF2_FGFRL1', 
      # 'FGF2_FGFR3', 'FGF2_FGFR2', 'FGF2_FGFR1', 'FGF2_CD44',
      # 'IGF2_IGF2R', 'IGF2_IGF1R', 
      # 'BMPR1A_BMPR2_BMP2', 'BMPR1B_BMPR2_BMP2', 'BMR1A_ACR2A_BMP2', 'BMR1B_AVR2A_BMP2',
      # 'CSF1_SLC7A1', 'CSF1_SIRPA', 
      'COL9A3_a2b1 complex', 'COL9A3_a1b1 complex', 'COL9A3_a10b1 complex', 
      'COL9A2_a2b1 complex', 'COL9A2_a1b1 complex', 'COL9A2_a10b1 complex',
      # 'COL8A1_a10b1 complex', 'COL8A1_a1b1 complex', 'COL8A1_a2b1 complex',
      # 'COL8A2_a10b1 complex', 'COL8A2_a1b1 complex', 'COL8A2_a2b1 complex',
      # 'COL6A1_a10b1 complex', 'COL6A1_a1b1 complex', 'COL6A1_a2b1 complex',
      # 'COL6A2_a10b1 complex', 'COL6A2_a1b1 complex', 'COL6A2_a2b1 complex',
      'COL2A1_a10b1 complex', 'COL2A1_a1b1 complex', 'COL2A1_a2b1 complex',
      # 'COL1A1_a10b1 complex', 'COL1A1_a1b1 complex', 'COL1A1_a2b1 complex',
      # 'COL1A2_a10b1 complex', 'COL1A2_a1b1 complex', 'COL1A2_a2b1 complex',
      'COL11A1_a10b1 complex', 'COL11A1_a1b1 complex', 'COL11A1_a2b1 complex',
      # 'COL11A2_a10b1 complex', 'COL11A2_a1b1 complex', 'COL11A2_a2b1 complex',
      'FN1_aVb5 complex', 'FN1_aVb1 complex', 
      # 'FN1_aVb3 complex', 
      'FN1_a5b1 complex', 'FN1_a3b1 complex', 'FN1_a2b1 complex', 'FN1_a10b1 complex'
      # 'TNF_SEMA4C', 'TNF_RIPK1', 'TNF_PTPRS', 'TNF_NOTCH1', 'TNF_FAS'
      )

sel.pair.LR <- 
  c('TGFB1_TGFbeta receptor1', 
    # 'TGFB1_TGFbeta receptor2', 'TGFB1_TGFBR3',
    'TGFB2_TGFbeta receptor1', 
    # 'TGFB2_TGFbeta receptor2', 'TGFB2_TGFBR3',
    # 'TGFB3_TGFbeta receptor1', 'TGFB3_TGFbeta receptor2', 'TGFB3_TGFBR3',
    # 'FGF7_FGFR1', 'FGF7_FGFR3', 'FGF7_FGFR2', 'FGF2_FGFRL1', 
    # 'FGF2_FGFR3', 'FGF2_FGFR2', 'FGF2_FGFR1', 'FGF2_CD44',
    # 'IGF2_IGF2R', 'IGF2_IGF1R', 
    # 'BMP2_BMPR1A_BMPR2', 'BMP2_BMPR1B_BMPR2', 'BMP2_BMR1A_ACR2A', 'BMP2_BMR1B_AVR2A',
    # 'CSF1_SLC7A1', 'CSF1_SIRPA', 
    'COL9A3_a2b1 complex', 'COL9A3_a1b1 complex', 'COL9A3_a10b1 complex', 
    'COL9A2_a2b1 complex', 'COL9A2_a1b1 complex', 'COL9A2_a10b1 complex',
    # 'COL8A1_a10b1 complex', 'COL8A1_a1b1 complex', 'COL8A1_a2b1 complex',
    # 'COL8A2_a10b1 complex', 'COL8A2_a1b1 complex', 'COL8A2_a2b1 complex',
    # 'COL6A1_a10b1 complex', 'COL6A1_a1b1 complex', 'COL6A1_a2b1 complex',
    # 'COL6A2_a10b1 complex', 'COL6A2_a1b1 complex', 'COL6A2_a2b1 complex',
    'COL2A1_a10b1 complex', 'COL2A1_a1b1 complex', 'COL2A1_a2b1 complex',
    # 'COL1A1_a10b1 complex', 'COL1A1_a1b1 complex', 'COL1A1_a2b1 complex',
    # 'COL1A2_a10b1 complex', 'COL1A2_a1b1 complex', 'COL1A2_a2b1 complex',
    'COL11A1_a10b1 complex', 'COL11A1_a1b1 complex', 'COL11A1_a2b1 complex',
    # 'COL11A2_a10b1 complex', 'COL11A2_a1b1 complex', 'COL11A2_a2b1 complex',
    'FN1_aVb5 complex', 'FN1_aVb1 complex', 
    # 'FN1_aVb3 complex', 
    'FN1_a5b1 complex', 'FN1_a3b1 complex', 'FN1_a2b1 complex', 'FN1_a10b1 complex'
    # 'TNF_SEMA4C', 'TNF_RIPK1', 'TNF_PTPRS', 'TNF_NOTCH1', 'TNF_FAS'
    )
pairs.LR <- data.frame(pair = sel.pair, LR = sel.pair.LR)

# normal
file.normal <- '/home/yzj/JingMA_NEW/res/CellPhoneDB_alone_subcelltype/out_Normal/LR_filter.txt'
df.normal <- read.delim(file.normal)
df.normal <- df.normal[df.normal$pair %in% sel.pair, ]
df.normal <- merge(df.normal, pairs.LR, by = 'pair')
df.normal$cell1 <- rep('1', nrow(df.normal))
df.normal$cell2 <- rep('2', nrow(df.normal))
for (i in rownames(df.normal)) {
    name.pair <- df.normal[i, 'pair']
    name.LR <- df.normal[i, 'LR']
    cells <- strsplit(df.normal[i, 'clusters'], split = '|', fixed = T)[[1]]
    if (name.pair == name.LR) {
        df.normal[i, 'cell1'] <- cells[1]
        df.normal[i, 'cell2'] <- cells[2]
    } else {
        df.normal[i, 'cell1'] <- cells[2]
        df.normal[i, 'cell2'] <- cells[1]
    }
}
df.normal$log10Pval <- -log10(df.normal$pvalue)
df.normal$LR <- factor(df.normal$LR, levels = sel.pair.LR)
df.normal$cell1[df.normal$cell1 == 'TC'] <- 'C0'
df.normal$cell2[df.normal$cell2 == 'TC'] <- 'C0'
df.normal <- df.normal[df.normal$cell1 %in% c("C2", "C1", "C0", "CSC"),]
df.normal <- df.normal[df.normal$cell2 %in% c("C2", "C1", "C0", "CSC"),]
df.normal$cell1 <- factor(df.normal$cell1,
                          levels = c("CSC", "C0", "C1", "C2"))
df.normal$cell2 <- factor(df.normal$cell2,
                          levels = c("CSC", "C0", "C1", "C2"))
df.normal$log2mean[df.normal$log2mean > 0.7] <- 0.7
df.normal$log2mean[df.normal$log2mean < -2.3] <- NA
# df.normal <- df.normal[!is.na(df.normal$log2mean),]
cell.labels <- c('CSPC', 'EC', 'IC', 'LC')
names(cell.labels) <- c('CSC', 'C0', 'C1', 'C2')

plot.dot.normal <- 
    ggplot(df.normal, aes(x = cell2, y = LR, size = log10Pval, color = log2mean)) + 
    geom_point(fill = 'cornsilk') + 
    facet_grid( ~ cell1, labeller = labeller(cell1 = cell.labels)) +
    scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(100), na.value = 'white') +
    scale_size_continuous(range = c(2.5,4.5), breaks = c(3, 5, 7)) +
    scale_x_discrete(breaks = c('CSC', 'C0', 'C1', 'C2'),
                     labels = c('CSPC', 'EC', 'IC', 'LC')) +
    labs(x = '', y = 'Interacting molecules', 
         color = expression(paste("log"[2], "(mean(molecule1, molecule2))")),
         size = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(panel.background = element_rect(color = 'gray', fill = 'transparent'),
          legend.position = 'bottom',
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 11, color = 'black'),
          axis.text.y = element_text(size = 11, color = 'black'), 
          axis.title.y = element_text(size = 14, face = 'bold', color = 'black')) + 
    guides(size = F)
path.cell.commu <- '/home/disk/drizzle/wgk/cell_commu'
ggsave(plot = plot.dot.normal, path = path.cell.commu, 
       filename = 'cell_commu_normal.pdf',
       height = 13.5, width = 15, units = 'cm')

# microtia
file.microtia <- '/home/yzj/JingMA_NEW/res/CellPhoneDB_alone_subcelltype/out_Microtia/LR_filter.txt'
df.microtia <- read.delim(file.microtia)
df.microtia <- df.microtia[df.microtia$pair %in% sel.pair, ]
df.microtia <- merge(df.microtia, pairs.LR, by = 'pair')
df.microtia$cell1 <- rep('1', nrow(df.microtia))
df.microtia$cell2 <- rep('2', nrow(df.microtia))
for (i in rownames(df.microtia)) {
    name.pair <- df.microtia[i, 'pair']
    name.LR <- df.microtia[i, 'LR']
    cells <- strsplit(df.microtia[i, 'clusters'], split = '|', fixed = T)[[1]]
    if (name.pair == name.LR) {
        df.microtia[i, 'cell1'] <- cells[1]
        df.microtia[i, 'cell2'] <- cells[2]
    } else {
        df.microtia[i, 'cell1'] <- cells[2]
        df.microtia[i, 'cell2'] <- cells[1]
    }
}
df.microtia$log10Pval <- -log10(df.microtia$pvalue)
df.microtia$LR <- factor(df.microtia$LR, levels = sel.pair.LR)
df.microtia$cell1[df.microtia$cell1 == 'TC'] <- 'C0'
df.microtia$cell2[df.microtia$cell2 == 'TC'] <- 'C0'
df.microtia <- df.microtia[df.microtia$cell1 %in% c("C2", "C1", "C0", "CSC"),]
df.microtia <- df.microtia[df.microtia$cell2 %in% c("C2", "C1", "C0", "CSC"),]
df.microtia$cell1 <- factor(df.microtia$cell1,
                          levels = c("CSC", "C0", "C1", "C2"))
df.microtia$cell2 <- factor(df.microtia$cell2,
                          levels = c("CSC", "C0", "C1", "C2"))
df.microtia$log2mean[df.microtia$log2mean > 0.7] <- 0.7
df.microtia$log2mean[df.microtia$log2mean < -2.3] <- NA
# df.microtia <- df.microtia[!is.na(df.microtia$log2mean),]

plot.dot.microtia <- 
    ggplot(df.microtia, aes(x = cell2, y = LR, size = log10Pval, color = log2mean)) + 
    geom_point(fill = 'cornsilk') + 
    facet_grid( ~ cell1, labeller = labeller(cell1 = cell.labels)) +
    scale_colour_gradientn(
        colours = colorRampPalette(rev(brewer.pal(n = 5, name = "RdYlBu")))(100), na.value = 'white') +
    scale_size_continuous(range = c(2.5,4.5), breaks = c(3, 5, 7)) +
    scale_x_discrete(breaks = c('CSC', 'C0', 'C1', 'C2'),
                     labels = c('CSPC', 'EC', 'IC', 'LC')) +
    labs(x = '', y = 'Interacting molecules', 
         color = expression(paste("log"[2], "(mean(molecule1, molecule2))")),
         size = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme(panel.background = element_rect(color = 'gray', fill = 'transparent'),
          legend.position = 'bottom',
          # axis.text.y = element_text(size = 10, color = 'black'), 
          # axis.title.y = element_text(size = 12, face = 'bold', color = 'black'),
          axis.text.y = element_blank(),
          axis.title.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.text.x = element_text(angle = 90, vjust = 1, hjust = 1, size = 11, 
                                     color = 'black')) + 
  guides(color = F)
ggsave(plot = plot.dot.microtia, path = path.cell.commu, 
       filename = 'cell_commu_microtia.pdf',
       height = 13, width = 9.5, units = 'cm')

