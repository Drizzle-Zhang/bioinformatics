library(ggplot2)

file.SNP <- '/mdshare/node9/zy/Reg_brain/SNP_aceDHS/results_snp.txt'


file.Gene <- '/mdshare/node9/zy/Reg_brain/gene_bed/results_disease_gene.txt'
df.Gene <- read.delim(file.Gene)
df.Gene <- df.Gene[df.Gene$pvalue < 1,]
df.Gene$logp <- -log10(df.Gene$pvalue)
sel.disease <- c('Chromosome 17p Deletion Syndrome', 'Gaucher Disease, Type 2 (disorder)',
                 'Diastematomyelia', 'Tethered Cord Syndrome', 'Acrania',
                 'Exencephaly', 'Iniencephaly', 'Spinal Cord Myelodysplasia',
                 'Craniorachischisis', 'Neural Tube Defects', 'HOLOPROSENCEPHALY 4 (disorder)',
                 'Digitorenocerebral Syndrome', 'Spastic paraplegia 17',
                 'Cerebellar Ataxia, Early Onset', 'Cerebellar Ataxia, Late Onset',
                 'KRABBE DISEASE, ATYPICAL, DUE TO SAPOSIN A DEFICIENCY',
                 'KLEEFSTRA SYNDROME 1', 'MENTAL RETARDATION, X-LINKED 3',
                 'Corticostriatal-Spinal Degeneration', 'Cerebellar Degenerations, Primary',
                 'AICARDI-GOUTIERES SYNDROME 3', 'SCAPHOCEPHALY, MAXILLARY RETRUSION, AND MENTAL RETARDATION',
                 'Spinocerebellar Degeneration', 'Rachischisis', 'Spina Bifida', 'Status Dysraphicus')
df.Gene <- df.Gene[df.Gene$disease_name %in% sel.disease,]
df.Gene$disease_name <- factor(df.Gene$disease_name, levels = rev(sel.disease))

plot_bar <- 
    ggplot(df.Gene, aes(x = disease_name, y = logp)) + 
    geom_bar(stat = 'identity', color = 'dimgray', fill = 'dimgray', 
             position=position_dodge(0.8)) + 
    coord_flip() + 
    labs(x = 'Disease', y = expression(paste("-log"[10], "(", italic("P"), "-value)"))) +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.text.x = element_text(size = 11, color = 'black', family = 'Arial'),
          axis.text.y = element_text(size = 8, color = 'black', family = 'Arial'),
          axis.title = element_text(size = 12, family = 'Arial'))
ggsave(filename = 'Gene_plot.png',
       path = '/mdshare/node9/zy/Reg_brain/plot', plot = plot_bar,
       units = 'cm', height = 14, width = 18)
