library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.out.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)


topN <- 50
path.tree <- '/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan'

level = 'genus'
# file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot/Top30_', 
#               level, '.csv')
file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/abundance/relative_abundance/', 
              level, '.csv')
genus <- read.table(file = file, header = TRUE, row.names = 1,
                    sep = '\t', stringsAsFactors = FALSE,
                    check.names = FALSE)
rownames(genus) <- paste0('g__', rownames(genus))
level = 'species'
file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/abundance/relative_abundance/', 
              level, '.csv')
species <- read.table(file = file, header = TRUE, row.names = 1,
                    sep = '\t', stringsAsFactors = FALSE,
                    check.names = FALSE)
rownames(species) <- paste0('s__', rownames(species))

path.input <- paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot_',
                    level)
path.plot <- '/home/drizzle_zhang/microbiome/result/Figs/'

sel.species <- c('s__Lactobacillus_gasseri', 's__Lactobacillus_murinus', 
                 's__Helicobacter_sp._MIT_03-1614')
sel.genus <- c('g__Lactobacillus', 'g__Clostridium_sensu_stricto_1',
               'g__Helicobacter', 'g__Turicibacter', 'g__Streptococcus', 
               'g__Ruminococcaceae_UCG-014', 'g__Prevotella_1', 'g__Prevotella_9',
               'g__Treponema_2')

# mean abundance
merge.abundance <- rbind(genus[sel.genus, ], species[sel.species, ])
mean.abundance <- rowMeans(merge.abundance)
sort(mean.abundance)
sort.abundance <- c('g__Lactobacillus', 's__Lactobacillus_murinus',
                    's__Lactobacillus_gasseri', 'g__Ruminococcaceae_UCG-014', 
                    'g__Prevotella_9', 'g__Helicobacter', 
                    's__Helicobacter_sp._MIT_03-1614', 'g__Treponema_2',
                    'g__Turicibacter', 'g__Clostridium_sensu_stricto_1',
                    'g__Prevotella_1', 'g__Streptococcus')

list.vec.dose <- list(c(0, 1, 2, 3), c(0, 1), c(0, 2), c(0, 3))
vec.gender <- c('male', 'female')
series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
for (vec.dose in list.vec.dose) {
    # dose
    str_dose <- paste0(as.character(vec.dose), collapse = '')
    path.dose <- paste(path.tree, str_dose, sep = '/')
    # file.tree <- paste0(path.dose, '/tree_backbone_', topN, '.txt')
    # df.tree <- read.table(file.tree, sep = '.')
    # colnames(df.tree) <- c('phylum', 'class', 'order', 'family', 'genus', 'OTU_ID')
    # sel.species <- as.character(unique(df.tree[, colnames(df.tree) == level]))
    # sel.species <- unlist(lapply(sel.species, function(x) {substr(x, 4, nchar(x))}))
    
    df.combine <- data.frame(stringsAsFactors = F)
    for (gender in vec.gender) {
        # gender
        # gender = 'male'
        df.meta.gender <- df.meta[df.meta$Gender == gender, ]
        df.plot.point <- data.frame()
        for (sub.time in series.time) {
            # select meta
            sel.meta <- df.meta.gender[df.meta.gender$Time == sub.time,]
            sel.meta <- sel.meta[sel.meta$Dose %in% vec.dose,]
            row.names(sel.meta) <- sel.meta$Sample
            
            # select sample
            use.sample <- sel.meta$SampleName
            sub.genus <- genus[sel.genus, use.sample]
            sub.species <- species[sel.species, use.sample]
            sub.abundance <- rbind(sub.genus, sub.species)
            diff.in <- cbind(t(sub.abundance), sel.meta)
            types <- row.names(sub.abundance)
            df.out <- data.frame()
            interval <- 0.02
            for (type in types) {
                sub.in <- diff.in[, c(type, 'Group')]
                fold.change <- log10(mean(sub.in[sub.in$Group == 'Treat', type]) / 
                                         mean(sub.in[sub.in$Group == 'Control', type]))
                colnames(sub.in) <- c('species', 'Group')
                out.wilcox <-
                    wilcox.test(as.formula('species ~ Group'), data = sub.in)
                # a = sub.in[sub.in$Group == 'Treat', type]
                # b = sub.in[sub.in$Group == 'Control',type]
                # out.wilcox <- ks.test(a, b, data = sub.in)
                p.value <- out.wilcox$p.value
                p.value.round <- round(out.wilcox$p.value, 4)
                # if (is.na(fold.change)) {
                #     
                # }
                # if (fold.change < 0) {
                #     position <- fold.change - interval
                # } else {
                #     position <- fold.change + interval
                # }
                df.out <- rbind(df.out, 
                                data.frame(type = type, fold.change = fold.change, 
                                           p.value = p.value, 
                                           # position = position,
                                           time = sub.time, p.value.round = p.value.round))
            }
            df.plot.point <- rbind(df.plot.point, df.out)
        }
        df.plot.point$log.pval <- -log10(df.plot.point$p.value)
        # file.bubble <- paste0(path.input, '/', gender, '_Bubble_',
        #                       paste0(as.character(vec.dose), collapse = ''), '.txt')
        # df.sub <- read.delim(file.bubble)
        df.plot.point$Gender <- rep(gender, nrow(df.plot.point))
        df.combine <- rbind(df.combine, df.plot.point)
    }
    # df.combine$Taxonomy <- fct_inorder(df.combine$type)
    # df.combine$Taxonomy <- (df.combine$type)
    df.combine$logPval <- df.combine$log.pval
    df.combine$logPval[df.combine$p.value > 0.25] <- NA
    df.combine$fold.change[df.combine$fold.change > 0.95] <- 0.95
    df.combine$fold.change[df.combine$fold.change < -0.95] <- -0.95
    df.combine$direct.logp <- df.combine$fold.change*df.combine$logPval
    df.combine$direct.logp[df.combine$direct.logp > 0] <- 1
    df.combine$direct.logp[df.combine$direct.logp < 0] <- -1
    # delete some species
    # final.species <- data.frame(stringsAsFactors = F)
    # for (spe in sel.species) {
    #     sub.male <- df.combine[df.combine$type == spe & df.combine$Gender == 'male',]
    #     sub.female <- df.combine[df.combine$type == spe & df.combine$Gender == 'female',]
    #     if (abs(sum(sub.male$direct.logp, na.rm = T)) > 2 | abs(sum(sub.female$direct.logp, na.rm = T)) > 2) {
    #         final.species <- rbind(final.species, 
    #                                data.frame(type = spe, 
    #                                           score = sum(sub.male$direct.logp, na.rm = T)))
    #         # final.species <- c(final.species, spe)
    #     }
    # }
    # df.final <- merge(df.combine, final.species, by = 'type')
    # order_species <- rownames(class)
    # final.species <- order_species[order_species %in% df.final$type]
    # final.species <- final.species$type[order(final.species$score, decreasing = F)]
    df.final <- df.combine
    df.final$Taxonomy <- factor(df.final$type, levels = sort.abundance)
    plot.bubble <- 
        ggplot(data = df.final, 
               aes(x = time, y = Gender, size = logPval, color = fold.change)) + 
        geom_point(fill = 'cornsilk') + 
        scale_colour_gradient2(low = '#00008B', mid = "white", high = '#8B0000', na.value = 'white') +
        # scale_colour_gradient2(low = 'blue', mid = "white", high = 'red', na.value = 'white') + 
        scale_size_continuous(range = c(2, 5), breaks = c(1, 2, 3, 4)) +
        facet_grid(Taxonomy ~ ., scales = 'free', space = 'free') +
        labs(x = 'Time', y = 'Genus', color = 'logFC', 
             size = expression(paste("-log"[10], "(adj", italic("P"), "-value)"))) + 
        theme_bw() +
        theme(
            panel.background = element_rect(color = 'black', size = 1,
                                            fill = 'transparent'),
            panel.grid = element_blank(),
            strip.background = element_rect(
                color = 'transparent', fill = 'transparent'),
            strip.text.y = element_text(angle = 0, color = 'black', hjust = 0,
                                        size = 12, family = 'Arial'),
            axis.text.x = element_text(size = 10, family = 'Arial'),
            axis.title.x = element_text(size = 12, family = 'Arial'),
            legend.title = element_text(size = 12, family = 'Arial'),
            legend.text = element_text(size = 10, family = 'Arial'), 
            axis.text.y = element_blank(),
            axis.title.y =  element_blank(),
            axis.ticks.y = element_blank()
            )
    sub.genders <- rep(c("Male","Female"), times = length(final.species))
    sub.genus <- rep(rev(final.species), c(rep(2, length(final.species))))
    sort.row <- paste(sub.genders, sub.genus, sep = '_')
    df.status <- data.frame(status=sub.genders,
                            p = rep('white', length(sort.row)), 
                            row_name = sort.row, 
                            cluster = sub.genus)
    df.status$row_name <- factor(df.status$row_name, levels = rev(sort.row))
    df.status$cluster <- factor(df.status$cluster, levels = final.species)
    df.status$status <- factor(df.status$status, levels = c("Male","Female"))
    plot.Status <- 
        ggplot(data = df.status, aes(x=p,y=row_name,fill=status))+
        geom_tile() + 
        facet_grid(cluster ~ ., scales = 'free', space = 'free') +
        scale_y_discrete(position="right") +
        scale_fill_manual(breaks = c("Male","Female"), 
                          values = c("#BC80BD", "#5AB4AC")) +
        # scale_fill_manual(values = c(muted("blue"), muted("red"))) +
        xlab(NULL) + ylab(NULL) +
        theme(panel.background = element_rect(color = 'transparent',
                                              fill = 'transparent'),
              strip.text.y = element_blank(),
              axis.text = element_blank(),
              axis.ticks = element_blank(),
              panel.spacing = unit(0.1, "cm"),
              legend.title = element_text(size = 12, family = 'Arial'),
              legend.text = element_text(size = 10, family = 'Arial'), 
              legend.position = 'left',
              plot.margin = unit(c(0, -1, 0.6, 0.6),"cm")) +
        labs(fill = "Status")
    plot.final <- plot.Status + plot.bubble + plot_layout(widths = c(1, 30),
                                                           guides = 'collect')
    
    ggsave(plot = plot.final, path = path.plot, 
           filename = paste0('Bubble_', 
                             paste0(as.character(vec.dose), collapse = ''), '.png'),
           height = 16, width = 28, units = 'cm')
    
}


# abundance
df.abundance <- data.frame(stringsAsFactors = F)
for (gender in vec.gender) {
    # gender
    # gender = 'male'
    df.meta.gender <- df.meta[df.meta$Gender == gender, ]
    use.sample <- df.meta.gender$SampleName
    merge.abundance <- rbind(genus[sel.genus, use.sample], 
                             species[sel.species, use.sample])
    mean.abundance <- rowMeans(merge.abundance)*100
    df.abundance <- rbind(df.abundance, 
                          data.frame(species = sort.abundance, 
                                     abundance = mean.abundance[sort.abundance],
                                     gender = rep(gender, length(sort.abundance))))
}
df.abundance$species <- factor(df.abundance$species, levels = rev(sort.abundance))
df.abundance$gender <- factor(df.abundance$gender, levels = rev(c('male', 'female')))

plot.abundance <- 
    ggplot(data = df.abundance, 
                         aes(x = species, y = abundance, fill = gender))  +
    geom_bar(position = 'dodge', stat = 'identity') + 
    scale_fill_manual(breaks = c("male","female"), 
                      values = c("#BC80BD", "#5AB4AC")) +
    # scale_x_discrete(position = "top") +
    labs(y = 'Relative Abundance (%)') +
    scale_y_reverse() +
    coord_flip() +
    theme_bw() +
    theme(panel.background = element_rect(color = 'black', size = 1,
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          axis.title.x = element_text(size = 12, family = 'Arial'),
          axis.text.x = element_text(size = 10, color = 'black', family = 'Arial'),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          axis.title.y = element_blank(),
          # legend.text = element_text(size = 10, family = 'Arial'),
          legend.position = 'none',
          # legend.key.size = unit(1, 'cm'),
          legend.key = element_blank())
ggsave(plot = plot.abundance, path = path.plot, 
       filename = paste0('bar_abundance_', 
                         paste0(as.character(vec.dose), collapse = ''), '.png'),
       height = 16, width = 8, units = 'cm')

    


