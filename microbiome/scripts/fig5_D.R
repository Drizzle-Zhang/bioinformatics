library(RColorBrewer)
library(ggplot2)
library(tidyverse)
library(patchwork)
library(scales)

# meta file
meta.file <- '/home/drizzle_zhang/microbiome/result/meta_sample.txt'
df.meta <- read.delim(meta.file, stringsAsFactors = FALSE)


topN <- 50
path.tree <- '/home/drizzle_zhang/microbiome/result/6.Phylogenetic/graphlan'

level = 'genus'
file = paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/abundance/relative_abundance/', 
              level, '.csv')
genus <- read.table(file = file, header = TRUE, row.names = 1,
                    sep = '\t', stringsAsFactors = FALSE,
                    check.names = FALSE)

path.input <- paste0('/home/drizzle_zhang/microbiome/result/3.Community_Structure/barplot_',
                    level)
path.plot <- '/home/drizzle_zhang/microbiome/result/Figs/'

sel.genus <- c('Lactobacillus', 'Clostridium_sensu_stricto_1',
               'Helicobacter', 'Turicibacter', 'Streptococcus', 
               'Ruminococcaceae_UCG-014',
               'Prevotella_1', 'Prevotella_9', 'Treponema_2')

# mean abundance
merge.abundance <- genus[sel.genus, ]
mean.abundance <- rowMeans(merge.abundance)
sort(mean.abundance)
sort.abundance <- c('Lactobacillus', 'Ruminococcaceae_UCG-014',
                    'Prevotella_9', 'Helicobacter', 'Treponema_2',
                    'Turicibacter', 'Clostridium_sensu_stricto_1',
                    'Prevotella_1', 'Streptococcus')

list.vec.dose <- list(c(0, 1, 2, 3), c(0, 1), c(0, 2), c(0, 3))
vec.gender <- c('male', 'female')
series.time <- unique(df.meta$Time)
for (vec.dose in list.vec.dose) {
    # dose
    str_dose <- paste0(as.character(vec.dose), collapse = '')
    path.dose <- paste(path.tree, str_dose, sep = '/')

    df.combine <- data.frame(stringsAsFactors = F)
    for (gender in vec.gender) {
        # gender
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
            sub.abundance <- sub.genus
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
                p.value <- out.wilcox$p.value
                p.value.round <- round(out.wilcox$p.value, 4)
                df.out <- rbind(df.out, 
                                data.frame(type = type, fold.change = fold.change, 
                                           p.value = p.value, 
                                           time = sub.time, p.value.round = p.value.round))
            }
            df.plot.point <- rbind(df.plot.point, df.out)
        }
        df.plot.point$log.pval <- -log10(df.plot.point$p.value)
        df.plot.point$Gender <- rep(gender, nrow(df.plot.point))
        df.combine <- rbind(df.combine, df.plot.point)
    }
    df.combine$logPval <- df.combine$log.pval
    df.combine$logPval[df.combine$p.value > 0.31] <- NA
    df.combine$fold.change[df.combine$fold.change > 0.95] <- 0.95
    df.combine$fold.change[df.combine$fold.change < -0.95] <- -0.95
    df.combine$direct.logp <- df.combine$fold.change*df.combine$logPval
    df.combine$direct.logp[df.combine$direct.logp > 0] <- 1
    df.combine$direct.logp[df.combine$direct.logp < 0] <- -1
    df.final <- df.combine
    df.final$Taxonomy <- factor(df.final$type, levels = sort.abundance)
    plot.bubble <- 
        ggplot(data = df.final, 
               aes(x = time, y = Gender, size = logPval, color = fold.change)) + 
        geom_point(fill = 'cornsilk') + 
        scale_colour_gradient2(low = '#00008B', mid = "white", high = '#8B0000', na.value = 'white') +
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
    sub.genders <- rep(c("Male","Female"), times = nrow(df.final))
    sub.genus <- rep(rev(rownames(df.final)), c(rep(2, nrow(df.final))))
    sort.row <- paste(sub.genders, sub.genus, sep = '_')
    df.status <- data.frame(status=sub.genders,
                            p = rep('white', length(sort.row)), 
                            row_name = sort.row, 
                            cluster = sub.genus)
    df.status$row_name <- factor(df.status$row_name, levels = rev(sort.row))
    df.status$cluster <- factor(df.status$cluster, levels = rownames(df.final))
    df.status$status <- factor(df.status$status, levels = c("Male","Female"))
    plot.Status <- 
        ggplot(data = df.status, aes(x=p,y=row_name,fill=status))+
        geom_tile() + 
        facet_grid(cluster ~ ., scales = 'free', space = 'free') +
        scale_y_discrete(position="right") +
        scale_fill_manual(breaks = c("Male","Female"), 
                          values = c("#BC80BD", "#5AB4AC")) +
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
           height = 12, width = 28, units = 'cm')
    
    
    # abundance
    df.abundance <- data.frame(stringsAsFactors = F)
    for (gender in vec.gender) {
        # gender
        df.meta.gender <- df.meta[df.meta$Gender == gender & df.meta$Dose %in% vec.dose, ]
        use.sample <- df.meta.gender$SampleName
        merge.abundance <- genus[sel.genus, use.sample]
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
        labs(y = 'Average Relative\nAbundance (%)') +
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
              legend.position = 'none',
              legend.key = element_blank())
    ggsave(plot = plot.abundance, path = path.plot, 
           filename = paste0('bar_abundance_', 
                             paste0(as.character(vec.dose), collapse = ''), '.png'),
           height = 12, width = 5, units = 'cm')
    
}

