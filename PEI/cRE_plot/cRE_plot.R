setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(ggplot2)

# cRE type
df.count.organ <- read.delim('organ_count.txt', 
                             sep = '\t', stringsAsFactors = F)

ggplot(df.count.organ, aes(x = Organ, y = Count, fill = Type)) + 
    geom_bar(stat = "identity", position = "stack") + coord_flip()


# boxplot function
boxplot <- function(path, file_in, fig_out, x_name, y_name,
                    width = 30, height = 25) {
    df.count <- read.delim(paste(path, file_in, sep = '/'), 
                           sep = '\t', stringsAsFactors = F)
    df.boxplot <- data.frame()
        for (i in 1:dim(df.count)[1]) {
            rows <- df.count[i, y_name]
            organs <- strsplit(df.count[i, x_name], ',')[[1]]
            for (organ in organs) {
                df.boxplot <- rbind(df.boxplot,
                                    data.frame(rows = rows, organ = organ))
            }
            df.boxplot <- rbind(df.boxplot,
                                data.frame(rows = rows, organ = 'All'))
        }
    plot.out <- ggplot(df.boxplot, aes(x = organ, y = rows, fill = organ)) + 
        geom_boxplot() + coord_flip()
    ggsave(filename = fig_out, path = path, plot = plot.out, units = 'cm',
           width = 30, height = 25)
}

# DHS
# peaks number in different organ
boxplot('./DHS', 'count.txt', 'count.organ.png', 
        'Biosample.organ', 'Biosample.file.rows')
# peaks number in different life stage
boxplot('./DHS', 'count.txt', 'count.lifestage.png', 
        'Biosample.life.stage', 'Biosample.file.rows', width = 20, height = 5)

# H3K27ac
boxplot('./H3K27ac', 'count.txt', 'count.organ.png', 
        'Biosample.organ', 'Biosample.file.rows')
boxplot('./H3K27ac', 'count.txt', 'count.lifestage.png', 
        'Biosample.life.stage', 'Biosample.file.rows', width = 20, height = 5)
boxplot('./H3K27ac', 'distribution.txt', 'medium.organ.png', 
        'Biosample.organ', 'X50', width = 20, height = 5)


# H3K4me3
boxplot('./H3K4me3', 'count.txt', 'count.organ.png', 
        'Biosample.organ', 'Biosample.file.rows')
boxplot('./H3K4me3', 'count.txt', 'count.lifestage.png', 
        'Biosample.life.stage', 'Biosample.file.rows', width = 20, height = 5)




