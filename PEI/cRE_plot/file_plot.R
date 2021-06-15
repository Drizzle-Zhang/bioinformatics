setwd('/home/drizzle_zhang/driver_mutation/cRE_plot')
library(VennDiagram)

######### VeenDiagram
extend_df <- function(df_input) {
    df_extend <- data.frame()
    for (i in 1:dim(df_input)[1]) {
        organs <- strsplit(df_input$Biosample.organ[i], ',')
        for (organ in organs) {
            df_extend <- rbind(
                df_extend,
                data.frame(
                    Biosample.term.name = df_input$Biosample.term.name[i],
                    Biosample.life.stage = df_input$Biosample.life.stage[i],
                    Biosample.organ = organ,
                    organ_life.stage = paste(
                        organ, df_input$Biosample.life.stage[i], sep = '|')))
        }
    }
    return(df_extend)
}

df_dhs <- read.delim('metadata_DHS.tsv', sep = '\t', stringsAsFactors = F)
df_H3K27ac <- read.delim('metadata_H3K27ac.tsv', 
                         sep = '\t', stringsAsFactors = F)
df_H3K4me3 <- read.delim('metadata_H3K4me3.tsv', 
                         sep = '\t', stringsAsFactors = F)

df_dhs_extend <- extend_df(df_dhs)
df_H3K27ac_extend <- extend_df(df_H3K27ac)
df_H3K4me3_extend <- extend_df(df_H3K4me3)

# organ
organ_dhs <- unique(df_dhs_extend$Biosample.organ)
organ_H3K27ac <- unique(df_H3K27ac_extend$Biosample.organ)
organ_H3K4me3 <- unique(df_H3K4me3_extend$Biosample.organ)
venn.diagram(
    list(DHS = organ_dhs, H3K27ac = organ_H3K27ac, H3K4me3 = organ_H3K4me3),
    fill = c("red", "green", "blue"),
    alpha = c(0.5, 0.5, 0.5),
    cex = 2,
    cat.cex = 2,
    cat.dist = 0.1,
    cat.fontface = 4,
    fontfamily = 3,
    filename = 'Venn_organ.tiff'
)

# life stage
life_dhs <- unique(df_dhs_extend$organ_life.stage)
life_H3K27ac <- unique(df_H3K27ac_extend$organ_life.stage)
life_H3K4me3 <- unique(df_H3K4me3_extend$organ_life.stage)
venn.diagram(
    list(DHS = life_dhs, H3K27ac = life_H3K27ac, H3K4me3 = life_H3K4me3),
    fill = c("red", "green", "blue"),
    alpha = c(0.5, 0.5, 0.5),
    cex = 2,
    cat.cex = 2,
    cat.dist = 0,
    cat.fontface = 4,
    fontfamily = 3,
    filename = 'Venn_organ_life.tiff'
)

# term 
term_dhs <- df_dhs$Biosample.term.name
term_H3K27ac <- df_H3K27ac$Biosample.term.name
term_H3K4me3 <- df_H3K4me3$Biosample.term.name
venn.diagram(
    list(DHS = term_dhs, H3K27ac = term_H3K27ac, H3K4me3 = term_H3K4me3),
    fill = c("red", "green", "blue"),
    alpha = c(0.5, 0.5, 0.5),
    cex = 2,
    cat.cex = 2,
    cat.dist = 0,
    cat.fontface = 4,
    fontfamily = 3,
    filename = 'Venn_term.tiff'
)

