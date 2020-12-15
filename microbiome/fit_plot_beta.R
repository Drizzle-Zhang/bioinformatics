library(ggplot2)

# gender
type.distance <- 'bray_curtis'
str.dose <- '01'
use.dim <- '5'

gender <- 'male'
path.plot <- paste0(
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
    type.distance, '_', gender)
file.male <- paste0(path.plot, paste0('/diastance_', use.dim, '_', str.dose, '.txt'))
df.male <- read.delim(file.male)

gender <- 'female'
path.plot <- paste0(
    '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
    type.distance, '_', gender)
file.female <- paste0(path.plot, paste0('/diastance_', use.dim, '_', str.dose, '.txt'))
df.female <- read.delim(file.female)

df.rbind <- rbind(df.male, df.female)

# series.time <- unique(df.rbind$Time)
# mod <- 'all'

series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
mod <- 'sel'
df.rbind <- df.rbind[df.rbind$Time %in% series.time,]

baseline.male <- df.rbind[(df.rbind$Gender == 'male') & (df.rbind$Time == -1), 'Distance']
baseline.female <- df.rbind[(df.rbind$Gender == 'female') & (df.rbind$Time == -1), 'Distance']
df.plot <- df.rbind[df.rbind$Time != -1, ]
df.plot$Color <- rep('0', dim(df.plot)[1])
df.plot[df.plot$Gender == 'male', 'Color'] <- 
    rep('#F8766D', dim(df.plot[df.plot$Gender == 'male', ])[1])
df.plot[df.plot$Gender == 'female', 'Color'] <- 
    rep('#00BFC4', dim(df.plot[df.plot$Gender == 'female', ])[1])

plot.fit <- 
    ggplot(data = df.plot, aes(x = Time, y = Distance, color = Gender)) +
    geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
    geom_point(color = df.plot$Color) +
    geom_hline(yintercept = baseline.male, color = '#F8766D', linetype = 5) +
    geom_hline(yintercept = baseline.female, color = '#00BFC4', linetype = 5) +
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent')) +
    ylab('Silhouette Distance')
file.fit <- 
    paste0('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
           'distance_compare_', mod, '_',
           type.distance, '_', str.dose, '_', use.dim, '.png')
ggsave(filename = file.fit, plot = plot.fit, units = 'cm', height = 10, width = 20)

# dose
type.distance <- 'bray_curtis'
str.dose <- c('01', '02', '03')
use.dim <- '5'

gender <- 'male'
df.male <- data.frame()
for (sub.dose in str.dose) {
    path.plot <- paste0(
        '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
        type.distance, '_', gender)
    file.male <- paste0(path.plot, paste0('/diastance_', use.dim, '_', sub.dose, '.txt'))
    df.sub <- read.delim(file.male)
    df.sub$Dose <- rep(substr(sub.dose, 2, 2), dim(df.sub)[1])
    df.male <- rbind(df.male, df.sub)
}

df.male <- df.male[df.male$Time %in% series.time,]

baseline.male.1 <- df.male[(df.male$Dose == '1') & (df.male$Time == -1), 'Distance']
baseline.male.2 <- df.male[(df.male$Dose == '2') & (df.male$Time == -1), 'Distance']
baseline.male.3 <- df.male[(df.male$Dose == '3') & (df.male$Time == -1), 'Distance']
df.plot.male <- df.male[df.male$Time != -1, ]
df.plot.male$Color <- rep('0', dim(df.plot.male)[1])
df.plot.male[df.plot.male$Dose == '1', 'Color'] <- 
    rep('#F8766D', dim(df.plot.male[df.plot.male$Dose == '1', ])[1])
df.plot.male[df.plot.male$Dose == '2', 'Color'] <- 
    rep('#00BA38', dim(df.plot.male[df.plot.male$Dose == '2', ])[1])
df.plot.male[df.plot.male$Dose == '3', 'Color'] <- 
    rep('#619CFF', dim(df.plot.male[df.plot.male$Dose == '3', ])[1])

plot.fit <- 
    ggplot(data = df.plot.male, aes(x = Time, y = Distance, color = Dose)) +
    geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
    geom_point(color = df.plot.male$Color) +
    geom_hline(yintercept = baseline.male.1, color = '#F8766D', linetype = 5) +
    geom_hline(yintercept = baseline.male.2, color = '#00BA38', linetype = 5) +
    geom_hline(yintercept = baseline.male.3, color = '#619CFF', linetype = 5) +
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent')) +
    ylab('Silhouette Distance')
file.fit <- 
    paste0('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
           'distance_compare_dose_', mod, '_', 
           type.distance, '_', gender, '_', use.dim, '.png')
ggsave(filename = file.fit, plot = plot.fit, units = 'cm', height = 10, width = 20)

# # normalization
# df.norm.male <- df.plot.male
# df.norm.male[df.norm.male$Dose == '1', 'Distance'] <- 
#     df.plot.male[df.plot.male$Dose == '1', 'Distance'] - baseline.male.1
# df.norm.male[df.norm.male$Dose == '2', 'Distance'] <- 
#     df.plot.male[df.plot.male$Dose == '2', 'Distance'] - baseline.male.2
# df.norm.male[df.norm.male$Dose == '3', 'Distance'] <- 
#     df.plot.male[df.plot.male$Dose == '3', 'Distance'] - baseline.male.3
# 
# ggplot(data = df.norm.male, aes(x = Time, y = Distance, color = Dose)) +
#     geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
#     geom_point(color = df.plot.male$Color) +
#     theme(panel.background = element_rect(color = 'gray',
#                                           fill = 'transparent')) +
#     ylab('Silhouette Distance')


gender <- 'female'
df.female <- data.frame()
for (sub.dose in str.dose) {
    path.plot <- paste0(
        '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
        type.distance, '_', gender)
    file.female <- paste0(path.plot, paste0('/diastance_', use.dim, '_', sub.dose, '.txt'))
    df.sub <- read.delim(file.female)
    df.sub$Dose <- rep(substr(sub.dose, 2, 2), dim(df.sub)[1])
    df.female <- rbind(df.female, df.sub)
}

df.female <- df.female[df.female$Time %in% series.time,]

baseline.female.1 <- df.female[(df.female$Dose == '1') & (df.female$Time == -1), 'Distance']
baseline.female.2 <- df.female[(df.female$Dose == '2') & (df.female$Time == -1), 'Distance']
baseline.female.3 <- df.female[(df.female$Dose == '3') & (df.female$Time == -1), 'Distance']
df.plot.female <- df.female[df.female$Time != -1, ]
df.plot.female$Color <- rep('0', dim(df.plot.female)[1])
df.plot.female[df.plot.female$Dose == '1', 'Color'] <- 
    rep('#F8766D', dim(df.plot.female[df.plot.female$Dose == '1', ])[1])
df.plot.female[df.plot.female$Dose == '2', 'Color'] <- 
    rep('#00BA38', dim(df.plot.female[df.plot.female$Dose == '2', ])[1])
df.plot.female[df.plot.female$Dose == '3', 'Color'] <- 
    rep('#619CFF', dim(df.plot.female[df.plot.female$Dose == '3', ])[1])

plot.fit <- 
    ggplot(data = df.plot.female, aes(x = Time, y = Distance, color = Dose)) +
    geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
    geom_point(color = df.plot.female$Color) +
    geom_hline(yintercept = baseline.female.1, color = '#F8766D', linetype = 5) +
    geom_hline(yintercept = baseline.female.2, color = '#00BA38', linetype = 5) +
    geom_hline(yintercept = baseline.female.3, color = '#619CFF', linetype = 5) +
    theme(panel.background = element_rect(color = 'gray',
                                          fill = 'transparent')) +
    ylab('Silhouette Distance')
file.fit <- 
    paste0('/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
           'distance_compare_dose_', mod, '_'
           , type.distance, '_', gender, '_', use.dim, '.png')
ggsave(filename = file.fit, plot = plot.fit, units = 'cm', height = 10, width = 20)

# # normalization
# df.norm.female <- df.plot.female
# df.norm.female[df.norm.female$Dose == '1', 'Distance'] <- 
#     df.plot.female[df.plot.female$Dose == '1', 'Distance'] - baseline.female.1
# df.norm.female[df.norm.female$Dose == '2', 'Distance'] <- 
#     df.plot.female[df.plot.female$Dose == '2', 'Distance'] - baseline.female.2
# df.norm.female[df.norm.female$Dose == '3', 'Distance'] <- 
#     df.plot.female[df.plot.female$Dose == '3', 'Distance'] - baseline.female.3
# 
# ggplot(data = df.norm.female, aes(x = Time, y = Distance, color = Dose)) +
#     geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
#     geom_point(color = df.plot.female$Color) +
#     theme(panel.background = element_rect(color = 'gray',
#                                           fill = 'transparent')) +
#     ylab('Silhouette Distance')



