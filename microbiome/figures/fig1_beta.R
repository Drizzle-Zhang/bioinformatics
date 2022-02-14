library(ggplot2)

# data combination
type.distance <- 'bray_curtis'
use.dim <- '5'

vec.gender <- c('male', 'female')
vec.linetype <- c(1, 2)
vec.point <- c(16, 21)
vec.str.dose <- c('01', '02', '03')
vec.color <- c('firebrick', '#33A02C', '#1F78B4')

df.combine <- data.frame(stringsAsFactors = F)
for (i in 1:length(vec.gender)) {
    gender <- vec.gender[i]
    for (j in 1:length(vec.str.dose)) {
        str.dose <- vec.str.dose[j]
        path.plot <- paste0(
            '/home/drizzle_zhang/microbiome/result/5.Beta_Diversity/PCoA/', 
            type.distance, '_', gender)
        file.male <- paste0(path.plot, paste0('/diastance_', use.dim, '_', str.dose, '.txt'))
        df.sub <- read.delim(file.male, stringsAsFactors = F)
        df.sub$Gender[df.sub$Gender == 'female'] <- 'Female'
        df.sub$Gender[df.sub$Gender == 'male'] <- 'Male'
        df.sub$Dose <- rep(str.dose, nrow(df.sub))
        df.sub$LineType <- rep(vec.linetype[i], nrow(df.sub))
        df.sub$PointShape <- rep(vec.point[i], nrow(df.sub))
        df.sub$Color <- rep(vec.color[j], nrow(df.sub))
        df.combine <- rbind(df.combine, df.sub)
    }
}

series.time <- c(-1, 1,  5,  9, 17, 21, 25, 29, 33, 41, 49, 60, 68, 84)
mod <- 'sel'
df.combine <- df.combine[df.combine$Time %in% series.time,]

# baseline
linewidth.base <- 0.2
df.baseline <- df.combine[df.combine$Time == -1,]
baseline.male.01 <- df.baseline[(df.baseline$Gender == 'Male') & (df.baseline$Dose == '01'), 'Distance']
baseline.male.02 <- df.baseline[(df.baseline$Gender == 'Male') & (df.baseline$Dose == '02'), 'Distance']
baseline.male.03 <- df.baseline[(df.baseline$Gender == 'Male') & (df.baseline$Dose == '03'), 'Distance']
baseline.female.01 <- df.baseline[(df.baseline$Gender == 'Female') & (df.baseline$Dose == '01'), 'Distance']
baseline.female.02 <- df.baseline[(df.baseline$Gender == 'Female') & (df.baseline$Dose == '02'), 'Distance']
baseline.female.03 <- df.baseline[(df.baseline$Gender == 'Female') & (df.baseline$Dose == '03'), 'Distance']
df.plot <- df.combine[df.combine$Time != -1, ]

# plot
plot.fit <- 
    ggplot(data = df.plot, aes(x = Time, y = Distance, color = Dose, 
                               shape = Gender, linetype = Gender)) +
    geom_smooth(method = lm, formula = y ~ poly(x, 3), se = F) +
    geom_point() +
    scale_linetype_manual(breaks = c('Male', 'Female'), values = vec.linetype) + 
    scale_color_manual(breaks = vec.str.dose, values = vec.color,
                       labels = c('0.0 vs. 0.375', '0.0 vs. 1.5',
                                  '0.0 vs. 6.0')) + 
    labs(color = 'Dosage (mg/kg IAA)') + 
    scale_shape_manual(breaks = c('Male', 'Female'), values = vec.point) + 
    geom_hline(yintercept = baseline.male.01, color = 'firebrick', 
               linetype = 1, size = linewidth.base) +
    geom_hline(yintercept = baseline.male.02, color = '#33A02C', 
               linetype = 1, size = linewidth.base) +
    geom_hline(yintercept = baseline.male.03, color = '#1F78B4', 
               linetype = 1, size = linewidth.base) +
    geom_hline(yintercept = baseline.female.01, color = 'firebrick', 
               linetype = 5, size = linewidth.base) +
    geom_hline(yintercept = baseline.female.02, color = '#33A02C', 
               linetype = 5, size = linewidth.base) +
    geom_hline(yintercept = baseline.female.03, color = '#1F78B4', 
               linetype = 5, size = linewidth.base) +
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          panel.grid = element_blank(),
          legend.key.height = unit(20, "pt"),
          legend.key.width = unit(30, "pt")) +
    ylab('Silhouette Distance') + xlab('Time (week)')
file.fit <- '/home/drizzle_zhang/microbiome/result/Figs/Beta.png'
ggsave(filename = file.fit, plot = plot.fit, units = 'cm', height = 10, width = 22)


# wald test
library("AER")
df.plot$Gender <- as.numeric(as.factor(df.plot$Gender))
df.plot$Dose <- as.numeric(as.factor(df.plot$Dose))
fm1 <- lm(Distance ~ Time + I(Time^2) + I(Time^3), data = df.plot)
fm2 <- lm(Distance ~ Time + I(Time^2) + I(Time^3) + Gender, data = df.plot)
fm3 <- lm(Distance ~ Time + I(Time^2) + I(Time^3) + Dose, data = df.plot)
print(waldtest(fm2, fm1))
print(waldtest(fm3, fm1))

series.time <- c(1,  9, 17, 25, 33, 41, 49, 60, 68, 84)
series.time <- c(33, 41, 49, 60, 68, 84)
df.plot <- df.combine[df.combine$Time %in% series.time, ]
# df.plot <- rbind(df.plot, df.plot)
# df.plot <- rbind(df.plot, df.plot)
df.plot$Gender <- as.numeric(as.factor(df.plot$Gender))
df.plot$Dose <- as.numeric(as.factor(df.plot$Dose))
fm1 <- lm(Distance ~ Time + I(Time^2) + I(Time^3), data = df.plot)
# 0.025
fm2 <- lm(Distance ~ Time + I(Time^2) + I(Time^3) + Gender, data = df.plot)
# 0.632

fm3 <- lm(Distance ~ Time + I(Time^2) + I(Time^3) + Dose, data = df.plot)
print(waldtest(fm2, fm1))
print(waldtest(fm3, fm1))
wilcox.test(df.plot$Distance[df.plot$Gender == 'male'],
            df.plot$Distance[df.plot$Gender == 'female'], paired = T)



