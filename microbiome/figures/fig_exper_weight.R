library(ggplot2)
library(Rmisc)
library(cowplot)
library(stringr)

df_weight <- read.delim2('/home/drizzle_zhang/microbiome/experiment_data/Weight.txt')
df_weight <- df_weight[, c(1, 2, 3)]
df_weight$Dosage.mg.kg.IAA. <- as.character(df_weight$Dosage.mg.kg.IAA.)
df_weight$Gender <- unlist(lapply(df_weight$Dosage.mg.kg.IAA., 
                                  function(x) {strsplit(x, split = ') ', fixed = T)[[1]][1]}))
df_weight$Dosage <- unlist(lapply(df_weight$Dosage.mg.kg.IAA., 
                                  function(x) {strsplit(x, split = ') ', fixed = T)[[1]][2]}))
df_weight <- na.omit(df_weight[, c("Time.weeks.", "Weight", "Gender", "Dosage")])
tgc1 <- summarySE(df_weight, measurevar = "Weight", 
                  groupvars = c("Gender", "Dosage", "Time.weeks."))
tgc1$Gender <- str_replace_all(tgc1$Gender, fixed('(Female'), 'Female')
tgc1$Gender <- str_replace_all(tgc1$Gender, fixed('(Male'), 'Male')

plot.weight <- 
    ggplot(tgc1, aes(x = Time.weeks., y = Weight, color = Dosage)) + 
    geom_errorbar(aes(ymin = Weight - se, ymax = Weight + se), width = 1.5) +
    geom_line() +
    geom_point(size = 0.5) +
    facet_grid(Gender~., scales = 'free') +
    scale_color_manual(breaks = c('0.0', '0.375', '1.5', '6.0'),
                         values = c('dimgrey', 'firebrick', '#33A02C', '#1F78B4')) +
    labs(x = 'Time (week)', y = 'Weight (g)', color = 'Dosage (mg/kg IAA)') +
    theme_bw() + 
    theme(panel.background = element_rect(color = 'black',
                                          fill = 'transparent'),
          strip.background = element_rect(
              color = 'black', fill = 'transparent'),
          axis.text = element_text(colour = 'black'),
          legend.title = element_text(size = 9, family = 'Arial'),
          legend.text = element_text(size = 8, family = 'Arial'),
          panel.grid = element_blank())
file.weight <- '/home/drizzle_zhang/microbiome/result/Figs/Weight.png'
ggsave(filename = file.weight, plot = plot.weight, units = 'cm', height = 10, width = 14)

# female
df_female <- df_weight[df_weight$Gender == '(Female', ]
df_female$Time.weeks. <- factor(df_female$Time.weeks.)
df_female$Dosage <- factor(df_female$Dosage)

fm1 <- lm(Weight ~ Time.weeks., data = df_female)
fm2 <- lm(Weight ~ Time.weeks. + Dosage, data = df_female)
waldtest(fm2, fm1)

# male
df_male <- df_weight[df_weight$Gender == '(Male', ]
df_male$Time.weeks. <- factor(df_male$Time.weeks.)
df_male$Dosage <- factor(df_male$Dosage)

fm1 <- lm(Weight ~ Time.weeks., data = df_male)
fm2 <- lm(Weight ~ Time.weeks. + Dosage, data = df_male)
waldtest(fm2, fm1)
