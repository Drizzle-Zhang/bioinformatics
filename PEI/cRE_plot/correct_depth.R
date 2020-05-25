library(ggplot2)
library(VIM)
library(car)
library(lattice)
setwd('/home/drizzle_zhang/driver_mutation/cRE_plot/model_test')

file.in = './score.txt'
df.scores <- read.delim(file.in, sep = '\t', header = T,
                        row.names = 'peak_id', stringsAsFactors = F)

# original scatter plot
df.2 <- df.scores[,c('ENCFF273MVV', 'ENCFF804BNU')]
df.2 <- na.omit(df.2)
fit <- lm(ENCFF273MVV ~ ENCFF804BNU, data = df.2)
summary(fit)

ggplot(data = df.2, aes(x = ENCFF804BNU, y = ENCFF273MVV)) + 
    geom_smooth(method = lm, formula = y~(x)) + geom_point(size = 0.1)

# correct (normalization)
densityplot(~(ENCFF273MVV), data = df.scores)
densityplot(~(ENCFF273MVV)^-0.1681, data = df.scores)
summary(powerTransform(df.scores$ENCFF804BNU))
densityplot(~(ENCFF804BNU)^-0.3465, data = df.scores)

df.2$ENCFF273MVV_norm <- (df.2$ENCFF273MVV)^-0.1681
df.2$ENCFF804BNU_norm <- (df.2$ENCFF804BNU)^-0.3465
ggplot(data = df.2, aes(x = ENCFF804BNU_norm, y = ENCFF273MVV_norm)) + 
    geom_smooth(method = glm, formula = y~(x)) + geom_point(size = 0.1)
cor(df.2$ENCFF273MVV, df.2$ENCFF804BNU)

# optimize fitting result
boxTidwell(ENCFF273MVV_norm ~ ENCFF804BNU_norm, data = df.2)
df.2$ENCFF804BNU_boxtid <- (df.2$ENCFF804BNU_norm)^0.23244
fit <- lm(ENCFF273MVV_norm ~ ENCFF804BNU_boxtid, data = df.2)
summary(fit)
ggplot(data = df.2, aes(x = ENCFF804BNU_boxtid, y = ENCFF273MVV_norm)) + 
    geom_smooth(method = lm, formula = y~(x)) + geom_point(size = 0.1)
df.2$ENCFF804BNU_pred <- 0.908264*df.2$ENCFF804BNU_boxtid - 0.070924
cor(df.2$ENCFF273MVV_norm, df.2$ENCFF804BNU_pred)

outlierTest(fit)
influencePlot(fit)
df.2$hatvalue <- hatvalues(fit)
df.2$rstudent <- rstudent(fit)
cutoff.hatvalue <- 3*mean(df.2$hatvalue)
df.2.del <- df.2[(df.2$hatvalue < cutoff.hatvalue) & 
                 (df.2$rstudent < 2) & (df.2$rstudent > -2),]

ggplot(data = df.2.del, aes(x = ENCFF804BNU_boxtid, y = ENCFF273MVV_norm)) + 
    geom_smooth(method = lm, formula = y~(x)) + geom_point(size = 0.1)

fit <- lm(ENCFF273MVV_norm ~ ENCFF804BNU_boxtid, data = df.2.del)
summary(fit)

df.2.del$ENCFF804BNU_pred <- 0.966241*df.2.del$ENCFF804BNU_boxtid - 0.110385
cor(df.2.del$ENCFF273MVV_norm, df.2.del$ENCFF804BNU_pred)
ggplot() + geom_point(aes(x = ENCFF804BNU_boxtid, y = ENCFF273MVV_norm), 
                      data = df.2, size = 0.1) + 
    geom_smooth(aes(x = ENCFF804BNU_boxtid, y = ENCFF273MVV_norm), 
                data = df.2.del, method = lm, formula = y~(x), fullrange = T)

df.2$ENCFF804BNU_pred <- 0.966241*df.2$ENCFF804BNU_boxtid - 0.110385
df.2.res <- df.2[,c('ENCFF273MVV', 'ENCFF804BNU', 
                    'ENCFF273MVV_norm', 'ENCFF804BNU_pred')]



df.1 <- df.scores[,c('ENCFF273MVV', 'ENCFF097LEF')]
df.1 <- na.omit(df.1)

# correct (normalization)
summary(powerTransform(df.scores$ENCFF273MVV))
densityplot(~(ENCFF273MVV)^-0.1681, data = df.scores)
summary(powerTransform(df.scores$ENCFF097LEF))
densityplot(~(ENCFF097LEF)^-0.0977, data = df.scores)

df.1$ENCFF273MVV_norm <- (df.1$ENCFF273MVV)^-0.1681
df.1$ENCFF097LEF_norm <- (df.1$ENCFF097LEF)^-0.0977

boxTidwell(ENCFF273MVV_norm ~ ENCFF097LEF_norm, data = df.1)
df.1$ENCFF097LEF_boxtid <- (df.1$ENCFF097LEF_norm)^-0.96316
fit <- lm(ENCFF273MVV_norm ~ ENCFF097LEF_boxtid, data = df.1)
df.1$hatvalue <- hatvalues(fit)
df.1$rstudent <- rstudent(fit)
cutoff.hatvalue <- 3*mean(df.1$hatvalue)
df.1.del <- df.1[(df.1$hatvalue < cutoff.hatvalue) & 
                     (df.1$rstudent < 2) & (df.1$rstudent > -2),]

ggplot(data = df.1.del, aes(x = ENCFF097LEF_boxtid, y = ENCFF273MVV_norm)) + 
    geom_smooth(method = lm, formula = y~(x)) + geom_point(size = 0.1)

fit <- lm(ENCFF273MVV_norm ~ ENCFF097LEF_boxtid, data = df.1.del)
summary(fit)

df.1.del$ENCFF097LEF_pred <- -0.490383*df.1.del$ENCFF097LEF_boxtid + 1.200160
cor(df.1.del$ENCFF273MVV_norm, df.1.del$ENCFF097LEF_pred)
ggplot() + geom_point(aes(x = ENCFF097LEF_boxtid, y = ENCFF273MVV_norm), 
                      data = df.1, size = 0.1) + 
    geom_smooth(aes(x = ENCFF097LEF_boxtid, y = ENCFF273MVV_norm), 
                data = df.1.del, method = lm, formula = y~(x), fullrange = T)

