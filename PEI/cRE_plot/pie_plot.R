library(ggplot2)
dt = data.frame(A = c(1559, 2800, 1527), 
                B = c('Enhancer','Protein-Promoter','Other-Promoter'))
dt = data.frame(A = c(34802, 90660, 45211), 
                B = c('Enhancer','Protein-Promoter','Other-Promoter'))
dt = data.frame(A = c(26804, 29549, 26863), 
                B = c('Enhancer','Protein-Promoter','Other-Promoter'))
dt = data.frame(A = c(26806, 33027, 27934, 10268, 29663), 
                B = c('Enhancer', 'Insulator', 'Other-Promoter(Enhancer)',
                      'Protein-Promoter','Protein-Promoter(Enhancer)'))

# dt = dt[order(dt$A, decreasing = TRUE),]
myLabel = as.vector(dt$B)   
myLabel = paste(myLabel, "(", round(dt$A / sum(dt$A) * 100, 2), "%)", sep = "")   

ggplot(dt, aes(x = "", y = A, fill = B)) +
    geom_bar(stat = "identity", width = 1) +    
    coord_polar(theta = "y", direction = 1) + 
    labs(x = "", y = "", title = "") + 
    theme(axis.ticks = element_blank()) + 
    theme(legend.title = element_blank(), legend.position = "top") + 
    scale_fill_discrete(breaks = dt$B, labels = myLabel) + 
    theme(axis.text.x = element_blank()) + 
    geom_text(aes(y = (sum(A) - (A/2 + c(0, cumsum(A)[-length(A)]))), label = myLabel), size = 5) +
    theme(panel.grid=element_blank()) +    ## 去掉白色圆框和中间的坐标线
    theme(panel.border=element_blank()) +    ## 去掉最外层正方形的框框
    theme(panel.background=element_blank())   ## 去掉背景颜色
