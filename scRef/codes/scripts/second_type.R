source('/home/zy/my_git/scRef/main/scRef.v20.R')

############# regard sc-counts data as reference
path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_merged.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('Unclassified'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic <- OUT$mat_exp
label_Tasic <- OUT$label
ref.labels <-label_Tasic[, 1]
ref.mtx <- exp_Tasic
ref.dataset <- 'Tasic'

############### import unlabeled data
############### Campbell
path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Campbell'
file.data.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat.txt')
file.label.unlabeled <- paste0(path.input, dataset, '_exp_sc_mat_cluster_original.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Habib <- OUT$mat_exp
label_Habib <- OUT$label
exp_sc_mat <- exp_Habib
label_sc <- label_Habib

source('/home/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels, cluster.cell = 5, CPU = 8)
print(result.scref$run.time)
file.res <- paste0(path.output, 'scMAGIC_', ref.dataset, '_', dataset, '.Rdata')
saveRDS(result.scref, file = file.res)
result.scref <- readRDS(file.res)
info.cluster <- result.scref$info.cluster
info.cluster$percent.othertype <- info.cluster$percent.other - info.cluster$percent.Unassigned
info.cluster$num.type <- info.cluster$num.cell * info.cluster$percent.othertype
dim(exp_sc_mat)[2]
sum(info.cluster$num.type[info.cluster$percent.Unassigned < 0.2 & info.cluster$percent.othertype > 0])

############### Tasic2018
path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Tasic2018'
file.data.unlabeled <- paste0(path.input, dataset, '/cell_exp.txt')
file.label.unlabeled <- paste0(path.input, dataset, '/cell_meta.txt')
# OUT <- prepare.data(file.data.unlabeled, file.label.unlabeled, del.label = c('miss'))
# saveRDS(OUT, file = paste0(path.output, dataset, '.Rdata'))
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Tasic2018 <- OUT$mat_exp
label_Tasic2018 <- OUT$label
exp_sc_mat <- exp_Tasic2018
label_sc <- label_Tasic2018

source('/home/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = T, cluster.cell = 10,
                      min_cell = 10, CPU = 8)
file.res <- paste0(path.output, 'scMAGIC_', ref.dataset, '_', dataset, '.Rdata')
saveRDS(result.scref, file = file.res)
result.scref <- readRDS(file.res)
info.cluster <- result.scref$info.cluster
info.cluster$percent.othertype <- info.cluster$percent.other - info.cluster$percent.Unassigned
info.cluster$num.type <- info.cluster$num.cell * info.cluster$percent.othertype
dim(exp_sc_mat)[2]
sum(info.cluster$num.type[info.cluster$percent.Unassigned < 0.2 & info.cluster$percent.othertype > 0])


############### HochgernerA
file.in <- '/home/zy/scRef/sc_data/Hochgerner/GSE95315_10X_expression_data_v2.tab'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'HochgernerA'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Hochgerner <- OUT$mat_exp
label_Hochgerner <- OUT$label
exp_sc_mat <- exp_Hochgerner
label_sc <- label_Hochgerner

source('/home/zy/my_git/scRef/main/scRef.v21.R')
setwd('~/my_git/scRef')
result.scref <- SCREF(exp_sc_mat, ref.mtx, ref.labels,
                      type_ref = 'sc-counts', use.RUVseq = T, 
                      cluster.speed = T, cluster.cell = 5,
                      min_cell = 10, CPU = 8)
file.res <- paste0(path.output, 'scMAGIC_', ref.dataset, '_', dataset, '.Rdata')
saveRDS(result.scref, file = file.res)
result.scref <- readRDS(file.res)
info.cluster <- result.scref$info.cluster
info.cluster$percent.othertype <- info.cluster$percent.other - info.cluster$percent.Unassigned
info.cluster$num.type <- info.cluster$num.cell * info.cluster$percent.othertype
dim(exp_sc_mat)[2]
sum(info.cluster$num.type[info.cluster$percent.Unassigned < 0.2 & info.cluster$percent.othertype > 0])

############### Mizrak
path.input <- '/home/zy/scRef/sc_data/'
path.output <- '/home/zy/scRef/Benchmark/mouse_brain/'
dataset <- 'Mizrak'
OUT <- readRDS(paste0(path.output, dataset, '.Rdata'))
exp_Mizrak <- OUT$mat_exp
label_Mizrak <- OUT$label
exp_sc_mat <- exp_Mizrak
label_sc <- label_Mizrak
