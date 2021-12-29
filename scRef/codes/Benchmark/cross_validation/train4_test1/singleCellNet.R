setwd('/home/zy/scRef/try_data')

# download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_Park_MouseKidney_062118.rda", "sampTab_Park_MouseKidney_062118.rda")
# 
# download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expMatrix_Park_MouseKidney_Oct_12_2018.rda", "expMatrix_Park_MouseKidney_Oct_12_2018.rda")

# download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/expMatrix_TM_Raw_Oct_12_2018.rda", "expMatrix_TM_Raw_Oct_12_2018.rda")
# 
# download.file("https://s3.amazonaws.com/cnobjects/singleCellNet/examples/sampTab_TM_053018.rda", "sampTab_TM_053018.rda")

stPark = utils_loadObject("sampTab_Park_MouseKidney_062118.rda")
expPark = utils_loadObject("expMatrix_Park_MouseKidney_Oct_12_2018.rda")
dim(expPark)

expTMraw = utils_loadObject("expMatrix_TM_Raw_Oct_12_2018.rda")
dim(expTMraw)

stTM = utils_loadObject("sampTab_TM_053018.rda")
dim(stTM)

stTM<-droplevels(stTM)

commonGenes = intersect(rownames(expTMraw), genesPark)
length(commonGenes)

expTMraw = expTMraw[commonGenes,]

set.seed(100) #can be any random seed number
stList = splitCommon(sampTab=stTM, ncells=100, dLevel="newAnn")
stTrain = stList[[1]]
expTrain = as.matrix(expTMraw[,rownames(stTrain)])

system.time(class_info<-scn_train(stTrain = stTrain, expTrain = expTrain, dLevel = "newAnn", colName_samp = "cell"))
