#Source for the packages to be downloaded from
source("http://bioconductor.org/biocLite.R")

#Downloading necessary packages
biocLite('caret')
biocLite('ggplot2')
biocLite('kernlab')

#Loading the packages
library(caret)
library(ggplot2)
library(kernlab)

#The working direcory, where your files will be saved
setwd("/dir/to/workdir")
#loads the data into R
train_data<-read.table("train_table.txt",sep=" ",dec=".", na.strings = "NaN", header=T, stringsAsFactors = T)

#Loads the contig classes table, and put it as class in train_data
type<-read.table("contig.classes.txt",sep="\t",dec=".", na.strings = "NaN", header=F, stringsAsFactors = T)
train_data$Class<-type$V2
#sets the rownames to contigname
rownames(train_data)<-train_data$Name
train_data<-train_data[,2:9]

#Downsamples the data so there are equal amounts of the different classes
InTrain <- downSample(train_data[,1:7],train_data[,8], yname="Class", list=F)

#split the downsample data in 1/3 to be able to run on local computer
datapart <- createDataPartition(InTrain$Class, p=1/3, list=F)
trainSmall<-InTrain[datapart,]
testSmall<-InTrain[-datapart,]
#traincontrol with method repeted crossvalidation with 3 k-fold cross validation
#classprobs=True to give the statistics for each resample
trCtrl <- trainControl(method = "repeatedcv", repeats = 3, classProbs=T)

#set seed is to assure reproducibility (set.seed(1) will give the same result everytime)
set.seed(1)

#Build the model with svm together with scaling and center the data
#add the training control and metrics to be used
SvmModel<-train(x=trainSmall[,1:7], y=trainSmall$Class, method="svmRadial",
               preProc=c("center", "scale"), trControl=trCtrl, metric="ROC")



#Plot the variable importance for the model for each class
pdf(file="VarImportance")
plot(varImp(SvmModel))
dev.off()

#Test/Predict with the model leaving out the class type
svmPred <- predict(SvmModel,testSmall[names(testSmall) != "Class"])

#Confusion matrix to get statistics on how good the predictior peformed on testdata
ConfMatrix<-confusionMatrix(svmPred, testSmall$Class)

