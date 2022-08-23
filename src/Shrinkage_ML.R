################################
# nearest shrunken centroid
###############################
library(pamr)


##########################################################################################
#Usage
#########################################################################################
#Training
pamr.train(data, gene.subset=NULL, sample.subset=NULL,
           threshold = NULL, n.threshold = 30,
           scale.sd = TRUE, threshold.scale = NULL, se.scale = NULL, offset.percent = 50,
           hetero=NULL, prior = NULL, remove.zeros = TRUE, sign.contrast="both",
           ngroup.survival = 2)

#Cross-validation
pamr.cv(fit, data, nfold = NULL, folds = NULL,...)
################################################################################
#make a discovery matrix for pamr

discoveryMatrix.pamr <- gexprMatrix[ rowMedians(as.matrix(gexprMatrix)) >=(log2(1000)),]

results <- pamr.train(data=list(x=as.matrix(discoveryMatrix.pamr),y=groups$values),n.threshold=100)  #Groups are the labels of your samples

#Cross-validation
p <- pamr.cv(fit=results,data=list(x=as.matrix(discoveryMatrix.pamr),y=groups))

#Get the gene list
genes <- pamr.listgenes(results,data=list(x=as.matrix(discoveryMatrix.pamr),y=groups,genenames=rownames(discoveryMatrix.pamr),geneid=rownames(discoveryMatrix.pamr)),threshold=4.5,genenames=T)
#################################################
###############
#SVM classifier
###############
library(caret)

svmTuneGrid <- expand.grid(C = 2^(-8:-5))

#Train
mod0 <- train(y=as.factor(make.names(groups$values)),x=t(discoveryMatrix),
              method = "svmLinear",
              metric = "ROC",
              tuneGrid = svmTuneGrid,trControl = trainControl(method = "repeatedcv",
                                                              repeats = 100,number=5,
                                                              classProbs = TRUE,savePredictions=T,                                                          
                                                              summaryFunction = twoClassSummary))
getTrainPerf(mod0)
predictions <- predict.train(mod0,newdata=t(discoveryMatrix),type="raw")
confusionMatrix(predictions,reference=make.names(groups$values))
#----------------------------------------------------------------------




suppressWarnings(RNGversion("3.5.0"))
set.seed(120)

x <- matrix(rnorm(1000*20),ncol=20)
y <- sample(c(1:4),size=20,replace=TRUE)

mydata <- list(x=x,y=factor(y), geneid=as.character(1:nrow(x)),
               genenames=paste("g",as.character(1:nrow(x)),sep=""))

mytrain <-   pamr.train(mydata)
mycv <- pamr.cv(mytrain,mydata)

pamr.listgenes(mytrain, mydata, threshold=1.6)

genes <- pamr.listgenes(results,data=list(x=as.matrix(discoveryMatrix.pamr),y=groups,genenames=rownames(discoveryMatrix.pamr),geneid=rownames(discoveryMatrix.pamr)),threshold=4.5,genenames=T)
