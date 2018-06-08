rm(list=ls(all=TRUE))
setwd("/Users/syniu/ncbi/public/sra/")
library(e1071)
library(caret)
library(mlbench)
# Import training and testing dataset for SVM
positive_training <- read.delim("SimulatedPositiveTUMatrix.txt")
negative_training <- read.delim("SimulatedNegativeTUMatrix.txt")
positive_strand_testing <- read.delim("TargetPositiveTUMatrix.txt")
negative_strand_testing <- read.delim("TargetNegativeTUMatrix.txt")

pos_train_df <- as.data.frame(positive_training)
str(pos_train_df)
neg_train_df <- as.data.frame(negative_training)
str(neg_train_df)
pos_str_test_df <- as.data.frame(positive_strand_testing)
neg_str_test_df <- as.data.frame(negative_strand_testing)


# add labels in the last column, 1 is presence and 0 is absence of TU
pos_train_df$target_variable <- rep(1,nrow(pos_train_df))
neg_train_df$target_variable <- rep(0,nrow(neg_train_df))
train_data <- rbind(pos_train_df, neg_train_df)
#real_data <- rbind(pos_str_test_df, neg_str_test_df)

# Convert target variable to factors 
train_data[["target_variable"]] <- factor(train_data[["target_variable"]])

# Remove rows with NA
train_data <- train_data[complete.cases(train_data),]


# Dataset summarized details
summary(train_data)
#summary(real_data)

# Slice data into training and testing dataset
intrain <- createDataPartition(y = train_data$target_variable, p= 0.7, list = FALSE)
training <- train_data[intrain,]
testing <- train_data[-intrain,]

# 
# ### Feature Selection
# 
# # prepare training scheme
# control <- trainControl(method="repeatedcv", number=10, repeats=3)
# # train the model
# model <- train(target_variable~., data=train_data, method="lvq", preProcess="scale", trControl=control)
# # estimate variable importance
# importance <- varImp(model, scale=FALSE)
# # summarize importance
# print(importance)
# # plot importance
# plot(importance)
# 
# 
# # define the control using a random forest selection function
# control <- rfeControl(functions=rfFuncs, method="cv", number=10)
# # run the RFE algorithm
# results <- rfe(train_data[,1:17], train_data[,18], sizes=c(1:17), rfeControl=control)
# # summarize the results
# print(results)
# # list the chosen features
# predictors(results)
# # plot the results
# plot(results, type=c("g", "o"))
# 
# 
# featurePlot(x = train_data[, 1:17], 
#             y = train_data$target_variable ,
#             plot = "box",
#             strip=strip.custom(par.strip.text=list(cex=.7)),
#             scales = list(x = list(relation="free"), 
#                           y = list(relation="free")))
# 
# featurePlot(x = train_data[, 1:17], 
#             y = train_data$target_variable, 
#             plot = "density",
#             strip=strip.custom(par.strip.text=list(cex=.7)),
#             scales = list(x = list(relation="free"), 
#                           y = list(relation="free")))

# Training the Linear SVM model 

trctrl <- trainControl(method = "repeatedcv", number = 10, repeats = 3)

svm_Linear <- train(target_variable ~ GapLongest4 + GeneFoldChange + GapProportion4 + GapLongest3 +
                    GapProportion2 + GapProportion1 + ExpressionSD2 + ExpressionMean1 + GapLongest1
                    , method = "svmLinear", data = training,
                    trControl=trctrl,
                    preProcess = c("center", "scale"),
                    tuneLength = 10)

svm_Linear

# Testing set prediction
test_pred <- predict(svm_Linear, newdata = testing)
test_pred

# Confusion matrix
confusionMatrix(test_pred, testing$target_variable)

# Use to model to predict on real dataset
#real_predict <- predict(svm_Linear, newdata = real_data)
#real_predict
pos_str_test_df <- replace(pos_str_test_df, is.na(pos_str_test_df), 0)
SVMForward.Target.scaled.Predict <- predict(svm_Linear, newdata = pos_str_test_df)
neg_str_test_df <- replace(neg_str_test_df, is.na(neg_str_test_df), 0)
SVMReverse.Target.scaled.Predict <- predict(svm_Linear, newdata = neg_str_test_df)
### Result plot 
# TU.MergeNFinalTableNPlot-Qin_Chou071315
# predict on real data by forward and reverse separately, and start from line 40

# Read .NA file
SelectedRNAseqData<-read.table("Ecoli.NA", head=F) 

# Read gff file 
x=read.delim("GCF_000005845.2_ASM584v2_genomic.gff", skip=7, header = F) 
x1<-x[x$V3=="gene",]
x2<-within(x1[1,], y<-data.frame(do.call('rbind',strsplit(as.character(x1[1,]$V9), ';', fixed=TRUE))))
x3<-within(x2, y<-data.frame(do.call('rbind',strsplit(as.character(x2$y$X2), '=', fixed=TRUE))))
Ngene = data.frame(x3$y$X2, x3$V4, x3$V5, x3$V7)
for (i in 2:nrow(x1)){
  x2<-within(x1[i,], y<-data.frame(do.call('rbind',strsplit(as.character(x1[i,]$V9), ';', fixed=TRUE))))
  x3<-within(x2, y<-data.frame(do.call('rbind',strsplit(as.character(x2$y$X2), '=', fixed=TRUE))))
  x4<-data.frame(x3$y$X2, x3$V4, x3$V5, x3$V7)
  Ngene = rbind(Ngene,x4)
}
colnames(Ngene)=c("ID", "Start", "End", "Strand") # Set Column Names
Ngene=Ngene[order(Ngene[,2]),] #Order the list by start positions of each gene


#Select Genes 
SelectedForwardGenes=Ngene[Ngene[,4]=="+",] #Select all genes on forward strand
dim(SelectedForwardGenes)
SelectedReverseGenes=Ngene[Ngene[,4]=="-",] #Select all genes on reverse strand
dim(SelectedReverseGenes)
endIndexForwardTargetPredict=nrow(SelectedForwardGenes)-1

# Initialize the signal matrix for the forward Results
#filename = paste (as.character(Result[1,]),SVMForward.Target.scaled.Predict, sep = "")
#temp = read.table(SVMForward.Target.scaled.Predict)
MatrixRow = length(SVMForward.Target.scaled.Predict)
MatrixColumn = 1
MatrixForward = mat.or.vec(MatrixRow, MatrixColumn+2)

# Assign values to the signal matrix
for (i in 1:MatrixColumn){
  MatrixForward[,i] = as.matrix(SVMForward.Target.scaled.Predict)
}
mode(MatrixForward) = 'numeric'
# parse the frequency
for (i in 1:MatrixRow){
  positiveNumber = length(MatrixForward[i,MatrixForward[i,]==1])
  negativeNumber = MatrixColumn - positiveNumber
  if (positiveNumber>=negativeNumber){
    MatrixForward[i,MatrixColumn+1] = 1
  }else{
    MatrixForward[i,MatrixColumn+1] = -1
  }
  MatrixForward[i,MatrixColumn+2] = max(positiveNumber,negativeNumber)/MatrixColumn
}

# initialize the signal matrix for the reverse Results
MatrixRow = length(SVMReverse.Target.scaled.Predict)
MatrixReverse = mat.or.vec(MatrixRow, MatrixColumn+2)

# assign values to the signal matrix
for (i in 1:MatrixColumn){
  MatrixReverse[,i] = as.matrix(SVMReverse.Target.scaled.Predict)
}
mode(MatrixReverse) = 'numeric'

# Parse the frequency
for (i in 1:MatrixRow){
  positiveNumber = length(MatrixReverse[i,MatrixReverse[i,]==1])
  negativeNumber = MatrixColumn - positiveNumber
  if (positiveNumber>=negativeNumber){
    MatrixReverse[i,MatrixColumn+1] = 1
  }else{
    MatrixReverse[i,MatrixColumn+1] = -1
  }
  MatrixReverse[i,MatrixColumn+2] = max(positiveNumber,negativeNumber)/MatrixColumn
}



# Generate Target TU SVM format
GenerateOneTargetTU<-function(InputGeneList, i){
  OneTUList=data.frame()
  OneTUList=c("5'Gene", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][1,3], as.character(InputGeneList[i:(i+1),][1,4]));
  OneTUList=rbind(OneTUList, c("IntergenicRegion", (InputGeneList[i:(i+1),][1,3]+1), (InputGeneList[i:(i+1),][2,2]-1),as.character( InputGeneList[i:(i+1),][1,4])));
  OneTUList=rbind(OneTUList, c("3'Gene", as.numeric(as.character(InputGeneList[i:(i+1),][2,2])), InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][2,4])));
  OneTUList=rbind(OneTUList, c("FullTU", InputGeneList[i:(i+1),][1,2], InputGeneList[i:(i+1),][2,3], as.character(InputGeneList[i:(i+1),][1,4])));
  row.names(OneTUList)=NULL
  OneTUList=as.data.frame(OneTUList)
  OneTUList[,2]=as.numeric(as.character(OneTUList[,2]))
  OneTUList[,3]=as.numeric(as.character(OneTUList[,3]))
  colnames(OneTUList)=c("ID", "Start", "End", "Strand") # Set Column Names
  return(OneTUList)
}



# Generate the final TU tables
GenerateFinalTUTable<-function(TargetTUsPrediction, SelectedGenes, TargetTUsPredictionFreq){
  colnames(TargetTUsPrediction) = "TUresult"
  TargetTUName1Name2StartEnd=data.frame()
  j=1
  for(i in 1:(nrow(SelectedGenes)-1)){
    TargetTUName1Name2StartEnd=rbind(TargetTUName1Name2StartEnd, cbind(SelectedGenes[i,1], SelectedGenes[(i+1),1],GenerateOneTargetTU(SelectedGenes,i)[4,2:3]))
    j=j+1
  }
  TargetTUsPrediction = (cbind(TargetTUsPrediction, TargetTUName1Name2StartEnd))
  FinalTUAll=data.frame()
  FinalTUOne=data.frame()
  FinalTUOneNameStart=as.character(TargetTUsPrediction[1,2])
  FinalTUOneNameEnd=as.character("")
  FinalTUOnePosStart=TargetTUsPrediction[1,4]
  FinalTUOnePosEnd=TargetTUsPrediction[1,5]
  FinalTUOneFreq=0
  FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneFreq)
  flag=TargetTUsPrediction[1,1]
  geneNum=0
  for(i in c(1:nrow(TargetTUsPrediction))){
    if(TargetTUsPrediction[i,1]==1){  #TU result = 1
      FinalTUOne[,2]=paste(FinalTUOne[,2],TargetTUsPrediction[i,3], sep=" ")
      FinalTUOne[,4]=TargetTUsPrediction[i,5]
      FinalTUOne[,5] = FinalTUOne[,5] + TargetTUsPredictionFreq[i,1]
      flag=1
      geneNum = geneNum+1
      if (i==nrow(TargetTUsPrediction)){
        FinalTUOne[,5] = (FinalTUOne[,5])/(geneNum)
      }
    }else{                                   #TU result = -1
      FinalTUOne[,5] = (FinalTUOne[,5] + TargetTUsPredictionFreq[i,1])/(geneNum+1)
      FinalTUAll=rbind(FinalTUAll,FinalTUOne)
      FinalTUOne=NULL
      flag=-1
      geneNum = 1
      FinalTUOneFreq = TargetTUsPredictionFreq[i,1]
      FinalTUOneNameStart=as.character(TargetTUsPrediction[i+1,2])
      FinalTUOneNameEnd=as.character("")
      FinalTUOnePosStart=TargetTUsPrediction[i+1,4]
      FinalTUOnePosEnd=TargetTUsPrediction[i,5]
      FinalTUOne=data.frame(FinalTUOneNameStart,FinalTUOneNameEnd,FinalTUOnePosStart,FinalTUOnePosEnd,FinalTUOneFreq)
    }
  }
  if(FinalTUOne[,3] == TRUE && FinalTUOne[,4] == TRUE){
    FinalTUAll=rbind(FinalTUAll,FinalTUOne)
  }
  return(FinalTUAll)
}


ForwardTargetTUsPrediction <- as.matrix(MatrixForward[,MatrixColumn+1])
mode(ForwardTargetTUsPrediction) <- 'numeric'
ReverseTargetTUsPrediction <- as.matrix(MatrixReverse[,MatrixColumn+1]) 
mode(ReverseTargetTUsPrediction) <- 'numeric'
ForwardTargetTUsPredictionFreq <- as.matrix(MatrixForward[,MatrixColumn+2])
mode(ForwardTargetTUsPredictionFreq) <- 'numeric'
ReverseTargetTUsPredictionFreq <- as.matrix(MatrixReverse[,MatrixColumn+2]) 
mode(ReverseTargetTUsPredictionFreq) <- 'numeric'
ForwardTUtable<-cbind(GenerateFinalTUTable(ForwardTargetTUsPrediction,SelectedForwardGenes,ForwardTargetTUsPredictionFreq),Strand="+")
ReverseTUtable<-cbind(GenerateFinalTUTable(ReverseTargetTUsPrediction,SelectedReverseGenes,ReverseTargetTUsPredictionFreq),Strand="-")
head (ForwardTUtable)
head (ReverseTUtable)
AllTUtable<-merge(ForwardTUtable, ReverseTUtable, all=TRUE)
AllTUtable[,2] <- paste(AllTUtable[,1], AllTUtable[,2], sep="")
AllTUtable <- AllTUtable[,c(3,4,6,5,2)]
AllTUtable <- AllTUtable[order(AllTUtable$FinalTUOnePosStart),]
row.names(AllTUtable) <- paste("TU",formatC(c(1:nrow(AllTUtable)), width=4,flag="0"), sep="")
FinalTUTableFileName = "FinalTUTable_forPlot"
write.table(AllTUtable, file = FinalTUTableFileName, sep = "\t", col.names=FALSE)

system(paste("sed -e 's/\"//g' ",FinalTUTableFileName," > ",FinalTUTableFileName,".clean",sep=""))
system(paste("mv ",FinalTUTableFileName,".clean ",FinalTUTableFileName,sep=""))



PlotSelectedRangeData <- function(Start, End, RNAseqData, GeneList, ID){
  #wcc091015#par(mfcol=c(3,1), mar=c(0,4,0,2)+0.1, oma=c(5,0,3,0)+0.1)
  par(mfcol=c(2,1), mar=c(0,4,0,2)+0.1, oma=c(5,0,3,0)+0.1)
  #wcc091015#plot(c(Start:End), RNAseqData[Start:End,1],type="h",col="blue",xlab="",xaxt="n", ylab="log2 Forward Signal")
  plotExtension=800;
  Start=max(c(1, Start-plotExtension))
  End=End+plotExtension
  plot(c(Start:End), RNAseqData[Start:End,1],type="h",col="blue",xlab="",xaxt="n", ylab="log2 RNA-seq Signal")
  plot(x=c(Start, End), y=c(0,100), type="n", xlab="", ylab="Gene Annotation", xaxt="n", yaxt="n")
  if(length(GeneList[GeneList[,4]=="-", 2][(GeneList[GeneList[,4]=="-", 2]>=Start) & (GeneList[GeneList[,4]=="-", 3]<=End)])>0){
    rect(GeneList[GeneList[,4]=="-", 2][(GeneList[GeneList[,4]=="-", 2]>=Start) & (GeneList[GeneList[,4]=="-", 3]<=End)], 45, GeneList[GeneList[,4]=="-", 3][(GeneList[GeneList[,4]=="-", 2]>=Start) & (GeneList[GeneList[,4]=="-", 3]<=End)], 50, col=palette(), border=0 )
    text(GeneList[GeneList[,4]=="-", 2][(GeneList[GeneList[,4]=="-", 2]>=Start) & (GeneList[GeneList[,4]=="-", 3]<=End)], 45, GeneList[GeneList[,4]=="-", 1][(GeneList[GeneList[,4]=="-", 2]>=Start) & (GeneList[GeneList[,4]=="-", 3]<=End)], cex=1, pos=4, srt=315)#adj=c(0,0),pos=1, srt=45
  }
  if(length(GeneList[GeneList[,4]=="+", 2][(GeneList[GeneList[,4]=="+", 2]>=Start) & (GeneList[GeneList[,4]=="+", 3]<=End)])>0){
    rect(GeneList[GeneList[,4]=="+", 2][(GeneList[GeneList[,4]=="+", 2]>=Start) & (GeneList[GeneList[,4]=="+", 3]<=End)], 55, GeneList[GeneList[,4]=="+", 3][(GeneList[GeneList[,4]=="+", 2]>=Start) & (GeneList[GeneList[,4]=="+", 3]<=End)], 60, col=palette(), border=0 )
    text(GeneList[GeneList[,4]=="+", 2][(GeneList[GeneList[,4]=="+", 2]>=Start) & (GeneList[GeneList[,4]=="+", 3]<=End)], 55, GeneList[GeneList[,4]=="+", 1][(GeneList[GeneList[,4]=="+", 2]>=Start) & (GeneList[GeneList[,4]=="+", 3]<=End)], cex=1, pos=4, srt=315)#adj=c(0,0),pos=1, srt=45
  }
  text(Start, 57.5, "+", cex=1, pos=2)
  text(Start, 47.5, "-", cex=1, pos=2)
  text(End, 57.5, "+", cex=1, pos=4)
  text(End, 47.5, "-", cex=1, pos=4)
  #wcc091015#plot(x=c(Start:End),RNAseqData[Start:End,2],type="h",col="red",xlab="",xaxt="n", ylab="log2 Reverse Signal")
  mtext( ID, outer = TRUE )
}


log2SelectedRNAseqData <- log2(SelectedRNAseqData+1)
for(i in c(1:nrow(AllTUtable))){
  #FinalPerformanceFileName = paste(sub("SVMForward","",args[4], perl=TRUE),"FinalPerformance", sep="")
  pdf(file="result_table.pdf")
  PlotSelectedRangeData(AllTUtable$FinalTUOnePosStart[i], AllTUtable$FinalTUOnePosEnd[i], log2SelectedRNAseqData, Ngene, row.names(AllTUtable)[i])
  dev.off()
  # png(file=paste(sub("SVMForward","",args[1], perl=TRUE),"SignalPlot.",row.names(AllTUtable)[i],".png",sep=""), bg = "transparent")
  # PlotSelectedRangeData(AllTUtable$FinalTUOnePosStart[i], AllTUtable$FinalTUOnePosEnd[i], log2SelectedRNAseqData, Ngene, row.names(AllTUtable)[i])
  # dev.off()
}
