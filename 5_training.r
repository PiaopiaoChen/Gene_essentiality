library(randomForest)
library(caret)
library(pROC)
library(ggplot2)

dir0=""
setwd(dir0)
fea_data1 <- read.table("11_s31_features",sep="\t",header=F,fill=TRUE)
colnames(fea_data1) <- c("gene", "primaryIdentifier","chromosome","eff_len","hits","reads","max_l","eff_len_nei10000","nei_hits","up100_hits","up200_hits","bins_wo_tn","upstream","downstream","upstream_r","downstream_r")


gene_info <- read.table("all_gene_info",sep="\t",header=T,stringsAsFactors=FALSE,quote = "")
gene_info <- data.frame(gene=gene_info$ORF,exon_number = gene_info$exon_num)

fea_data1<- merge(fea_data1,gene_info,by="gene")
fea_data1 <- fea_data1[(!grepl("^Q",fea_data1$gene)) & (!grepl("^R",fea_data1$gene)),]
fea_data2 <- fea_data1[which(fea_data1[,4] > 300),]
rownames(fea_data2) <- c(1:nrow(fea_data2))


#calculate nei_r
nei_r <- NULL
for (i in 1:nrow(fea_data2)){
   if(fea_data2$nei_hits[i]!=0 && fea_data2$eff_len_nei10000[i] !=0){
      temp <- (fea_data2$hits[i]/fea_data2$eff_len[i])/(fea_data2$nei_hits[i]/fea_data2$eff_len_nei10000[i])
	}
	else{
	   temp <- 0
	}
	nei_r <- c(nei_r,temp)
}

#calculate features used in random forest  
fea_data3 <- data.frame(gene=fea_data2$gene, primaryIdentifier = fea_data2$primaryIdentifier, chromosome = fea_data2$chromosome, eff_len = fea_data2$eff_len,hits= fea_data2$hits,reads = fea_data2$reads, max_l=fea_data2$max_l,hit_r = round(fea_data2$hits/fea_data2$eff_len,3), reads_r=round(fea_data2$reads/fea_data2$eff_len,3), max_r= round(fea_data2$max_l/fea_data2$eff_len,3), nei_r=round(nei_r,3),up100=fea_data2$up100_hits, bins=fea_data2$bins_wo_tn,up_down_r= (fea_data2$upstream_r +fea_data2$downstream_r)/2)

fea_data3$primaryIdentifier <- factor(fea_data3$primaryIdentifier)
train <- sample(nrow(features), nrow(features)*0.5)
TrainSet <- features[train, ]
ValidSet <- features[-train, ]

trcontrol <- trainControl(method = 'cv', number=10, savePredictions = T, classProbs = TRUE, summaryFunction = twoClassSummary,returnResamp="all")
train_forest <-   train(primaryIdentifier ~ hits + reads + max_l  + hit_r + reads_r + max_r  + nei_r + up100 + bins, data=TrainSet, method = "rf", trControl = trcontrol,metric="ROC") 

importance <- varImp(train_forest, scale=FALSE)
importance[[1]]$Overall <- importance[[1]]$Overall / sum(importance[[1]]$Overall)
# plot importance
plot(importance)


# Predicting on Validation set
predValid1 <- predict(train_forest, ValidSet)   #type = "class" 
predValid2 <- predict(train_forest, ValidSet, type = "prob")
confusionMatrix(predict(train_forest,ValidSet),ValidSet$primaryIdentifier)
ValidSet_rf <- data.frame(gene=ValidSet$gene,rf=predValid1,rf_ES = predValid2$ES,primaryIdentifier=ValidSet$primaryIdentifier)





