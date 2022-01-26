setwd("C:/project/essentiality/transposon_16strains/hit_on_gene/5_random_forest/2_change")
library("castor")
library(Rcpp)

tree <- read_tree(file="2_1_tree.nwk")

trait2 <- read.table("1_all_with_missing_data_change",header=T, fill=TRUE,sep="\t")

trait <- data.frame(ORF=trait2$ORF,CLIB219=trait2$CLIB219,BJ14=trait2$BJ14,BJ20=trait2$BJ20,YPS128=trait2$YPS128,UWOPS83.787.3=trait2$UWOPS83.787.3,HN16=trait2$HN16,HN8=trait2$HN8,YJM454=trait2$YJM454,YJM326=trait2$YJM326,S288C=trait2$S288C,X322134S=trait2$X322134S,WE372=trait2$WE372,DBVPG4651=trait2$DBVPG4651,YJM975=trait2$YJM975,DBVPG1107=trait2$DBVPG1107,RM11=trait2$RM11)
#1 means essential; 2 means non-essential

# reconstruct all tip states via MPR
estimated_tip_states_all <- NULL
cost_all <- NULL
for (i in 1:nrow(trait)){

x <- as.integer(trait[i,2:ncol(trait)])
likelihoods = hsp_max_parsimony(tree, x, 2)$likelihoods
estimated_tip_states = max.col(likelihoods[1:16,])
estimated_tip_states_all <- rbind(estimated_tip_states_all,estimated_tip_states)
cost <- asr_max_parsimony(tree, estimated_tip_states)$total_cost
cost_all <- c(cost_all,cost)

}
estimated_tip_states_all <- data.frame(estimated_tip_states_all)
rownames(estimated_tip_states_all) <- trait[,1]
colnames(estimated_tip_states_all) <- colnames(trait)[2:ncol(trait)]

write.table(estimated_tip_states_all,"4_all_parsimony_estimation_missing_value",sep="\t",quote = FALSE, row.names = TRUE,col.names=TRUE)



##############estimate total cost

library("castor")
library(Rcpp)
library(grid)
library(vcd)
library(reshape2)
library(ggplot2)
  
tree <- read_tree(file="2_1_tree.nwk")

trait2 <- read.table("4_all_parsimony_estimation_missing_value",header=T, fill=TRUE)

trait <- data.frame(CLIB219=trait2$CLIB219,BJ14=trait2$BJ14,BJ20=trait2$BJ20,YPS128=trait2$YPS128,UWOPS83.787.3=trait2$UWOPS83.787.3,HN16=trait2$HN16,HN8=trait2$HN8,YJM454=trait2$YJM454,YJM326=trait2$YJM326,S288C=trait2$S288C,X322134S=trait2$X322134S,WE372=trait2$WE372,DBVPG4651=trait2$DBVPG4651,YJM975=trait2$YJM975,DBVPG1107=trait2$DBVPG1107,RM11=trait2$RM11)

cost_all <- NULL
for (i in 1:nrow(trait)){

x <- as.integer(trait[i,1:ncol(trait)])
cost <- asr_max_parsimony(tree, x)$total_cost
cost_all <- c(cost_all,cost)
}
table(cost_all)

#4753 total genes  5814-567=5247
total_change <- c(as.numeric(cost_all),rep(0,5247))

mean(total_change)* mean(total_change)/(var(total_change)- mean(total_change)) 
a <-  mean(total_change)
b <- var(total_change)
 (a*a)/(b-a)
 
#number of changes/number of genes
mean_rate <- sum(cost_all)/length(total_change)
expectation <-  dpois(0:7,lambda=mean_rate)*length(total_change)

expectation2 <-  round(dpois(0:7,lambda=mean_rate)*length(total_change))
expectation3 <- c(rep(0,expectation2[1]),rep(1,expectation2[2]),rep(2,expectation2[3]),rep(3,expectation2[4]),rep(4,expectation2[5]),rep(5,expectation2[6]),rep(6,expectation2[7]),rep(7,expectation2[8]))

library(EnvStats)
sigma0 <- sqrt(mean(total_change))
aa <- varTest(total_change, alternative = "greater", conf.level = 0.95, sigma0)

#with zero
mydata <- data.frame(table(total_change))
colnames(mydata) <- c("Var1","true")
mydata$exp <- expectation

colnames(mydata) <- c("num","true","expect")
all_data <- melt(mydata, id.vars=("num"))
all_data$value <- as.numeric(as.character(all_data$value))
all_data$num <- as.numeric(as.character(all_data$num))

ggplot(all_data,aes(x=num,y=value,fill=variable))+geom_bar(stat="identity",position = "dodge",width=.6)+theme(axis.text.x = element_text(angle = 45, hjust = 0.5, vjust = 0.5)) +  xlab("Number of changes")+ylab("Number of genes")+theme_bw()+theme(panel.border = element_blank(),panel.grid.major = element_blank(),panel.grid.minor = element_blank(),axis.line = element_line(colour = "black"))+ scale_fill_manual(values=c('#FF0000','#808080'))+scale_y_continuous(breaks=seq(0,5000,500))


 

