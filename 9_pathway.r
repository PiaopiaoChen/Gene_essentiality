setwd("")
mydata <- read.table("1_all_changed_gene",header=T)
strains <- c("RM11","DBVPG1107","YJM975","DBVPG4651","WE372","X322134S","S288C","YJM326","YJM454","HN8","HN16","UWOPS83.787.3","YPS128","BJ20","BJ14","CLIB219")
rownames(mydata) <- mydata$ORF
mydata <- mydata[-1]
colnames(mydata) <- strains
mydata2 <- data.frame(t(mydata))
gene <- colnames(mydata2)


data <- read.table(file = "11_pathway_hit_with_2_subunits",sep="\t")


temp4 <- NULL
for (i in 1:nrow(data)){
 temp <- unlist(strsplit(data[i,2],"_"))
 temp2 <- data.frame(t(combn(1:length(temp),2)))

 for (j in 1:nrow(temp2)){
     temp3 <- c(temp[temp2[j,1]],temp[temp2[j,2]])
	 temp4 <- rbind(temp4,temp3)
 }
}
complex_change <- data.frame(temp4)


#complex[which(complex$X1 == "YLR312W-A"),1] <- "YLR312W.A"
#complex[which(complex$X2 == "YLR312W-A"),2] <- "YLR312W.A"

########################################

data_all_complex <- read.table(file = "11_pathway_hit_all_gene_list")
data_all_complex2 <- unique(data_all_complex)
com <- data.frame(t(combn(1:nrow(data_all_complex2),2)))

temp_all <- NULL
for (i in 1:nrow(com)){
   temp <- c(data_all_complex2$V1[com[i,1]],data_all_complex2$V1[com[i,2]])
   temp_all <- rbind(temp_all,temp)
}
complex_all <- data.frame(temp_all)

correlation_matrix_others  <- NULL
complex_change_rm_redundant <- NULL
for (i in 1:nrow(complex_all)){
  if ((length(which(complex_all[i,1] == complex_change[,1]  &  complex_all[i,2] == complex_change[,2])) == 0)  & (length(which(complex_all[i,1] == complex_change[,2]  &  complex_all[i,2] == complex_change[,1])) == 0)){
     correlation_matrix_others <- rbind(correlation_matrix_others,complex_all[i,])
  }else{
     complex_change_rm_redundant <- rbind(complex_change_rm_redundant,complex_all[i,])
  }
}
correlation_matrix_complex_gene_pairs <- data.frame(correlation_matrix_others)
complex_change_rm_redundant <- data.frame(complex_change_rm_redundant)

write.table(complex_change_rm_redundant,file="11_1_same_subunits_pair",append=FALSE,row.names =FALSE,quote = FALSE,sep="\t")
write.table(correlation_matrix_complex_gene_pairs,file="11_2_pathway_pair",append=FALSE,row.names =FALSE,quote = FALSE,sep="\t")


############# not_complex_gene

data_all_gene <- read.table(file = "1_all_changed_gene",header=T)
all_changed_list <- data_all_gene$ORF
changed_not_complex <- setdiff(all_changed_list,data_all_complex2$V1)

com <- data.frame(t(combn(1:length(changed_not_complex),2)))

temp_all <- NULL
for (i in 1:nrow(com)){
   temp <- c(changed_not_complex[com[i,1]],changed_not_complex[com[i,2]])
   temp_all <- rbind(temp_all,temp)
}
not_complex_all <- data.frame(temp_all)
write.table(not_complex_all,file="11_4_not_pathway_pair",append=FALSE,row.names =FALSE,quote = FALSE,sep="\t")

all_three <- rbind(complex_change_rm_redundant,correlation_matrix_complex_gene_pairs,not_complex_all)

###last group one complex one not
com <- data.frame(t(combn(1:567,2)))
temp_all <- NULL
for (i in 1:nrow(com)){
   temp <- c(all_changed_list[com[i,1]],all_changed_list[com[i,2]])
   temp_all <- rbind(temp_all,temp)
}
changed_all <- data.frame(temp_all)


mix_group  <- NULL
for (i in 1:nrow(changed_all)){
  if ((length(which(changed_all[i,1] == all_three[,1]  &  changed_all[i,2] == all_three[,2])) == 0)  & (length(which(changed_all[i,1] == all_three[,2]  &  changed_all[i,2] == all_three[,1])) == 0)){
     mix_group <- rbind(mix_group,changed_all[i,])
  }
}
mix_group <- data.frame(mix_group)
write.table(mix_group,file="11_3_one_pathway_one_not_pair",append=FALSE,row.names =FALSE,quote = FALSE,sep="\t")




####correlation calculated in perl 9_get_correlation.pl

cor1 <- read.table(file = "11_1_same_subunits_pair_correlation")
cor2 <- read.table(file = "11_2_pathway_pair_correlation")
cor3 <- read.table(file = "11_3_one_pathway_one_not_pair_correlation")
cor4 <- read.table(file = "11_4_not_pathway_pair_correlation")

mean(cor1$V3)
mean(cor2$V3)
mean(cor3$V3)
mean(cor4$V3)

#cor_data <- c(cor1$V3,cor2$V3,cor3$V3,cor4$V3)
#cor_sample <- c(rep("g1",nrow(cor1)),rep("g2",nrow(cor2)),rep("g3",nrow(cor3)),rep("g4",nrow(cor4)))
#cordata <- data.frame(cor_data,cor_sample)
#cordata$cor_sample <- factor(cordata$cor_sample,levels=unique(cordata$cor_sample))

cordata <- rbind(cor1,cor2,cor3,cor4)
cordata$cor_sample <- c(rep("g1",nrow(cor1)),rep("g2",nrow(cor2)),rep("g3",nrow(cor3)),rep("g4",nrow(cor4)))
cordata$cor_sample <- factor(cordata$cor_sample,levels=unique(cordata$cor_sample))
colnames(cordata) <- c("gene1","gene2","cor_data","cor_sample")

write.table(cordata,file="11_pathway_correlation",append=FALSE,row.names =FALSE,quote = FALSE,sep="\t")
#11_pathway_corr_group1_4  5.2 X 4.8
library(ggplot2)
ggplot(cordata,aes(cor_sample,cor_data,color=cor_sample,group=cor_sample))+geom_boxplot(width=0.5) +labs(x="complex membership", y="Correlation") + theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+theme(panel.background = element_blank(), axis.line = element_line(),legend.position="none")+scale_y_continuous(limits=c(-0.8,0.8),breaks=seq(-0.8,0.8,0.4))

wilcox.test(cordata[which(cordata$cor_sample == "g1"),3], cordata[which(cordata$cor_sample == "g2"),3])
#p-value = 2.506e-13
wilcox.test(cordata[which(cordata$cor_sample == "g1"),3], cordata[which(cordata$cor_sample == "g3"),3])
#p-value = 2.387e-11
wilcox.test(cordata[which(cordata$cor_sample == "g1"),3], cordata[which(cordata$cor_sample == "g4"),3])
#p-value = 1.164e-10

##############################


#11_pie_pathway_yes
mydata <- data.frame(group = c("Pathway genes & change", "Pathway genes & unchange"),value = c(110,664))
library(ggplot2)
ggplot(mydata, aes(x="", y=value, fill=group))+geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  axis.text.x=element_blank(),
  plot.title=element_text(size=14, face="bold")
  )
#11_pie_pathway_not
mydata2 <- data.frame(group = c("Non-pathway genes & change","Non-pathway genes & unchange"),value = c(457, 4583))
ggplot(mydata2, aes(x="", y=value, fill=group))+geom_bar(width = 1, stat = "identity")+ coord_polar("y", start=0)+theme_minimal()+
  theme(
  axis.title.x = element_blank(),
  axis.title.y = element_blank(),
  panel.border = element_blank(),
  panel.grid=element_blank(),
  axis.ticks = element_blank(),
  axis.text.x=element_blank(),
  plot.title=element_text(size=14, face="bold")
  )
