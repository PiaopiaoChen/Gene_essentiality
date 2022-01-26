
library("castor")
library(Rcpp)
library(grid)
library(vcd)
library(reshape2)
library(ggplot2)
 

setwd("")
all_interaction_data <- read.table("12_genetic_interaction",header=F, sep="\t")

all_interaction_data$merge <- apply(all_interaction_data[,c(1:2)],1,function(x) paste(sort(x),collapse='_'))
all_interaction_rmd <- all_interaction_data[!duplicated(all_interaction_data$merge),1:4]

all_gene_list <- c(all_interaction_rmd$V1,all_interaction_rmd$V2)

all_gene_hit_num <- data.frame(table(all_gene_list))
colnames(all_gene_hit_num)<- c("gene","nums")


#changed_gene_list
  
tree <- read_tree(file="2_1_tree.nwk")

trait2 <- read.table("4_all_parsimony_estimation_missing_value",header=T, fill=TRUE)

trait <- data.frame(CLIB219=trait2$CLIB219,BJ14=trait2$BJ14,BJ20=trait2$BJ20,YPS128=trait2$YPS128,UWOPS83.787.3=trait2$UWOPS83.787.3,HN16=trait2$HN16,HN8=trait2$HN8,YJM454=trait2$YJM454,YJM326=trait2$YJM326,S288C=trait2$S288C,X322134S=trait2$X322134S,WE372=trait2$WE372,DBVPG4651=trait2$DBVPG4651,YJM975=trait2$YJM975,DBVPG1107=trait2$DBVPG1107,RM11=trait2$RM11)
rownames(trait) <- rownames(trait2)
cost_all <- NULL
for (i in 1:nrow(trait)){
x <- as.integer(trait[i,1:ncol(trait)])
cost <- asr_max_parsimony(tree, x)$total_cost
cost_all <- c(cost_all,cost)
}

changed_gene_list_with_cost <- data.frame(gene=rownames(trait),cost=cost_all, change="Changed")

###unchanged_gene_list
genes_change_or_not <- read.table("1_all_gene_list",header=F, fill=TRUE,sep="\t")
colnames(genes_change_or_not) <- c("gene","anno","change")
genes_change_or_not <- genes_change_or_not[,1:3]

unchanged_gene_list_with_cost <- data.frame(gene=genes_change_or_not[genes_change_or_not$change=="Unchanged",1],cost=0, change="Unchanged")

######
all_gene_cost <- rbind(changed_gene_list_with_cost,unchanged_gene_list_with_cost)

all_gene_cost_hit <- merge(all_gene_cost,all_gene_hit_num,by="gene",all.x=TRUE)
all_gene_cost_hit <- merge(all_gene_cost_hit,genes_change_or_not,by=c("gene","change"),all.x=TRUE)
all_gene_cost_hit$change <- as.factor(all_gene_cost_hit$change)


all_gene_cost_hit$anno2 <- "ok"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Unchanged" & all_gene_cost_hit$anno == "nonessential"),6] <- "monomorphic_nonessential_genes"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Unchanged" & all_gene_cost_hit$anno == "essential"),6] <- "monomorphic_essential_enes"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Changed"),6] <- "polymorphic_genes"

all_gene_cost_hit$anno2 <- factor(all_gene_cost_hit$anno2, levels=c("monomorphic_essential_enes","polymorphic_genes","monomorphic_nonessential_genes"))

#7.2 7  12_genetic_monomorphic_polymorphic_fig3b  #6X6

ggplot(all_gene_cost_hit,aes(x=anno2,y=nums,group=anno2,color=anno2))+geom_boxplot(width=0.5,outlier.shape = NA,notch = TRUE)+ coord_cartesian(ylim = c(0,800))+labs(y="Number of interactions") + theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+theme(panel.background = element_blank(), axis.line = element_line(),legend.position="none")

aa <- wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno2 == "monomorphic_essential_enes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno2 == "polymorphic_genes"),4])
aa$p.value
# 2.654e-05

bb <- wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno2 == "monomorphic_nonessential_genes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno2 == "polymorphic_genes"),4])
bb$p.value

cc <- wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno2 == "monomorphic_essential_enes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno2 == "monomorphic_nonessential_genes"),4])
cc$p.value

all_gene_cost_hit$anno3 <- "ok"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Unchanged" & all_gene_cost_hit$anno == "nonessential"),7] <- "monomorphic_nonessential_genes"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Unchanged" & all_gene_cost_hit$anno == "essential"),7] <- "monomorphic_essential_genes"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Changed" & all_gene_cost_hit$anno == "nonessential"),7] <- "polymorphic_nonessential_genes"
all_gene_cost_hit[which(all_gene_cost_hit$change == "Changed" & all_gene_cost_hit$anno == "essential"),7] <- "polymorphic_essential_genes"
all_gene_cost_hit$anno3 <- factor(all_gene_cost_hit$anno3, levels=c("monomorphic_essential_genes","polymorphic_essential_genes","polymorphic_nonessential_genes","monomorphic_nonessential_genes"))


#12_genetic_figS5
ggplot(all_gene_cost_hit,aes(x=anno3,y=nums,group=anno3,color=anno3))+geom_boxplot(width=0.5,outlier.shape = NA,notch = TRUE)+ coord_cartesian(ylim = c(0,800))+labs(y="Number of interactions") + theme(axis.title = element_text(size=10), axis.text.x = element_text(size=10), axis.text.y = element_text(size=10))+theme(panel.background = element_blank(), axis.line = element_line(),legend.position="none")

wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "monomorphic_essential_genes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "polymorphic_essential_genes"),4])
# 1.133e-06

dd <- wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "polymorphic_nonessential_genes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "monomorphic_nonessential_genes"),4])
dd$p.value
#p-value < 2.2e-16

wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "polymorphic_essential_genes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "polymorphic_nonessential_genes"),4])
# p-value = 2.718e-05

wilcox.test(all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "polymorphic_essential_genes"),4], all_gene_cost_hit[which(all_gene_cost_hit$anno3 == "monomorphic_nonessential_genes"),4])
#p-value = 0.02502


write.table(all_gene_cost_hit,file="12_genetic_interaction_result",append=FALSE,row.names =FALSE,quote = FALSE,sep="\t")

