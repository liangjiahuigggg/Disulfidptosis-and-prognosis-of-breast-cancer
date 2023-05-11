#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggpubr")


#???ð?
library(limma)
library(ggpubr)
crgCluFile="CRGcluster.txt"        #ͭ???????͵Ľ????ļ?
geneCluFile="geneCluster.txt"      #???????͵Ľ????ļ?
scoreFile="risk.all.txt"           #?????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\36.clusterRisk")     #???ù???Ŀ¼

#??ȡ?????ļ?
crgClu=read.table(crgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
score=read.table(scoreFile, header=T, sep="\t", check.names=F, row.names=1)
score$riskScore[score$riskScore>quantile(score$riskScore,0.99)]=quantile(score$riskScore,0.99)

#?ϲ?????
twoCluster=cbind(crgClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(score))
data=cbind(score[sameSample,,drop=F], twoCluster[sameSample,,drop=F])


#######ͭ???????͵?????ͼ########
#???ñȽ???
data$CRGcluster=factor(data$CRGcluster, levels=levels(factor(data$CRGcluster)))
group=levels(factor(data$CRGcluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#????ͼ????ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$CRGcluster)))]
	
#????boxplot
boxplot=ggboxplot(data, x="CRGcluster", y="riskScore", color="CRGcluster",
			      xlab="CRGcluster",
			      ylab="Risk score",
			      legend.title="CRGcluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#????ͼƬ
pdf(file="CRGcluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######ͭ???????͵?????ͼ########


#######???????͵?????ͼ########
#???ñȽ???
data$geneCluster=factor(data$geneCluster, levels=levels(factor(data$geneCluster)))
group=levels(factor(data$geneCluster))
comp=combn(group, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#??????ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data$geneCluster)))]
	
#????boxplot
boxplot=ggboxplot(data, x="geneCluster", y="riskScore", color="geneCluster",
			      xlab="geneCluster",
			      ylab="Risk score",
			      legend.title="geneCluster",
			      palette=bioCol,
			      add = "jitter")+ 
	stat_compare_means(comparisons = my_comparisons)
	
#????ͼƬ
pdf(file="geneCluster.pdf", width=5, height=4.5)
print(boxplot)
dev.off()
#######???????͵?????ͼ########


######?⿡??ʦ΢??: eduBio

