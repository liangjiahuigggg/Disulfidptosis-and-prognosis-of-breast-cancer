#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("reshape2")
#install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(ggpubr)

expFile="cuproGeneExp.txt"         #?????????ļ?
geneCluFile="geneCluster.txt"      #???????͵Ľ????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\33.CRGclusterDiff")       #???ù???Ŀ¼

#??ȡ?????????ļ?
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=t(data)

#??ȡ???????͵Ľ????ļ?
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
sameSample=intersect(row.names(data), row.names(geneClu))
expClu=cbind(data[sameSample,,drop=F], geneClu[sameSample,,drop=F])

#??????ת????ggplot2?????ļ?
data=melt(expClu, id.vars=c("geneCluster"))
colnames(data)=c("geneCluster", "Gene", "Expression")

#????ͼ?ε???ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"geneCluster"])))]

#????????ͼ
p=ggboxplot(data, x="Gene", y="Expression", color = "geneCluster",
	     xlab="",
	     ylab="Gene expression",
	     legend.title="geneCluster",
	     palette = bioCol,
	     width=0.8)
p=p+rotate_x_text(45)
p1=p+stat_compare_means(aes(group=geneCluster),
	      symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", " ")),
	      label = "p.signif")

#????????ͼ
pdf(file="boxplot.pdf", width=7, height=5)
print(p1)
dev.off()



