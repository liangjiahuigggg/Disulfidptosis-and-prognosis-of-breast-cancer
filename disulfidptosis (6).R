#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")
#BiocManager::install("GSEABase")
#BiocManager::install("GSVA")

#install.packages("ggpubr")


#???ð?
library(reshape2)
library(ggpubr)
library(limma)
library(GSEABase)
library(GSVA)

expFile="merge.txt"               #?????????ļ?
clusterFile="CRGcluster.txt"      #???ͽ????ļ?
gmtFile="immune.gmt"              #???߻??????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\25.ssGSEA")     #???ù???Ŀ¼

#??ȡ?????????ļ?,?????????ļ?????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)

#??ȡ???????ļ?
geneSets=getGmt(gmtFile, geneIdType=SymbolIdentifier())

#ssGSEA????
ssgseaScore=gsva(data, geneSets, method='ssgsea', kcdf='Gaussian', abs.ranking=TRUE)
#??ssGSEA???ֽ??н???
normalize=function(x){
  return((x-min(x))/(max(x)-min(x)))}
ssgseaScore=normalize(ssgseaScore)
#????ssGSEA???ֽ???
ssgseaOut=rbind(id=colnames(ssgseaScore), ssgseaScore)
write.table(ssgseaOut,file="ssGSEA.result.txt",sep="\t",quote=F,col.names=F)

#??ȡ???͵Ľ????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)

#???ݺϲ?(????????????ϸ?????ֽ??кϲ?)
ssgseaScore=t(ssgseaScore)
sameSample=intersect(row.names(ssgseaScore), row.names(cluster))
ssgseaScore=ssgseaScore[sameSample,,drop=F]
cluster=cluster[sameSample,,drop=F]
scoreCluster=cbind(ssgseaScore, cluster)

#??????ת????ggplot2?????ļ?
data=melt(scoreCluster, id.vars=c("CRGcluster"))
colnames(data)=c("CRGcluster", "Immune", "Fraction")

#????????ͼ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(levels(factor(data[,"CRGcluster"])))]
p=ggboxplot(data, x="Immune", y="Fraction", color="CRGcluster", 
     ylab="Immune infiltration",
     xlab="",
     legend.title="CRGcluster",
     palette=bioCol)
p=p+rotate_x_text(50)

#????ͼ??
pdf(file="boxplot.pdf", width=8, height=6.5)
p+stat_compare_means(aes(group=CRGcluster),symnum.args=list(cutpoints = c(0, 0.001, 0.01, 0.05, 1), symbols = c("***", "**", "*", "")),label = "p.signif")
dev.off()




