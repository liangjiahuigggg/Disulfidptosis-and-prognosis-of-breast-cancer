#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("ggplot2")


#???ð?
library(limma)
library(ggplot2)

expFile="cuproGeneExp.txt"       #?????????ļ?
clusterFile="CRGcluster.txt"     #???͵Ľ????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\26.PCA")     #???ù???Ŀ¼

#??ȡ?????ļ?,?????????ļ?????????
rt=read.table(expFile, header=T, sep="\t", check.names=F)
rt=as.matrix(rt)
rownames(rt)=rt[,1]
exp=rt[,2:ncol(rt)]
dimnames=list(rownames(exp),colnames(exp))
data=matrix(as.numeric(as.matrix(exp)),nrow=nrow(exp),dimnames=dimnames)
data=avereps(data)
data=data[rowMeans(data)>0,]
data=t(data)

#PCA????
data.pca=prcomp(data, scale. = TRUE)
pcaPredict=predict(data.pca)
write.table(pcaPredict, file="PCA.result.txt", quote=F, sep="\t")

#??ȡ?????ļ?
cluster=read.table(clusterFile, header=T, sep="\t", check.names=F, row.names=1)
CRGcluster=as.vector(cluster[,1])

#???÷??͵???ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
crgCluCol=bioCol[1:length(levels(factor(CRGcluster)))]

#????ͼ??
PCA=data.frame(PC1=pcaPredict[,1], PC2=pcaPredict[,2], CRGcluster=CRGcluster)
PCA.mean=aggregate(PCA[,1:2], list(CRGcluster=PCA$CRGcluster), mean)
pdf(file="PCA.pdf", width=6, height=4.6)
ggplot(data = PCA, aes(PC1, PC2)) + geom_point(aes(color = CRGcluster)) +
	scale_colour_manual(name="CRGcluster", values =crgCluCol)+
    theme_bw()+
    theme(plot.margin=unit(rep(1.5,4),'lines'))+
    annotate("text",x=PCA.mean$PC1, y=PCA.mean$PC2, label=PCA.mean$CRGcluster, cex=7)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

