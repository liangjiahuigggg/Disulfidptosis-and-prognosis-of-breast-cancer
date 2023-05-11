#install.packages("pheatmap")


library(pheatmap)        #???ð?
expFile="uniSigExp.txt"            #?????????ļ?
geneCluFile="geneCluster.txt"      #???????͵Ľ????ļ?
crgCluFile="CRGcluster.txt"        #ͭ???????͵Ľ????ļ?
cliFile="clinical.txt"             #?ٴ??????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\32.geneHeatmap")     #???ù???Ŀ¼

#??ȡ?????????ļ???��???????ļ?
exp=read.table(expFile, header=T, sep="\t", check.names=F, row.names=1)
crgClu=read.table(crgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
sameSample=intersect(row.names(exp), row.names(crgClu))
exp=exp[sameSample,,drop=F]
expData=cbind(exp, geneCluster=geneClu[sameSample,], CRGcluster=crgClu[sameSample,])
Project=gsub("(.*?)\\_.*", "\\1", rownames(expData))
rownames(expData)=gsub("(.*?)\\_(.*?)", "\\2", rownames(expData))
expData=cbind(expData, Project)

#?ϲ??ٴ?????
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
cli[,"Age"]=ifelse(cli[,"Age"]=="unknow", "unknow", ifelse(cli[,"Age"]>65,">65","<=65"))
sameSample=intersect(row.names(expData), row.names(cli))
expData=expData[sameSample,,drop=F]
cli=cli[sameSample,,drop=F]
data=cbind(expData, cli)

#??ȡ??ͼ????
data=data[order(data$geneCluster),]
Type=data[,((ncol(data)-2-ncol(cli)):ncol(data))]
data=t(data[,1:(ncol(expData)-3)])

#??Type <- Type[,c(-3,-5)]????ͼע?͵???ɫ
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
ann_colors=list()
CRGcol=bioCol[1:length(levels(factor(Type$CRGcluster)))]
names(CRGcol)=levels(factor(Type$CRGcluster))
ann_colors[["CRGcluster"]]=CRGcol
GENEcol=bioCol[1:length(levels(factor(Type$geneCluster)))]
names(GENEcol)=levels(factor(Type$geneCluster))
ann_colors[["geneCluster"]]=GENEcol

#??ͼ???ӻ?
pdf("heatmap.pdf", height=6, width=8)
pheatmap(data,
         annotation=Type,
         annotation_colors = ann_colors,
         color = colorRampPalette(c(rep("blue",5), "white", rep("red",5)))(50),
         cluster_cols =F,
         cluster_rows =T,
         scale="row",
         show_colnames=F,
         show_rownames=F,
         fontsize=6,
         fontsize_row=6,
         fontsize_col=6)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

