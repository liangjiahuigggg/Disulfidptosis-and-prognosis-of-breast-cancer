#install.packages("ggpubr")


#引用包
library(ggpubr)
library(reshape2)

tmbFile="TMB.txt"             #肿瘤突变负荷文件
riskFile="risk.all.txt"       #风险文件
cluFile="geneCluster.txt"     #基因分型的结果文件
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\48.riskTMB")      #设置工作目录

#读取输入文件
tmb=read.table(tmbFile, header=T, sep="\t", check.names=F, row.names=1)       #读取TMB数据文件
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)     #读取风险文件
clu=read.table(cluFile, header=T, sep="\t", check.names=F, row.names=1)       #读取基因分型的结果文件

#合并数据
tmb=as.matrix(tmb)
tmb=log2(tmb+1)
sameSample=intersect(row.names(tmb), row.names(risk))
tmb=tmb[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]
rownames(clu)=gsub("(.*?)\\_(.*?)", "\\2", rownames(clu))
clu=clu[sameSample,,drop=F]
data=cbind(risk, tmb, clu)
data=data[,c("riskScore", "risk", "geneCluster", "TMB")]

#设置比较组
data$risk=factor(data$risk, levels=c("low", "high"))
risk=levels(factor(data$risk))
comp=combn(risk, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}

#绘制箱线图
boxplot=ggboxplot(data, x="risk", y="TMB", fill="risk",
		          xlab="",
		          ylab="Tumor Burden Mutation",
		          legend.title="Risk",
		          palette=c("#0088FF", "#FF5555") )+ 
	    stat_compare_means(comparisons = my_comparisons)
pdf(file="boxplot.pdf",width=5,height=4.5)
print(boxplot)
dev.off()

#相关性图形
length=length(levels(factor(data$geneCluster)))
bioCol=c("#0066FF","#FF9900","#FF0000","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
p1=ggplot(data, aes(riskScore, TMB)) + 
		  xlab("Risk score")+ylab("Tumor Burden Mutation")+
		  geom_point(aes(colour=geneCluster))+
		  scale_color_manual(values=bioCol[1:length])+ 
		  geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
		  stat_cor(method = 'spearman', aes(x =riskScore, y =TMB))
#相关性图形
pdf(file="cor.pdf", width=6, height=4.6)
print(p1)
dev.off()




