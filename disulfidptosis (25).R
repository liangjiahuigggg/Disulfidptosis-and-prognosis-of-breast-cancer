#if (!requireNamespace("BiocManager", quietly = TRUE))
#    install.packages("BiocManager")
#BiocManager::install("limma")

#install.packages("tidyverse")
#install.packages("ggplot2")
#install.packages("ggExtra")
#install.packages("ggpubr")


#???ð?
library(limma)
library(reshape2)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggExtra)

immFile="CIBERSORT-Results.txt"     #????ϸ???????Ľ????ļ?
riskFile="risk.all.txt"             #?????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\43.immuneCor")     #???ù???Ŀ¼

#??ȡ????ϸ???????Ľ????ļ??????????ݽ???????
immune=read.table(immFile, header=T, sep="\t", check.names=F, row.names=1)
immune=immune[immune[,"P-value"]<0.05,]
data=as.matrix(immune[,1:(ncol(immune)-3)])
rownames(data)=gsub("(.*?)\\_(.*?)", "\\2", rownames(data))

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(data), row.names(risk))
data=data[sameSample,,drop=F]
risk=risk[sameSample,,drop=F]

#??????????ϸ??????ѭ?????õ????մ?????????ϸ??????????
for(i in colnames(data)[1:ncol(data)]){
	x=as.numeric(risk[,"riskScore"])
	x[x>quantile(x,0.99)]=quantile(x,0.99)
	y=as.numeric(data[,i])
	if(sd(y)<0.01){next}
	cor=cor.test(x, y, method="spearman")
	#??pvalueС??0.05??????ϸ????????????ɢ??ͼ?Ļ???
	if(cor$p.value<0.05){
		outFile=paste0("cor.", i, ".pdf")
		df1=as.data.frame(cbind(x,y))
		p1=ggplot(df1, aes(x, y)) + 
				  xlab("Risk score") + ylab(i)+
				  geom_point() + geom_smooth(method="lm",formula = y ~ x) + theme_bw()+
				  stat_cor(method = 'spearman', aes(x =x, y =y))
		p2=ggMarginal(p1, type="density", xparams=list(fill = "orange"), yparams=list(fill = "blue"))
		#??????ͼ??
		pdf(file=outFile, width=5.2, height=5)
		print(p2)
		dev.off()
	}
}

#??????????ϸ???????Է???
outTab=data.frame()
risk=risk[,3:(ncol(risk)-1),drop=F]
for(immune in colnames(data)){
	if(sd(data[,immune])<0.01){next}
	for(gene in colnames(risk)){
		x=as.numeric(data[,immune])
		y=as.numeric(risk[,gene])
		y[y>quantile(y,0.99)]=quantile(y,0.99)
		corT=cor.test(x,y,method="spearman")
		cor=corT$estimate
		pvalue=corT$p.value
		text=ifelse(pvalue<0.001,"***",ifelse(pvalue<0.01,"**",ifelse(pvalue<0.05,"*","")))
		outTab=rbind(outTab,cbind(Gene=gene, Immune=immune, cor, text, pvalue))
	}
}

#????????????ͼ
outTab$cor=as.numeric(outTab$cor)
outTab[,"Gene"]=factor(outTab[,"Gene"], levels=colnames(risk))
pdf(file="geneImmuneCor.pdf", width=7, height=5.5)
ggplot(outTab, aes(Gene, Immune)) + 
	geom_tile(aes(fill = cor), colour = "grey", size = 1)+
	scale_fill_gradient2(low = "#5C5DAF", mid = "white", high = "#EA2E2D") + 
	geom_text(aes(label=text),col ="black",size = 3) +
	theme_minimal() + 
	theme(axis.title.x=element_blank(), axis.ticks.x=element_blank(), axis.title.y=element_blank(),
	      axis.text.x = element_text(angle = 45, hjust = 1, size = 9, face = "bold"),   #x??????
	      axis.text.y = element_text(size = 9, face = "bold")) +       #y??????
	labs(fill =paste0("***  p<0.001","\n", "**  p<0.01","\n", " *  p<0.05","\n", "\n","Correlation")) +   #????ͼ??
	scale_x_discrete(position = "bottom")      #X????????????ʾ??ͼ?ε??·?
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

