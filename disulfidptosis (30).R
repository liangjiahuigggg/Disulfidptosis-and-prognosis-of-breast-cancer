#install.packages("ggplot2")
#install.packages("ggpubr")


#???ð?
library(plyr)
library(ggplot2)
library(ggpubr)

riskFile="risk.all.txt"     #?????ļ?
cliFile="MSI.txt"           #΢???ǲ??ȶ????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\49.MSI")     #???ù???Ŀ¼

#??ȡ?????ļ?
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)
cli=read.table(cliFile, header=T, sep="\t", check.names=F, row.names=1)
sameSample=intersect(row.names(risk), row.names(cli))
rt=cbind(risk[sameSample,,drop=F], cli[sameSample,,drop=F])
rt$MSI=factor(rt$MSI, levels=c("MSS", "MSI-L", "MSI-H"))
rt$risk=factor(rt$risk, levels=c("low", "high"))

#????ͼ?ε???ɫ
bioCol=c("#0066FF","#FF0000","#FF9900","#6E568C","#7CC767","#223D6C","#D20A13","#FFD121","#088247","#11AA4D")
bioCol=bioCol[1:length(unique(rt[,"MSI"]))]

#ͳ?Ƹߵͷ????鲡????Ŀ
rt1=rt[,c("MSI", "risk")]
df=as.data.frame(table(rt1))
#?????ߵͷ??????İٷ???
df=ddply(df, .(risk), transform, percent = Freq/sum(Freq) * 100)
#?ٷֱ?λ??
df=ddply(df, .(risk), transform, pos = (cumsum(Freq) - 0.5 * Freq))
df$label=paste0(sprintf("%.0f", df$percent), "%")
#???ưٷֱ???״ͼ
p=ggplot(df, aes(x = factor(risk), y = percent, fill = MSI)) +
	   geom_bar(position = position_stack(), stat = "identity", width = .7) +
	   scale_fill_manual(values=bioCol)+
	   xlab("riskScore")+ ylab("Percent weight")+  guides(fill=guide_legend(title="MSI"))+
	   geom_text(aes(label = label), position = position_stack(vjust = 0.5), size = 3) +
	   #coord_flip()+
	   theme_bw()
pdf(file="barplot.pdf", width=4, height=5)
print(p)
dev.off()

#???ñȽ???
rt2=rt[,c("MSI", "riskScore")]
type=levels(factor(rt2[,"MSI"]))
comp=combn(type, 2)
my_comparisons=list()
for(i in 1:ncol(comp)){my_comparisons[[i]]<-comp[,i]}
#????????ͼ
boxplot=ggboxplot(rt2, x="MSI", y="riskScore", fill="MSI",
		          xlab="",
		          ylab="Risk score",
		          legend.title="MSI",
		          palette=bioCol
		          )+ 
	stat_compare_means(comparisons=my_comparisons)
pdf(file="boxplot.pdf",width=4,height=4.5)
print(boxplot)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

