#install.packages("ggplot2")
#install.packages("ggalluvial")


#???ð?
library(ggalluvial)
library(ggplot2)
library(dplyr)

crgCluFile="CRGcluster.txt"       #ͭ???????͵Ľ????ļ?
geneCluFile="geneCluster.txt"     #???????͵Ľ????ļ?
riskFile="risk.all.txt"           #?????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\35.ggalluvial")     #???ù???Ŀ¼

#??ȡ?????ļ?
crgClu=read.table(crgCluFile, header=T, sep="\t", check.names=F, row.names=1)
geneClu=read.table(geneCluFile, header=T, sep="\t", check.names=F, row.names=1)
risk=read.table(riskFile, header=T, sep="\t", check.names=F, row.names=1)

#?ϲ?????
twoCluster=cbind(crgClu, geneClu)
rownames(twoCluster)=gsub("(.*?)\\_(.*?)", "\\2", rownames(twoCluster))
sameSample=intersect(row.names(twoCluster), row.names(risk))
rt=cbind(risk[sameSample,,drop=F], twoCluster[sameSample,,drop=F])

#׼??ɣ??ͼ?????ļ?
rt=rt[,c("CRGcluster", "geneCluster", "risk", "fustat")]
colnames(rt)=c("CRGcluster", "geneCluster", "Risk", "Fustat")
rt[,"Fustat"]=ifelse(rt[,"Fustat"]==0, "Alive", "Dead")
corLodes=to_lodes_form(rt, axes = 1:ncol(rt), id = "Cohort")

#????ɣ??ͼ
pdf(file="ggalluvial.pdf", width=6, height=5.5)
mycol=rep(c("#0066FF","#FF9900","#FF0000","#029149","#6E568C","#E0367A","#D8D155","#223D6C","#D20A13","#431A3D","#91612D","#FFD121","#088247","#11AA4D","#58CDD9","#7A142C","#5D90BA","#64495D","#7CC767"),15)
ggplot(corLodes, aes(x = x, stratum = stratum, alluvium = Cohort,fill = stratum, label = stratum)) +
  	 scale_x_discrete(expand = c(0, 0)) +  
  	 #??????????ɫ??forward˵????????ɫ??ǰ??????״ͼһ?£?backward˵????????ɫ??????????״ͼһ?¡?
  	 geom_flow(width = 2/10,aes.flow = "forward") + 
	 geom_stratum(alpha = .9,width = 2/10) +
	 scale_fill_manual(values = mycol) +
	 #size=3??????????С
	 geom_text(stat = "stratum", size = 3,color="black") +
	 xlab("") + ylab("") + theme_bw() + 
	 theme(axis.line = element_blank(),axis.ticks = element_blank(),axis.text.y = element_blank()) + #ȥ????????
	 theme(panel.grid =element_blank()) + 
	 theme(panel.border = element_blank()) + 
	 ggtitle("") + guides(fill = FALSE)                            
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

