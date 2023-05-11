inputFile="cnvMatrix.txt"     #?????ļ?
setwd("C:\\Users\\ljh\\Desktop\\多组学同死亡\\150.cuproOmics资料\\150.cuproOmics资料0\\14.CNVfreq")      #???ù???Ŀ¼

rt=read.table(inputFile, header=T, sep="\t", check.names=F, row.names=1)    #??ȡ?????ļ?
GAIN=rowSums(rt> 0)         #??ȡ?????????ӵ???Ʒ??Ŀ
LOSS=rowSums(rt< 0)         #??ȡ??????ȱʧ????Ʒ??Ŀ
GAIN=GAIN/ncol(rt)*100      #???㿽???????ӵİٷ???
LOSS=LOSS/ncol(rt)*100      #???㿽????ȱʧ?İٷ???
data=cbind(GAIN, LOSS)
data=data[order(data[,"GAIN"],decreasing = T),]

#????ͼ??
data.max = apply(data, 1, max)
pdf(file="CNVfreq.pdf", width=8, height=5.5)
cex=1.3
par(cex.lab=cex, cex.axis=cex, font.axis=2, las=1, xpd=T)
bar=barplot(data.max, col="grey80", border=NA,
            xlab="", ylab="CNV.frequency(%)", space=1.5,
            xaxt="n", ylim=c(0,1.2*max(data.max)))
points(bar,data[,"GAIN"], pch=20, col=2, cex=3)
points(bar,data[,"LOSS"], pch=20, col=3, cex=3)
legend("top", legend=c('GAIN','LOSS'), col=2:3, pch=20, bty="n", cex=2, ncol=2)
par(srt=45)
text(bar, par('usr')[3]-0.2, rownames(data), adj=1, cex=1)
dev.off()


######??????ѧ??: https://www.biowolf.cn/
######?γ?��??1: https://shop119322454.taobao.com
######?γ?��??2: https://ke.biowolf.cn
######?γ?��??3: https://ke.biowolf.cn/mobile
######?⿡??ʦ????: seqbio@foxmail.com
######?⿡??ʦ΢??: eduBio

