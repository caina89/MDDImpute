library(ggplot2)
library(plyr)
load("utilities.R")
cols=c("#f4a582","#ca0020","#404040")

## get data 
phenos=read.table("sigphenos.txt")
phenos=c(as.character(phenos$x),"LifetimeMDD")
tests=c(1:10)
data=NULL
prss=c("downsample.GPpsy.50k", "downsample.GPpsy.100k", "downsample.SoftImpAll.50k","downsample.SoftImpAll.100k", "GPpsy", "SoftImpAll")
for (pheno in phenos){
	print(pheno)
	for (prs in prss){
		for (test in tests){
			prsfile=paste0(wdir,"/", pheno,"/",prs,".test",test,".prs.summary")
			if (file.exists(prsfile)){
				d=tryElse(read.table(paste0(wdir,"/", pheno,"/",prs,".test",test,".prs.summary"),header=T,as.is=T),data.frame(matrix(NA, nrow = 1, ncol = ncol(data), dimnames=list(NULL, colnames(data)))))
				d$pheno=pheno
				d$prs=prs
				d$test=test
				data=rbind(data,d)
			}
		}
	}
}
data$PRSPheno="GPpsy"
data$PRSPheno[grep("SoftImpAll",data$prs)]="SoftImpAll"
data$SampleSize="Full"
data$SampleSize[grep("50k",data$prs)]="DS_50K"
data$SampleSize[grep("100k",data$prs)]="DS_100K"

## get PRS Pleiotropy per fold 
mdd=data[which(data$pheno=="LifetimeMDD"),]
data$prs.test=paste(data$prs,data$test,sep=":")
mdd$prs.test=paste(mdd$prs,mdd$test,sep=":")
data$MDD.R2=mdd$PRS.R2[match(data$prs.test,mdd$prs.test)]
data$ratio=data$PRS.R2/data$MDD.R2
write.table(data,paste0(wdir,"/Downsampling.PRSPleio.txt"),sep="\t",quote=F,row.names=T,col.names=T)

## get PRS Pleiotropy across folds 
dsummary=ddply(data,.(PRSPheno,SampleSize,pheno),summarise,meanp=mean(P),meanratio=mean(ratio),sdratio=sd(ratio),meanr2=mean(PRS.R2),cvr2=sd(PRS.R2/mean(PRS.R2)),meanthresh=mean(Threshold),meanlogp=mean(-log10(Threshold)), sdthresh=sd(Threshold),sdlogp=sd(-log10(Threshold)), sdsnp=sd(Num_SNP),meansnp=mean(Num_SNP))
dsummary$lower=dsummary$meanratio-dsummary$sdratio*1.96
dsummary$upper=dsummary$meanratio+dsummary$sdratio*1.96
write.table(dsummary,paste0(wdir,"/Downsampling.PRSPleio.summary.txt"),sep="\t",quote=F,row.names=T,col.names=T)

## testing and plotting 

## GPpsy
gppsy=dsummary[which(dsummary$PRSPheno=="GPpsy"&dsummary$SampleSize=="Full"),]
gppsy=gppsy[order(gppsy$meanratio,decreasing=T),] ## order by PRS Pleitropy
gppsy=gppsy[which(gppsy$meanp<0.05),] ## only keep phenos with mean PRS prediction P value < 0.05 
gppsy$orderpleio=seq(1:nrow(gppsy))
gpheno=gppsy$pheno

## for 50K downsampling 
g1=gppsy[,c("pheno","orderpleio","meanratio")]
g50=dsummary[which(dsummary$PRSPheno=="GPpsy"&dsummary$SampleSize=="DS_50K"),]
g50=g50[match(gpheno,g50$pheno),]
g1$dsmeanratio=g50$meanratio
g1$SampleSize="DS_50K"

## test for correlation between full and downsampled (50K) PRS pleiotropy measures for GPpsy 
cor.test(g1$meanratio,g1$dsmeanratio)
#         Pearson's product-moment correlation
# data:  g1$meanratio and g1$dsmeanratio
# t = 24.472, df = 61, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9226803 0.9711893
# sample estimates:
#       cor 
# 0.9526597 

## test for difference between full and downsampled (50K) PRS pleiotropy measures for GPpsy 
t.test(g1$meanratio,g1$dsmeanratio,paired=T,alternative="two.sided") 
#         Paired t-test
# data:  g1$meanratio and g1$dsmeanratio
# t = -8.0562, df = 62, p-value = 3.162e-11
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.12784747 -0.07701512
# sample estimates:
# mean of the differences 
#              -0.1024313 

## for 100K downsampling 
g2=gppsy[,c("pheno","orderpleio","meanratio")]
g100=dsummary[which(dsummary$PRSPheno=="GPpsy"&dsummary$SampleSize=="DS_100K"),]
g100=g100[match(gpheno,g100$pheno),]
g2$dsmeanratio=g100$meanratio
g2$SampleSize="DS_100K"

## test for correlation between full and downsampled (100K) PRS pleiotropy measures for GPpsy 
cor.test(g2$meanratio,g2$dsmeanratio)
#         Pearson's product-moment correlation
# data:  g2$meanratio and g2$dsmeanratio
# t = 53.558, df = 61, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9826987 0.9936771
# sample estimates:
#       cor 
# 0.9895337 

## test for difference between full and downsampled (100K) PRS pleiotropy measures for GPpsy 
t.test(g2$meanratio,g2$dsmeanratio,paired=T,alternative="two.sided")
#         Paired t-test
# data:  g2$meanratio and g2$dsmeanratio
# t = 0.35008, df = 62, p-value = 0.7275
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.006766895  0.009640266
# sample estimates:
# mean of the differences 
#             0.001436685

## prepare plotting data 
forplot1=rbind(g1,g2)
forplot1$SampleSize=factor(forplot1$SampleSize, levels=c("DS_50K","DS_100K"))
forplot2=dsummary[which(dsummary$PRSPheno=="GPpsy"),c("meanratio","SampleSize","pheno")]
forplot2=forplot2[which(forplot2$pheno%in%gpheno),]
forplot2$order=gppsy$orderpleio[match(forplot2$pheno,gppsy$pheno)]
forplot2$SampleSize=factor(forplot2$SampleSize, levels=c("DS_50K","DS_100K","Full"))

## SoftImpAll
soft=dsummary[which(dsummary$PRSPheno=="SoftImpAll"&dsummary$SampleSize=="Full"),]
soft=soft[order(soft$meanratio,decreasing=T),]
soft=soft[which(soft$meanp<0.05),]
soft$orderpleio=seq(1:nrow(soft))
spheno=soft$pheno

## for 50K downsampling 
s1=soft[,c("pheno","orderpleio","meanratio")]
s50=dsummary[which(dsummary$PRSPheno=="SoftImpAll"&dsummary$SampleSize=="DS_50K"),]
s50=s50[match(spheno,s50$pheno),]
s1$dsmeanratio=s50$meanratio
s1$SampleSize="DS_50K"

## test for difference between full and downsampled (50K) PRS pleiotropy measures for SoftImpAll 
cor.test(s1$meanratio,s1$dsmeanratio)
#         Pearson's product-moment correlation
# data:  s1$meanratio and s1$dsmeanratio
# t = 34.353, df = 61, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9590592 0.9849237
# sample estimates:
#       cor 
# 0.9751159

## test for difference between full and downsampled (50K) PRS pleiotropy measures for SoftImpAll 
t.test(s1$meanratio,s1$dsmeanratio,paired=T,alternative="two.sided")
#         Paired t-test
# data:  s1$meanratio and s1$dsmeanratio
# t = -0.1277, df = 62, p-value = 0.8988
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.011030573  0.009705839
# sample estimates:
# mean of the differences 
#           -0.0006623668

## 100K downsampling 
s2=soft[,c("pheno","orderpleio","meanratio")]
s100=dsummary[which(dsummary$PRSPheno=="SoftImpAll"&dsummary$SampleSize=="DS_100K"),]
s100=s100[match(spheno,s100$pheno),]
s2$dsmeanratio=s100$meanratio
s2$SampleSize="DS_100K"

## test for correlation between full and downsampled (100K) PRS pleiotropy measures for SoftImpAll 
cor.test(s2$meanratio,s2$dsmeanratio)
#         Pearson's product-moment correlation
# data:  s2$meanratio and s2$dsmeanratio
# t = 67.789, df = 61, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.9891225 0.9960329
# sample estimates:
#       cor 
# 0.9934281

## test for difference between full and downsampled (100K) PRS pleiotropy measures for SoftImpAll 
t.test(s2$meanratio,s2$dsmeanratio,paired=T,alternative="two.sided")
#         Paired t-test
# data:  s2$meanratio and s2$dsmeanratio
# t = 0.3896, df = 62, p-value = 0.6982
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#  -0.004183885  0.006209598
# sample estimates:
# mean of the differences 
#             0.001012857 

## prepare plotting data 
forplot3=rbind(s1,s2)
forplot3$SampleSize=factor(forplot3$SampleSize, levels=c("DS_50K","DS_100K"))
forplot4=dsummary[which(dsummary$PRSPheno=="SoftImpAll"),c("meanratio","SampleSize","pheno")]
forplot4=forplot4[which(forplot4$pheno%in%spheno),]
forplot4$order=soft$orderpleio[match(forplot4$pheno,soft$pheno)]
forplot4$SampleSize=factor(forplot4$SampleSize, levels=c("DS_50K","DS_100K","Full"))

## plot (Extended Data Figure 5)
pdf("Downsample.PRSPleio.pdf",height=6,width=8,useDingbats=FALSE)
p1=ggplot(forplot1,aes(meanratio,dsmeanratio,colour=SampleSize))+geom_point()+geom_smooth(method="lm")+theme_bw()
p1=p1+ylim(0,1)+xlim(0,1)+geom_abline(intercept = 0, slope = 1)+ggtitle("GPpsy")+scale_colour_manual(values=cols[1:2])
p1=p1+xlab("FullSample/nPRS Pleiotropy") + ylab("Downsampled\nPRS Pleiotropy")
p1=p1+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p2=ggplot(forplot2,aes(order,meanratio,colour=SampleSize))+geom_point()+geom_line()+theme_bw()
p2=p2+scale_colour_manual(values=cols)+ylim(0,1)+xlim(0,65)+ggtitle("GPpsy")
p2=p2+xlab("Phenotype") + ylab("PRS Pleiotropy")
p2=p2+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p3=ggplot(forplot3,aes(meanratio,dsmeanratio,colour=SampleSize))+geom_point()+geom_smooth(method="lm")+theme_bw()
p3=p3+ylim(0,1)+xlim(0,1)+geom_abline(intercept = 0, slope = 1)+ggtitle("SoftImpAll")+scale_colour_manual(values=cols[1:2])
p3=p3+xlab("FullSample/nPRS Pleiotropy") + ylab("Downsampled\nPRS Pleiotropy")
p3=p3+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p4=ggplot(forplot4,aes(order,meanratio,colour=SampleSize))+geom_point()+geom_line()+theme_bw()
p4=p4+scale_colour_manual(values=cols)+ylim(0,1)+xlim(0,65)+ggtitle("SoftImpAll")
p4=p4+xlab("Phenotype") + ylab("PRS Pleiotropy")
p4=p4+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
multiplot(p1,p3,p2,p4,cols=2)
dev.off()

## GPpsy: do sample size affect p value thresholds  
forplot1=dsummary[which(dsummary$PRSPheno=="GPpsy"),]
forplot1=forplot1[which(forplot1$pheno%in%gpheno),]
forplot1$order=gppsy$orderpleio[match(forplot1$pheno,gppsy$pheno)]
forplot1$SampleSize=factor(forplot1$SampleSize, levels=c("DS_50K","DS_100K","Full"))

## plot (Supplementary Figure 10)
pdf("GPpsy.Pthreshold.pdf",height=8,width=9,useDingbats=FALSE)
p1=ggplot(forplot1,aes(pheno,meanthresh,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p1=p1+theme_bw()+theme(axis.text.x=element_blank())
p1=p1+scale_colour_manual(values=cols)
p1=p1+xlab("Phenotype") + ylab("Mean (P Threshold)")
p1=p1+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p3=ggplot(forplot1,aes(pheno,sdthresh,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p3=p3+theme_bw()+theme(axis.text.x=element_blank())
p3=p3+scale_colour_manual(values=cols)
p3=p3+xlab("Phenotype") + ylab("SD (P Threshold)")
p3=p3+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p5=ggplot(forplot1,aes(pheno,meansnp,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p5=p5+theme_bw()+theme(axis.text.x=element_blank())
p5=p5+scale_colour_manual(values=cols)
p5=p5+xlab("Phenotype") + ylab("Mean (NSNPs)")
p5=p5+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p6=ggplot(forplot1,aes(pheno,sdsnp,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p6=p6+theme_bw()+theme(axis.text.x=element_blank())
p6=p6+scale_colour_manual(values=cols)
p6=p6+xlab("Phenotype") + ylab("SD (NSNPs)")
p6=p6+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p7=ggplot(forplot1,aes(pheno,cvr2,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p7=p7+theme_bw()+theme(axis.text.x=element_blank())
p7=p7+scale_colour_manual(values=cols)
p7=p7+xlab("Phenotype") + ylab("SD (R2/mean(R2))")
p7=p7+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p8=ggplot(forplot1,aes(pheno,sdratio,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p8=p8+theme_bw()+theme(axis.text.x=element_blank())
p8=p8+scale_colour_manual(values=cols)
p8=p8+xlab("Phenotype") + ylab("SD (PRS Pleiotropy)")
p8=p8+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
multiplot(p1,p5,p7,p3,p6,p8,cols=2)
dev.off()

## SoftImpAll: do sample size affect p value thresholds 
forplot1=dsummary[which(dsummary$PRSPheno=="SoftImpAll"),]
forplot1=forplot1[which(forplot1$pheno%in%spheno),]
forplot1$order=soft$orderpleio[match(forplot1$pheno,soft$pheno)]
forplot1$SampleSize=factor(forplot1$SampleSize, levels=c("DS_50K","DS_100K","Full"))

## plot (Supplementary Figure 11)
pdf("SoftImpAll.Pthreshold.pdf",height=8,width=9,useDingbats=FALSE)
p1=ggplot(forplot1,aes(pheno,meanthresh,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p1=p1+theme_bw()+theme(axis.text.x=element_blank())
p1=p1+scale_colour_manual(values=cols)
p1=p1+xlab("Phenotype") + ylab("Mean (P Threshold)")
p1=p1+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p3=ggplot(forplot1,aes(pheno,sdthresh,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p3=p3+theme_bw()+theme(axis.text.x=element_blank())
p3=p3+scale_colour_manual(values=cols)
p3=p3+xlab("Phenotype") + ylab("SD (P Threshold)")
p3=p3+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p5=ggplot(forplot1,aes(pheno,meansnp,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p5=p5+theme_bw()+theme(axis.text.x=element_blank())
p5=p5+scale_colour_manual(values=cols)
p5=p5+xlab("Phenotype") + ylab("Mean (NSNPs)")
p5=p5+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p6=ggplot(forplot1,aes(pheno,sdsnp,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p6=p6+theme_bw()+theme(axis.text.x=element_blank())
p6=p6+scale_colour_manual(values=cols)
p6=p6+xlab("Phenotype") + ylab("SD (NSNPs)")
p6=p6+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p7=ggplot(forplot1,aes(pheno,cvr2,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p7=p7+theme_bw()+theme(axis.text.x=element_blank())
p7=p7+scale_colour_manual(values=cols)
p7=p7+xlab("Phenotype") + ylab("SD (R2/mean(R2))")
p7=p7+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
p8=ggplot(forplot1,aes(pheno,sdratio,group=SampleSize,colour=SampleSize))+geom_point()+geom_smooth()
p8=p8+theme_bw()+theme(axis.text.x=element_blank())
p8=p8+scale_colour_manual(values=cols)
p8=p8+xlab("Phenotype") + ylab("SD (PRS Pleiotropy)")
p8=p8+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
multiplot(p1,p5,p7,p3,p6,p8,cols=2)
dev.off()







