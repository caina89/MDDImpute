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

mdd=data[which(data$pheno=="LifetimeMDD"),]
data$prs.test=paste(data$prs,data$test,sep=":")
mdd$prs.test=paste(mdd$prs,mdd$test,sep=":")
data$MDD.R2=mdd$PRS.R2[match(data$prs.test,mdd$prs.test)]
data$ratio=data$PRS.R2/data$MDD.R2
write.table(data,paste0(wdir,"/Downsampling.PRSPleio.txt"),sep="\t",quote=F,row.names=T,col.names=T)

library(plyr)
dsummary=ddply(data,.(PRSPheno,SampleSize,pheno),summarise,meanp=mean(P),meanratio=mean(ratio),sdratio=sd(ratio),meanr2=mean(PRS.R2),cvr2=sd(PRS.R2/mean(PRS.R2)),meanthresh=mean(Threshold),meanlogp=mean(-log10(Threshold)), sdthresh=sd(Threshold),sdlogp=sd(-log10(Threshold)), sdsnp=sd(Num_SNP),meansnp=mean(Num_SNP))
dsummary$lower=dsummary$meanratio-dsummary$sdratio*1.96
dsummary$upper=dsummary$meanratio+dsummary$sdratio*1.96
write.table(dsummary,paste0(wdir,"/Downsampling.PRSPleio.summary.txt"),sep="\t",quote=F,row.names=T,col.names=T)
