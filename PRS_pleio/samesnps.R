library(ggplot2)
library(plyr)
load("utilities.R")
cols=c( 'orange3', 'orange1', 'grey', "black", 'blue1', 6, 'blue4', 'cornflowerblue', 'cyan', 'purple3', 'purple1', 'pink1', 'pink4' )

## wdir = prs output directory 

quantphenos=read.table("pheno.quant.txt",as.is=T,header=F)
binaryphenos=read.table("pheno.binary.txt",as.is=T,header=F)
phenos=c(quantphenos$V1,binaryphenos$V1)

## PRS Pleiotropy calculated for best p value thresholds  

adata=read.table("bestp.PRSPleio.raw.txt"),header=T,as.is=T)
asummary=ddply(adata,.(pheno,prs),summarise,meanp=mean(P),meanratio=mean(ratio),sdratio=sd(ratio))

## PRS Pleiotropy using PGC29 SNPs (SNPs clumped for prediction of LifetimeMDD at P threshold 0.5) ld LD pruned SNPs (r2<0.2) 

snps=c("pgc29SNPs", "ld2") 

for (snp in snps){ 
    ## get data 
    tests=c(1:10)
    data=NULL
    prss=c("Soft-ImpAll","Soft-ImpOnly","Auto-ImpAll","Auto-ImpOnly","MTAG.FamilyHistory","MTAG.Envs","MTAG.All","GPpsy","LifetimeMDD")
    for (pheno in phenos){
      print(pheno)
      for (prs in prss){
        for (test in tests){
          prsfile=paste0(wdir,"/", pheno,"/",prs,".test",test,".",snp,".prs.summary")
          if (file.exists(prsfile)){
            d=tryElse(read.table(paste0(wdir,"/", pheno,"/",prs,".test",test,".",snp,".prs.summary"),header=T,as.is=T),data.frame(matrix(NA, nrow = 1, ncol = ncol(data), dimnames=list(NULL, colnames(data)))))
            d$pheno=pheno
            d$prs=prs
            d$test=test
            data=rbind(data,d)
          }
        }
      }
    }
    mdd=data[which(data$pheno=="LifetimeMDD"),]
    data$prs.test=paste(data$prs,data$test,sep=":")
    mdd$prs.test=paste(mdd$prs,mdd$test,sep=":")
    data$MDD.R2=mdd$PRS.R2[match(data$prs.test,mdd$prs.test)]
    data$ratio=data$PRS.R2/data$MDD.R2
    write.table(data,paste0(wdir,"/",snp,".PRSPleio.txt"),sep="\t",quote=F,row.names=T,col.names=T)
    dsummary=ddply(pgcdata,.(pheno,prs),summarise,meanp=mean(P),meanratio=mean(ratio),sdratio=sd(ratio),meanr2=mean(PRS.R2),cvr2=sd(PRS.R2/mean(PRS.R2)),meanthresh=mean(Threshold),meanlogp=mean(-log10(Threshold)), sdthresh=sd(Threshold),sdlogp=sd(-log10(Threshold)), sdsnp=sd(Num_SNP),meansnp=mean(Num_SNP))
    dsummary$lower=dsummary$meanratio-dsummary$sdratio*1.96
    dsummary$upper=dsummary$meanratio+dsummary$sdratio*1.96
    # write.table(dsummary,paste0(wdir,"/SoftImpAll.downsampled.txt"),sep="\t",quote=F,row.names=T,col.names=T)
    write.table(dsummary,paste0(wdir,"/",snp,".PRSPleio.summary.txt"),sep="\t",quote=F,row.names=T,col.names=T)

    ## get phenotypes that have mean PRS prediction P value < 0.05 at best p value thresholds  
    prs=unique(dsummary$prs)
    odata=NULL
    for (p in prs){
      pdata=dsummary[which(dsummary$prs==p),]
      fdata=asummary[which(asummary$prs==p),]
      puse=fdata$pheno[which(fdata$meanp<0.05)]
      pdata=pdata[which(pdata$pheno%in%puse),]
      pdata=pdata[order(pdata$meanratio,decreasing=T),]
      pdata$order=seq(1,nrow(pdata))
      odata=rbind(odata,pdata)
    }
    idata=odata[which(odata$prs %in% c("LifetimeMDD", "GPpsy", "Soft-ImpAll", "Soft-ImpOnly", "Auto-ImpAll", "Auto-ImpOnly")),]
    mdata=odata[which(odata$prs %in% c("LifetimeMDD", "GPpsy", "MTAG.Envs", "MTAG.FamilyHistory", "MTAG.All")),]
    mdata$prs[which(mdata$prs=="MTAG.FamilyHistory")]="MTAG.FamHist" ## shorten for plotting
    idata$prs=factor(idata$prs, levels=c("LifetimeMDD", "GPpsy", "Auto-ImpAll", "Auto-ImpOnly", "Soft-ImpAll", "Soft-ImpOnly"))
    mdata$prs=factor(mdata$prs, levels=c("LifetimeMDD", "GPpsy", "MTAG.Envs", "MTAG.All", "MTAG.FamilyHistory"))

    ## plot
    icols=cols[c( 4, 3, 1, 2, 10, 11 )]
    mcols=cols[c( 4, 3, 7, 9, 8 )]
    pdf(paste0(snp,".PRSPleio.pdf"),height=3,width=8,useDingbats=FALSE)
    p2=ggplot(idata,aes(order,meanratio,colour=prs))+geom_line()+theme_bw()
    p2=p2+scale_colour_manual(values=icols)+ylim(0,1)+xlim(0,100)+ggtitle("PGC29 SNPs (P < 0.5)")
    p2=p2+xlab("Phenotype") + ylab("PRS Pleiotropy")
    p2=p2+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
    p4=ggplot(mdata,aes(order,meanratio,colour=prs))+geom_line()+theme_bw()
    p4=p4+scale_colour_manual(values=mcols)+ylim(0,1)+xlim(0,100)+ggtitle("PGC29 SNPs (P < 0.5)")
    p4=p4+xlab("Phenotype") + ylab("PRS Pleiotropy")
    p4=p4+theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
    panel.background = element_blank(), axis.line = element_line(colour = "black"), legend.position="right")
    multiplot(p2,p4,cols=2)
    dev.off()
}


