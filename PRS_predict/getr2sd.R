library(fmsb)

liability_transform <-function(r2, k, p) {
    if(r2 == 0) {
        return (r2)
    } else {
        x = qnorm(1 - k)
        z = dnorm(x)
        i = z / k
        cc = k * (1 - k) * k * (1 - k) / (z * z * p * (1 - p))
        theta = i * ((p - k) / (1 - k)) * (i * ((p - k) / (1- k)) - x)
        e = 1 - p^(2 * p) * (1 - p)^(2 * (1 - p))
        r2_L  = cc * e * r2 / (1 + cc * e * theta * r2)
        return (r2_L)
    }
}

## load data
data=read.table("testpheno.txt",header=F,as.is=T)
covar=read.table("testcovar.txt",header=F,as.is=T)
best=read.table("test.prs.best"),header=T,as.is=T)

## prevalences 
k=0.08 ## population prevalence, 0.08 is K for CONVERGE
p=as.numeric(table(data$V3)["1"]/sum(table(data$V3))) ## sample prevalence 

## get PRS r2 
best$phenotyped=data$V3[match(best$FID,data$V1)]
best=merge(best,covar,by="FID",all=F)
best=best[which(best$In_Regression=="Yes"),]
use0=best[,c("phenotyped", paste0("PC", 1:10))]
use1=best[,c("phenotyped", "PRS", paste0("PC", 1:10))]
res0 <- glm(phenotyped ~ ., family=binomial(),data=use0)
res1 <- glm(phenotyped ~ ., family=binomial(),data=use1)
a=NagelkerkeR2(res1)$R2
a_L=liability_transform(a,k,p) 
b=NagelkerkeR2(res0)$R2
b_L=liability_transform(b,k,p)
r2=a-b ## Nagelkerke's R2 
r2_L=a_L-b_L ## Nagelkerke's R2 on liability scale
P=pchisq(deviance(res0) - deviance(res1),df.residual(res0) - df.residual(res1),lower.tail = F) ## P value 

## bootstrap to get SD 
r2_boot=NULL
r2L_boot=NULL
for( i in 1:1e3 ){
    is = sample(1:nrow(use0),replace=T)
    b0 = use0[is,]
    b1 = use1[is,]
    bres0 <- glm(phenotyped ~ ., family=binomial(),data=b0)
    bres1 <- glm(phenotyped ~ ., family=binomial(),data=b1)
    a = NagelkerkeR2(bres1)$R2
    a_L = liability_transform(a,k,p)
    b = NagelkerkeR2(bres0)$R2
    b_L = liability_transform(b,k,p)
    r2_boot[i]=a-b
    r2L_boot[i]=a_L-b_L
}
r2_sd=sd(r2_boot) ## SD for R2
r2L_sd=sd(r2L_boot) ## SD for R2 on liability scale 
