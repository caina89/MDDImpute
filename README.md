# MDDImpute

Code to reproduce analyses in MDD phenotype imputation paper: Dahl A., et al, Phenotype integration improves power and preserves specificity in biobank-based genetic studies of MDD, BioRxiv 2022 (https://www.biorxiv.org/content/10.1101/2022.08.15.503980v1.abstract)

## Phenotype Imputation 

Code used for SoftImpute imputation can be found at https://github.com/andywdahl/mdd-impute 

## GWAS 

All GWAS in the paper are performed using PLINK2

```
#GWAS on imputed/observed quantitative phenotypes 
plink2 --bfile $bfile --pheno $pfile --linear hide-covar --variance-standardize --covar $covar --out $outfile 
#GWAS on observed binary phenotypes 
plink2 --bfile $bfile --pheno $pfile --1 --logistic hide-covar --covar-variance-standardize --ci 0.95 --covar $covar --out $outfile 
```

## MTAG 

MTAG is a weighted multi-GWAS meta-analysis approach that takes summary statistics from single-trait GWAS and outputs trait-specific association statistics, weighing each of the input GWAS by their genetic correlations with each other. Input file format for MTAG: 

```
snpid	chr	bpos	a1	a2	freq	beta	betase	pval	n
rs58276399	1	731718	C	T	0.1108	0.0232769769044932	0.0150276	0.121372	104776
rs141242758	1	734349	C	T	0.1108	0.02384347168445	0.0150138	0.112334	105072
rs28544273	1	751343	A	T	0.1209	0.0165915950929196	0.0143098	0.24639	107746
rs28527770	1	751756	C	T	0.1211	0.0174567405994606	0.0142902	0.221795	107823
rs3115860	1	753405	C	A	0.1277	0.0269534695602576	0.0138726	0.0520048	109006
```

In this paper we report 6 different MTAG runs; command line for one example (MTAG.GPpsy, where the input GWAS are LifetimeMDD and GPpsy) is shown below 

```
conda activate python2
python mtag.py \
--se_name betase --stream_stdout --n_min 0 --use_beta_se  \
--sumstats LifetimeMDD.mtag.txt,GPpsy.mtag.txt \
--out LifetimeMDD.GPpsy.mtag --std_betas  
``` 

## Summary statistics 

All summary statistics of GWAS described in the paper (SoftImpute, AutoComplete, MTAG, external cohorts) are available on figshare: https://doi.org/10.6084/m9.figshare.19604335.v1 

## PRS Predictions 

PRS predictions are run with PRSice (v2) and we take the best p value thresholds for constructing PRS for each prediction run (floating p value threshold across all runs. 

Summary statistics for all GWAS (available on figshare) have the following format (.ma format also used in SbayesR). All GWAS (in UKBiobank) are performed only on SNPs with INFO score > 0.9 and MAF > 0.05. 

```
SNP	A1	A2	freq	b	se	p	N
rs12726255	G	A	0.139585019432	0.0156010662286512	0.0226	0.4903	42455
rs80057011	C	A	0.093396254858	0.0366951051877892	0.0268	0.1704	42455
rs12063663	C	T	0.092396254858	0.0316970867432546	0.0268	0.2375	42455
rs4081333	T	C	0.0926773784007	0.0295976364389436	0.0267	0.2681	42455
rs12064046	G	T	0.091	0.0307994715535284	0.0268	0.2507	42455
```

We use the following PRSice command line, switching --binary-target T/F for binary and quantitative target phenotypes. We performed PRS in 10-fold cross validation in UKBiobank for in-sample PRS analysis, so phenotype/genotype files are in folds, to do this without having separate files one can use --remove/--keep options 

```
a=$GWAS ## summary statistics for SoftImpute, Autocomplete or MTAG GWAS 
b=$targetpheno ## 217 target phenotypes in UKBiobank as detailed in Supplementary Table 1 of the paper 
for test in {1..10}
do 
PRSice_linux --base $a.ma --cov $covar --num-auto 22 \
--a1 A1 --a2 A2 --pvalue p --snp SNP --stat b --beta --base-maf freq:0.05 \
--clump-kb 250kb --clump-p 1.000000 --clump-r2 0.100000 \
--interval 5e-05 --upper 0.5 --lower 5e-08 --num-auto 22 \
--bar-levels 0.001,0.05,0.1,0.2,0.3,0.4,0.5,1 \
--out $outdir/$b/$a.test$test.prs \
--pheno $phendir/test$test.$b.phen \
--target $genodir/tenfold/allchr.test$test --binary-target F \
--thread 1
done 
```



