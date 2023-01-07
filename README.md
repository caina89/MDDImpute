# MDDImpute

Code to reproduce analyses in MDD phenotype imputation paper: Dahl A., et al, Phenotype integration improves power and preserves specificity in biobank-based genetic studies of MDD, BioRxiv 2022 (https://www.biorxiv.org/content/10.1101/2022.08.15.503980v1.abstract)

## Phenotype Imputation 

Code used for SoftImpute imputation can be found at https://github.com/andywdahl/mdd-impute 

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

