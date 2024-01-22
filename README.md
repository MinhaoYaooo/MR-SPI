# Introduction

`MR.SPI` is an R package to conduct two-sample Mendelian randomization by first selecting valid IVs and then perform post-selection inference. The inputs of `MR.SPI` are summary statistics which can be obtained in Genome-Wide Association Studies (GWAS). The outputs of `MR.SPI` are the estimation of causal effect and its corresponding confidence interval, and the set of valid IVs.

# Installation
You can install the development version of `MR.SPI` from Github via the `devtools` package.
```
devtools::install_github("MinhaoYaooo/MR-SPI")
```

# Simulation Example

We first set the following parameters:

```
n1 <- 10000  # number of sample 1
n2 <- 20000 # number of sample 2
beta <- 0.5 # magnitude of causal effect
tau <- 0.25 # IV strength
gamma <- tau*c(1,1,1,1,1,-1,-1,-1,-1,-1)  # SNP-exposure effect
pi <- tau*c(0,0,0,0,0,0,1,1,1,1) # invalidity of the candidate IVs
```

Since `MR.SPI` is used for two-sample Mendelian randomization, we create two independent samples whose genotypes follow the same distribution:

```
pSNP <- runif(10, 0.05, 0.95)  # generate allele frequency

Z1 <- sapply(pSNP, rbinom, n = n1, size = 2) # generate raw genotype of sample 1
Z1 <- scale(Z1)

Z2 <- sapply(pSNP, rbinom, n = n2, size = 2) # generate raw genotype of sample 2
Z2 <- scale(Z2)

D1 <- Z1 %*% gamma + rnorm(n1, 0, 1) # generate the exposure of sample 1

epsilonSigma = matrix(c(1,0.8,0.8,1),2,2)
epsilon = MASS::mvrnorm(n2,rep(0,2),epsilonSigma)
D2 <- Z2 %*% gamma + epsilon[,1] # generate the exposure of sample 2
Y2 <- D2*beta + Z2 %*% pi + epsilon[,2] # generate the outcome of sample 2
```

After generating two samples and their corresponding phenotypes, we calculate the summary statistics as the inputs:

```
library(bigstatsr)

GWAS1 <- big_univLinReg(X = as_FBM(Z1), y.train = D1) # calculate summary statistics of the exposure in sample 1.
GWAS2 <- big_univLinReg(X = as_FBM(Z2), y.train = Y2) # calculate summary statistics of the outcome in sample 2.
```

Then with the summary statistics, we can estimate the causal effect with `MR.SPI` and the corresponding standard confidence interval:

```
library(igraph)

GammaHat = as.numeric(GWAS2$estim); gammaHat = as.numeric(GWAS1$estim) 
se_Gamma <- as.numeric(GWAS2$std.err); se_gamma <- as.numeric(GWAS1$std.err);

mr.spi.standard <- MR.SPI(gammaHat, GammaHat, se_gamma, se_Gamma, n1, n2)
```

If we wish to construct a confidence interval that is robust to finite-sample IV selection error, then we set `robust=TRUE`:

```
library(intervals)

mr.spi.robust <-MR.SPI(gammaHat, GammaHat, se_gamma, se_Gamma, n1, n2, robust = TRUE)
```

# Real Data Example

In this example, we analyze the causal effect of TREM2 on Alzheimer's disease (AD). The exposure data is extracted from the proteomics data in UKB-PPP, which can be downloaded [HERE](https://www.biorxiv.org/content/10.1101/2022.06.17.496443v1.supplementary-material):

```
library(ieugwasr)
library(TwoSampleMR)
library(MR.SPI)
library(dplyr)
library(stringr)
library(data.table)

protein.dat <- openxlsx::read.xlsx("SNPs.xlsx", sheet = 1)                                       # read the proteomics data
protein.list <- unique(protein.dat$Assay.Target)                                                 # filter the unique proteins
info <- str_split_fixed(protein.dat$`Variant.ID.(CHROM:GENPOS.(hg37):A0:A1:imp:v1)`,":",6)       # extract the position information
protein.dat.process <- protein.dat[,c(2,3,9,11,12,13,14,15)]                                     # extract useful columns
protein.dat.process$A1 <- info[,4]; protein.dat.process$A2 <- info[,3]                           # add allele information
colnames(protein.dat.process) <- c("chr","pos","id.exposure","SNP","eaf.exposure","beta.exposure","se.exposure","pval.exposure",
                                   "effect_allele.exposure","other_allele.exposure")             # rename the columns

exp_dat <- protein.dat.process[protein.dat.process$id.exposure=='TREM2',]                        # extract the summary statistics for TREM2
```
Then, we extract the outcome data for AD, which can be downloaded [HERE](https://ctg.cncr.nl/software/summary_statistics):
```
AD <- fread('AD_sumstats_Jansenetal_2019sept.txt')
AD <- format_data(AD, type = 'outcome', 
                  snp_col = 'SNP',
                  beta_col = 'BETA',
                  se_col = 'SE',
                  eaf_col = 'EAF',
                  effect_allele_col = 'A1',
                  other_allele_col = 'A2',
                  pval_col = 'P',
                  samplesize_col = 'Nsum', 
                  chr_col = 'CHR',
                  pos_col = 'BP')

out_dat <- AD[AD$SNP %in% exp_dat$SNP,]
dat <- harmonise_data(exp_dat, out_dat, action=1)
```

