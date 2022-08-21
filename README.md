# MR-SPI

`MR.SPI` is an R package to conduct two-sample Mendelian randomization by first selecting valid IVs and then perform post-selection inference. The inputs of `MR.SPI` are summary statistics which can be obtained in Genome-Wide Association Studies (GWAS). The outputs of `MR.SPI` are the estimation of causal effect and its corresponding confidence interval, and a set of valid IVs <img src="https://render.githubusercontent.com/render/math?math=\hat{V}">.

# Installation
You can install the development version of `MR.SPI` from Github via the `devtools` package.
```
devtools::install_github("MinhaoYaooo/MR-SPI")
```

# Example

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

GWAS1 <- big_univLinReg(X = as_FBM(Z1), y.train = D1) # calculate summary statistics of D~Z
GWAS2 <- big_univLinReg(X = as_FBM(Z2), y.train = Y2) # calculate summary statistics of Y~Z
```

Then with the summary statistics, we can estimate the causal effect with `MR.SPI` and the corresponding conditional confidence interval:

```
library(igraph)

GammaHat = as.numeric(GWAS2$estim); gammaHat = as.numeric(GWAS1$estim) 
se_Gamma <- as.numeric(GWAS2$std.err); se_gamma <- as.numeric(GWAS1$std.err);

mr.spi.con <- MR.SPI(gammaHat, GammaHat, se_gamma, se_Gamma, n1, n2, max_clique = FALSE, unif = FALSE)
```

If we wish to construct the uniformly valid confidence interrval, then we can set `unif=TRUE`:

```
library(intervals)

mr.spi.sample <-MR.SPI(gammaHat, GammaHat, se_gamma, se_Gamma, n1, n2, max_clique = FALSE, unif = TRUE, Sampling = TRUE)
```
