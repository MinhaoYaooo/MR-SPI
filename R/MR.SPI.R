#' @title  Mendelian randomization that first selects valid IVs and then perform post-selection inference
#'
#' @description This package is used for conduting two-sample Mendelian randomization with indepnendent SNPs by selecting valid IVs that satisfy the core assumptions and inferring the causal effect.
#' @param gamma a numeric vector of GWAS summary statistics of the exposure
#' @param Gamma a numeric vector of GWAS summary statistics of the outcome
#' @param se_gamma a numeric vector of standard errors of gamma
#' @param se_Gamma a numeric vector of standard errors of Gamma
#' @param n1 the sample size of GWAS summary statistics  of the exposure
#' @param n2 the sample size of GWAS summary statistics  of the outcome
#' @param freq a numeric vector of allele frequencies to adjust for the regression coefficients and the standard errors, with default NULL
#' @param tuning a numeric scalar value tuning parameter for selecting the valid IVs, with default 1
#' @param max_clique an option to replace the majority and plurality voting procedures with finding maximal clique in the IV voting matrix, with default TRUE
#' @param unif an option to produce uniformly valid confidence interval, with default FALSE
#' @param alpha a numeric scalar value between 0 and 1 indicating the significance level for the confidence interval, with default 0.05
#' @param a grid size for constructing beta grids, with default 0.6
#' @param Sampling use the sampling method if TRUE; else use the proposed searching method, with default TRUE (only applicable if unif=TRUE)
#' @param rho a numeric scalar denoting thresholding level used in the sampling property (only applicable if unif=TRUE)
#' @param M sampling times, with default value 1000 (only applicable if unif=TRUE)
#' @param prop proportion of intervals kept when sampling, with default value 0.1 (only applicable if unif=TRUE)
#' @param filtering Filtering the resampled data or not, with default value TRUE (only applicable if unif=TRUE)
#' @return \code{MR.SPI} return a list containing the following:
#'     \item{\code{betaHat}}{a numeric scalar denoting the estimate of causal effect.}
#'     \item{\code{beta.sdHat}}{a numeric scalar denoting the estimated standard error of betaHat.}
#'     \item{\code{ci}}{a two dimensional numeric vector denoting the 1-alpha confidence intervals for betaHat with lower and upper endpoints.}
#'     \item{\code{SHat}}{a numeric vector denoting the set of relevant IVs}
#'     \item{\code{VHat}}{a numeric vector denoting the set of valid IVs}
#'     \item{\code{voting.mat}}{a matrix of voting. if (i,j)-th entry=1, then i-th SNP and j-th SNP vote for each other to be a valid IV.}
#' @export
MR.SPI <- function(gamma, Gamma, se_gamma, se_Gamma, n1, n2, freq=NULL, tuning=1, max_clique=TRUE, unif=FALSE,
                   alpha=0.05, a=0.6, Sampling=TRUE, rho=NULL, M=1000, prop=0.1, filtering=TRUE){

  if(!is.null(freq)){
    gamma <- gamma * sqrt(2*freq*(1-freq))
    Gamma <- Gamma * sqrt(2*freq*(1-freq))
    se_gamma <- se_gamma * sqrt(2*freq*(1-freq))
    se_Gamma <- se_Gamma * sqrt(2*freq*(1-freq))
  }

  pz = length(Gamma)

  V_gamma <- matrix(nrow = pz, ncol = pz)
  V_Gamma <- matrix(nrow = pz, ncol = pz)

  for (i in 1:pz) {
    for (j in i:pz) {
      if(i==j){
        V_gamma[i,i] <- se_gamma[i]^2
        V_Gamma[i,i] <- se_Gamma[i]^2
      }else{
        V_gamma[i,j] <- V_gamma[j,i] <- gamma[i]*gamma[j] / n1
        V_Gamma[i,j] <- V_Gamma[j,i] <- Gamma[i]*Gamma[j] / n2
      }
    }
  } # Construct V_gamma_hat and V_Gamma_hat based on the formula

  SHat <- (1:pz)[(abs(gamma) >= (sqrt(log(n1)*tuning*diag(V_gamma))))]

  SHat.bool = rep(FALSE,pz); SHat.bool[SHat] = TRUE

  nCand = length(SHat)
  VHats.bool = matrix(FALSE,nCand,nCand); colnames(VHats.bool) = rownames(VHats.bool) = SHat

  for(j in SHat) {
    beta.j <- Gamma[j] / gamma[j]
    pi.j = (Gamma * gamma[j]  - gamma * Gamma[j]) / gamma[j]
    sigmasq.j = diag(V_Gamma)+(gamma/gamma[j])^2*V_Gamma[j,j]-2*gamma/gamma[j]*V_Gamma[,j] +
      beta.j^2*(diag(V_gamma)+(gamma/gamma[j])^2*V_gamma[j,j]-2*gamma/gamma[j]*V_gamma[,j])
    PHat.bool.j = abs(pi.j) <= sqrt(sigmasq.j*tuning^2*log(min(n1,n2)))
    VHat.bool.j = PHat.bool.j * SHat.bool
    VHats.bool[as.character(SHat),as.character(j)] = VHat.bool.j[SHat]
  }

  VHats.boot.sym<-VHats.bool
  for(i in 1:dim(VHats.boot.sym)[1]){
    for(j in 1:dim(VHats.boot.sym)[2]){
      VHats.boot.sym[i,j]<-min(VHats.bool[i,j],VHats.bool[j,i])
    }
  }

  diag(VHats.boot.sym) <- 1
  # Voting method
  VM= apply(VHats.boot.sym,1,sum)
  VM.m = rownames(VHats.boot.sym)[VM > (0.5 * length(SHat))] # Majority winners
  VM.p = rownames(VHats.boot.sym)[max(VM) == VM]

  if (max_clique) {
    voting.graph <- igraph::as.undirected(igraph::graph_from_adjacency_matrix(VHats.boot.sym))
    max.clique <- igraph::largest.cliques(voting.graph)
    VHat <- unique(igraph::as_ids(Reduce(c,max.clique))) # take the union if multiple max cliques exist
    VHat <- sort(as.numeric(VHat))
  } else{
    V.set<-NULL
    for(index in VM.p){
      V.set<-union(V.set,names(which(VHats.boot.sym[index,]==1)))
    }
    VHat<-NULL
    for(index in V.set){
      VHat<-union(VHat,names(which(VHats.boot.sym[,index]==1)))
    }
    VHat=sort(as.numeric(VHat))
  }

  if (max_clique) {
    max.clique.mat <- matrix(0,nrow = length(max.clique),ncol = length(max.clique[[1]]))
    CI.temp <- matrix(0,nrow = length(max.clique), ncol = 2)
    beta.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    betavar.temp <- matrix(0,nrow = length(max.clique), ncol = 1)
    for (i in 1:length(max.clique)) {
      temp <- SHat[sort(as.numeric(max.clique[[i]]))]
      max.clique.mat[i,] <- temp
      betaHat = sum(gamma[temp] * Gamma[temp]) / sum(gamma[temp]^2)
      betaVarHat = var_TSHT(gamma[temp], Gamma[temp], V_gamma[temp,temp], V_Gamma[temp,temp])
      ci = c(betaHat - stats::qnorm(1-alpha/2) * sqrt(betaVarHat),betaHat + stats::qnorm(1-alpha/2) * sqrt(betaVarHat))
      CI.temp[i,] <- ci
      beta.temp[i,] <- betaHat
      betavar.temp[i,] <- betaVarHat
    }
    uni<- intervals::Intervals(CI.temp)
    ###### construct the confidence interval by taking a union
    CI.union<-as.matrix(intervals::interval_union(uni))

  } else{
    betaHat <- sum(gamma[VHat] * Gamma[VHat]) / sum(gamma[VHat]^2)

    betaVarHat <- var_TSHT(gamma[VHat], Gamma[VHat], V_gamma[VHat,VHat], V_Gamma[VHat,VHat])

    ci = c(betaHat-stats::qnorm(1-alpha/2)*sqrt(betaVarHat), betaHat+stats::qnorm(1-alpha/2)*sqrt(betaVarHat))
  }

  if(unif){

    var.beta = diag(V_Gamma)/gamma^2 + diag(V_gamma)*Gamma^2/gamma^4
    var.beta = var.beta[VHat]
    CI.init = matrix(NA, nrow=length(VHat), ncol=2)
    CI.init[,1] = (Gamma/gamma)[VHat] - sqrt(log(n1)*var.beta)
    CI.init[,2] = (Gamma/gamma)[VHat] + sqrt(log(n1)*var.beta)
    uni = intervals::Intervals(CI.init)
    CI.init.union = as.matrix(intervals::interval_union(uni))
    beta.grid = grid.CI(CI.init.union, grid.size=(min(n1,n2))^{-a})

    if(Sampling){
      ## Sampling Method
      CI.sampling = Searching.CI.sampling(gamma, Gamma, V_gamma, V_Gamma, InitiSet=VHat, n1=n1, n2=n2, beta.grid = beta.grid,
                                          alpha=alpha, rho=rho, M=M, prop=prop, filtering=filtering)
      CI=CI.sampling$CI
      rule=CI.sampling$rule
    }else{
      ## Searching Method
      CI.searching = Searching.CI(gamma, Gamma, V_gamma, V_Gamma, InitiSet = VHat, alpha=alpha,
                                  beta.grid = beta.grid)
      CI=CI.searching$CI
      rule=CI.searching$rule
    }

  }


  if(unif){
    returnList <- list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat),ci=CI,SHat=SHat,VHat=VHat,voting.mat=VHats.boot.sym)
  }else{
    if (max_clique) {
      returnList <- list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat),ci=CI.union,SHat=SHat,VHat=VHat,voting.mat=VHats.boot.sym)
    } else {
      returnList <- list(betaHat=betaHat,beta.sdHat = sqrt(betaVarHat),ci=ci,SHat=SHat,VHat=VHat,voting.mat=VHats.boot.sym)
    }
  }

  cat(paste0(pz,' SNPs used as candidate instruments. MR.SPI identifies ',length(SHat),' relevant instruments and ',length(VHat),' valid instruments.\n'))
  cat(paste0('Estimated Causal Effect by MR.SPI: ', round(betaHat,3),'\n'))
  if(unif){
    cat(paste0('Uniform Confidence Interval: (',round(returnList$ci[1],3),' , ',round(returnList$ci[2],3),')\n\n'))
  }else{
    cat(paste0('Conditional Confidence Interval: (',round(returnList$ci[1],3),' , ',round(returnList$ci[2],3),')\n\n'))
  }

  return(returnList)
}



########## Helpful Functions ##########

var_TSHT <- function(gamma, Gamma, V_gamma, V_Gamma){

  Var_Gg <- t(gamma) %*% V_Gamma %*% gamma + t(Gamma) %*% V_gamma %*% Gamma
  Var_gg <- 4 * t(gamma) %*% V_gamma %*% gamma
  Cov_Gg <- 2 * t(Gamma) %*% V_gamma %*% gamma

  E_Gg <- sum(gamma*Gamma)
  E_gg <- sum(gamma^2) + sum(diag(V_gamma))

  sig2 <- (E_Gg/E_gg)^2 * (Var_Gg/E_Gg^2 + Var_gg / E_gg^2 - 2*Cov_Gg/(E_Gg*E_gg) )

  return(sig2)
}


grid.CI <- function(CI.matrix, grid.size){
  d = dim(CI.matrix)[1]
  grid.seq = NULL
  for(l in 1:d) grid.seq = c(grid.seq, seq(CI.matrix[l, 1], CI.matrix[l, 2], by=grid.size))
  return(grid.seq)
}

Searching.CI <- function(gamma, Gamma, V_gamma, V_Gamma, alpha=0.05, InitiSet, beta.grid){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V_Gamma)[1]

  ## new rho method
  Tn = stats::qnorm(1-alpha/(2*pz))

  ## valid grid
  valid.grid = rep(NA, n.beta)
  for(j in 1:n.beta){
    b = beta.grid[j]
    temp = sqrt(diag(V_Gamma + b^2*V_gamma))
    valid.grid[j] = sum(abs(Gamma[InitiSet] - b*gamma[InitiSet]) < Tn*temp[InitiSet])
  }

  ## select beta
  if(length(beta.grid[which(valid.grid > threshold.size)])==0){
    rule=FALSE
    warning("Rule Fails. SS will give misleading CIs, SEs, and p-values.")
    sel.index = which(valid.grid==max(valid.grid))
  }else{
    rule=TRUE
    sel.index = which(valid.grid>threshold.size)
  }
  CI = t(as.matrix(c(min(beta.grid[sel.index]), max(beta.grid[sel.index]))))

  return(list(CI=CI, rule=rule))
}

Searching.CI.sampling <- function(gamma, Gamma, V_gamma, V_Gamma, InitiSet, n1, n2,
                                  beta.grid, alpha=0.05, rho=NULL, M=1000, prop=0.1, filtering=TRUE){
  threshold.size = length(InitiSet)/2
  n.beta = length(beta.grid)
  pz = dim(V_Gamma)[1]

  ## new rho method
  Tn = stats::qnorm(1-alpha/(2*pz))

  ## Covariance Matrix
  Cov1 <- cbind(V_Gamma, matrix(0, pz, pz))
  Cov2 <- cbind(matrix(0, pz, pz), V_gamma)
  Cov.total <- rbind(Cov1,Cov2)

  valid.grid.sample = matrix(NA, nrow=M, ncol=n.beta)
  Gen.mat = MASS::mvrnorm(M, rep(0, 2*pz), Cov.total)
  if(is.null(rho)) rho = (log(min(n1,n2))/M)^(1/(2*length(InitiSet)))/6 # initial rho if not specified

  if(filtering){
    temp = abs(t(t(Gen.mat[,(pz+1):(2*pz)]) / sqrt(diag(V_gamma))))
    temp = apply(temp, MARGIN=1, FUN=max)
    temp = (temp <= stats::qnorm(1-alpha/(2*pz)))
    Gen.mat = Gen.mat[temp,]
    M = sum(temp)
  }

  while(rho < 0.5){
    for(m in 1:M){
      Gamma.sample = Gamma - Gen.mat[m, 1:pz]
      gamma.sample = gamma - Gen.mat[m, (pz+1):(2*pz)]
      for(j in 1:n.beta){
        b = beta.grid[j]
        temp = sqrt(diag(V_Gamma + b^2*V_gamma))
        valid.grid.sample[m, j] = sum(abs(Gamma.sample[InitiSet] - b*gamma.sample[InitiSet])<
                                        rho*Tn*temp[InitiSet])
      }
    }
    CI = matrix(NA, nrow=M, ncol=2)
    for(m in 1:M){
      if(length(which(valid.grid.sample[m, ] > threshold.size))>0){
        CI[m, 1] = min(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
        CI[m, 2] = max(beta.grid[which(valid.grid.sample[m, ] > threshold.size)])
      }
    }
    CI = CI[!rowSums(is.na(CI)), , drop=FALSE] # CI = na.omit(CI)

    ## check CI dims, stop iterations if CI dim is big enough
    if(dim(as.matrix(CI))[1] >= prop*M) break
    rho = 1.25 * rho # increase rho with iterations
  }

  rule = TRUE
  if(dim(as.matrix(CI))[1] < prop*M){
    warning("Sampling Criterion not met, trasfer to Searching Method.")
    CI.searching = Searching.CI(gamma, Gamma, V_gamma, V_Gamma, InitiSet, beta.grid)
    rule = CI.searching$rule
    CI = CI.searching$CI
  }else{
    uni = intervals::Intervals(CI)
    CI = as.matrix(intervals::interval_union(uni))
    CI = t(as.matrix(c(min(CI[,1]), max(CI[,2]))))
  }

  return(list(CI=CI, rule=rule))
}

