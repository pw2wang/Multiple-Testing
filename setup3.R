#setwd("\\\\files.math.uwaterloo.ca/y7chi/ResearchDesktop/ted/MCF")
setwd("/Users/tedwang/Desktop/research_paper/simulation_MCF/")
files.sources = paste0("./helper/", list.files("./helper"))
invisible(sapply(files.sources, source))
library(qvalue)
library(fdrDiscreteNull)
library(DiscreteFDR)
library(parallel)
library(pbapply)
library(ggplot2)
library(tidyverse)

### get pvalues ####################
MT_compare <- function(i, pi_0,m, alpha,mu){
  #mu = 10
  #m = 5000
  #pi_0 = 0.9
  #alpha = 0.05
  
  lambda1 = mu
  n.test = m
  N <- rpois(n.test, lambda1) + 1 # total of the binomial
  num_null = pi_0 * m
  num_non_null = round((1-pi_0) * m)
  true_null_idx <- c(1:(0.9*n.test))
  non_null_idx <- c((0.9*n.test+1):n.test)
  P.all <- c(rep(0.5,num_null), 0.5 + sample(c(-1,1),num_non_null,replace = TRUE)*runif(num_non_null,min=0.2,max=0.5))
  X <- sapply(1:n.test, function(x) rbinom(1,N[x],P.all[x]))  
  dat <- as.matrix(cbind(X,N-X))
  ResTmp = GeneralizedEstimatorsGrouped(data=dat,
                                        test_in = "Binomial Test",FET_via_in = "IndividualMarginals",
                                        grpby = 'quantileOfRowTotal', ngrp_in = 3,
                                        FDRlevel_in = 0.05) # , lambda_in = 0.5, epsilon_in = 1)  
  rawpvalues <- ResTmp$pvalues
  pvsupp <- ResTmp$pvalSupp
  rawpvalues[rawpvalues > 1] <- 1
  pvsupp <- lapply(pvsupp, function(x) {x[x>1] = 1;return(x)})
  
  plessvalues <- unlist(lapply(1:length(pvsupp), function(x) get_prev(pvsupp[[x]], rawpvalues[x])))
  eqProb = rawpvalues - plessvalues
  #rndpvalues <- plessvalues + runif(1)*eqProb
  #rndpvalues <- min(1,rndpvalues)
  rndpvalues <- sapply(1:length(rawpvalues), function(x) plessvalues[x] + runif(1)*eqProb[x])
  
  
  
  ##########################################################
  ###################        BH        #####################
  ##########################################################
  BH_detection <- lapply(
    alpha,
    FUN = function(x, pval)
      BHFDRApp(pval, x)[[1]][, 2],
    pval = rawpvalues
  )
  
  ##########################################################
  ###################        BHH        ####################
  ##########################################################
  # Heyse adjusted pvalue
  heysepvalues <- HeyseAdjFunc(rawpvalues, pvsupp)
  Heyse_detection <- lapply(
    alpha,
    FUN = function(x, pval)
      FDRAdjustedPval(pval, x)[[1]][, 2],
    pval = heysepvalues
  )
  
  ##########################################################
  ###################        AHSU        ###################
  ##########################################################
  pvec <- match.pvals(pvsupp, rawpvalues)
  sorted.pvals <- sort(pvec)
  stepf <- build.stepfuns(pvsupp)
  pv.list.all <- sort(unique(as.numeric(unlist(pvsupp))))
  pv.list.c.m <- lapply(
    alpha,
    FUN = function(x, ls)
      short.eff(ls, x / (1 + x)),
    ls = pv.list.all
  )
  c.m.s <- unlist(
    lapply(
      1:length(alpha),
      FUN = function(x, stepf, ls, alpha)
        find.cm(alpha[x], stepf, ls[[x]]),
      stepf = stepf,
      ls = pv.list.c.m,
      alpha = alpha
    )
  )
  obs.pvals <- lapply(
    c.m.s,
    FUN = function(x, pvals)
      pvals[pvals <= x],
    pvals = sorted.pvals
  )
  AHSU_detection <- lapply(
    1:length(alpha),
    FUN = function(x, pvals, c.m.s, stepf, alpha)
      max(which(
        kernel.ADBH.fast(stepf, pvals[[x]], c.m.s[x]) <= 1:length(pvals[[x]]) * alpha[x]
      )),
    pvals = obs.pvals,
    c.m.s = c.m.s,
    stepf = stepf,
    alpha = alpha
  )
  AHSU_detection <- lapply(
    1:length(alpha),
    FUN = function(x, pval, obs.pvals, indices)
      which(pval <= (obs.pvals[[x]])[indices[[x]]]),
    pval = rawpvalues,
    obs.pvals = obs.pvals,
    indices = AHSU_detection
  )
  ##########################################################
  ###################      Storey      #####################
  ##########################################################
  lambda_storey <- 0.5
  pi0_storey <- storeyPi0Est(lambda_storey, rawpvalues)
  storey_detection <- lapply(
    alpha,
    FUN = function(x, pval, pi_0)
      StoreyFDREst(pval, x, pi_0)[[1]],
    pval = rawpvalues,
    pi_0 = pi0_storey
  )
  
  
  
  ##########################################################
  ###############     mcf detection         ################
  ##########################################################
  lambda_storey <- 0.5
  pi0_habiger <- storeyPi0Est(lambda_storey, rndpvalues)
  habiger_res <- lapply(
    alpha,
    FUN = function(x, pval, pi_0)
      StoreyFDREst(pval, x, pi_0),
    pval = rndpvalues,
    pi_0 = pi0_habiger
  )
  habiger_thres <- lapply(
    habiger_res,
    FUN = function(x)
      x[[2]][1]
  )
  habiger_detection <- lapply(
    habiger_res,
    FUN = function(x)
      x[[1]]
  )
  
  p_org <- rawpvalues
  p_prev <- plessvalues
  
  B <- 200
  B_ecdf <- B * m
  replicated_randp <- replicate(B, runif(m, p_prev, p_org))
  randp_ecdf <- as.vector(replicated_randp)
  pi0 <- mean(apply(
    replicated_randp,
    2,
    FUN = function(x)
      (1 + sum(x > lambda_storey)) / ((1 - lambda_storey) *
                                        m)
  ))
  
  lambda_star <- unlist(habiger_thres)
  mcf_detection <-
    lapply(
      lambda_star,
      FUN = function(x, p_prev, p_org, randp_ecdf)
        mcf_detect(x, p_prev, p_org, randp_ecdf),
      p_prev = p_prev,
      p_org = p_org,
      randp_ecdf = randp_ecdf
    )
  
  
  #####################
  ####### Results 
  #####################
  
  BH_power <-
    unlist(lapply(
      BH_detection,
      FUN = function(x)
        get_power(x, non_null_idx)
    ))
  Heyse_power <-
    unlist(lapply(
      Heyse_detection,
      FUN = function(x)
        get_power(x, non_null_idx)
    ))
  AHSU_power <-
    unlist(lapply(
      AHSU_detection,
      FUN = function(x)
        get_power(x, non_null_idx)
    ))
  mcf_power <-
    unlist(lapply(
      mcf_detection,
      FUN = function(x)
        get_power(x, non_null_idx)
    ))
  storey_power <-
    unlist(lapply(
      storey_detection,
      FUN = function(x)
        get_power(x, non_null_idx)
    ))
  
  BH_fdr <-
    unlist(lapply(
      BH_detection,
      FUN = function(x)
        get_fdr(x, true_null_idx)
    ))
  Heyse_fdr <-
    unlist(lapply(
      Heyse_detection,
      FUN = function(x)
        get_fdr(x, true_null_idx)
    ))
  AHSU_fdr <-
    unlist(lapply(
      AHSU_detection,
      FUN = function(x)
        get_fdr(x, true_null_idx)
    ))
  
  mcf_fdr <-
    unlist(lapply(
      mcf_detection,
      FUN = function(x)
        get_fdr(x, true_null_idx)
    ))
  storey_fdr <-
    unlist(lapply(
      storey_detection,
      FUN = function(x)
        get_fdr(x, true_null_idx)
    ))
  #methods <- c("BH","Storey","HSU","Heyse","AHSU","MCF")
  res_FDR <-cbind(
    m,
    pi_0,
    alpha,
    BH_fdr,
    storey_fdr,
    Heyse_fdr,
    AHSU_fdr,
    mcf_fdr
  )
  
  res_power <- cbind(
    m,
    pi_0,
    alpha,
    BH_power,
    storey_power,
    Heyse_power,
    AHSU_power,
    mcf_power)
  
  res_total <- list(res_FDR, res_power)
  names(res_total) <- c("FDR", "power")
  return(res_total)
}
N <- 200
cl <- makeCluster(detectCores())  
clusterExport(cl,ls())
clusterEvalQ(cl,library(DiscreteFDR))
clusterEvalQ(cl,library(fdrDiscreteNull))
res_10 <- pblapply(1:N,MT_compare,pi_0 = 0.9,m=5000,alpha=seq(0.01,0.25,0.01),mu=10,cl = cl)
res_15 <- pblapply(1:N,MT_compare,pi_0 = 0.9,m=5000,alpha=seq(0.01,0.25,0.01),mu=15,cl = cl)

stopCluster(cl)
gc()

FDR_res_10 <- lapply(1:N, function(x) res_10[[x]]$FDR)
power_res_10 <- lapply(1:N, function(x) res_10[[x]]$power)
FDR_res_10 <- as.data.frame(Reduce("+", FDR_res_10) / length(FDR_res_10))
power_res_10 <- as.data.frame(Reduce("+", power_res_10) / length(power_res_10))

FDR_res_15 <- lapply(1:N, function(x) res_15[[x]]$FDR)
power_res_15 <- lapply(1:N, function(x) res_15[[x]]$power)
FDR_res_15 <- as.data.frame(Reduce("+", FDR_res_15) / length(FDR_res_15))
power_res_15 <- as.data.frame(Reduce("+", power_res_15) / length(power_res_15))

t <- gather(FDR_res_25, type, value, ends_with("fdr"))
ggplot(t, aes(x = alpha, y = value, colour = type)) + 
  geom_line() + labs(title="Power results")


t <- gather(power_res_25, type, value, ends_with("power"))

ggplot(t, aes(x = alpha, y = value, colour = type)) + 
  geom_line()+ ggtitle("FDR results")
save.image("mcf_simulation1.RData")


