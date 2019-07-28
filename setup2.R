setwd("\\\\files.math.uwaterloo.ca/y7chi/ResearchDesktop/ted/MCF")
files.sources = paste0("./helper/", list.files("./helper"))
invisible(sapply(files.sources, source))
library(qvalue)
library(fdrDiscreteNull)
#library(DiscreteFDR)
library(parallel)
library(pbapply)
library(ggplot2)

### get pvalues ####################
MT_compare <- function(i, pi_0,m, alpha,mu){
  ### data set-up ###############
  #m = 1000
  #pi_0 = 0.9
  #alpha = c(0.05,0.06)
  #alpha = 0.05
  lambda1 <- mu # or 20
  lambda2 <- mu # or 20
  ### fixed across repeated simulations ########################
  n.test = m
  N1 <- rpois(n.test, lambda1) # total of each binomial
  N2 <- rpois(n.test, lambda2)
  p1.Null <- runif(n.test*pi_0, min=0.1,max=0.9)
  p2.Null <- p1.Null
  p1.Nonnull <- 0.5 * runif(as.numeric(as.character(n.test*(1-pi_0))),min=0.1,max=0.9)
  p2.Nonnull <- p1.Nonnull + runif(as.numeric(as.character(n.test*(1-pi_0))),min=0.2,max=0.5)
  p1 <- c(p1.Null, p1.Nonnull)
  p2 <- c(p2.Null, p2.Nonnull)
  true_null_idx <- c(1:(0.9*n.test))
  non_null_idx <- c((0.9*n.test+1):n.test)
  
  ############
  N1.temp = N1
  N2.temp = N2
  p1.temp = p1
  p2.temp = p2
  #####################
  
  ### repeated simulations ##################
  #alpha = seq(0.01,0.1,by=0.01)
  #n.alpha = length(alpha.all)
  
  p.org <- rep(0,n.test) 
  p.next <- rep(0,n.test)
  #p.min <- rep(0,n.test)
  Z1 <- rep(0,n.test)
  Z2 <- rep(0,n.test)
  Z3 <- rep(0,n.test)
  Z4 <- rep(0,n.test)
  a <- sapply(1:length(N1), function(x) rbinom(1,N1[x],p1[x]))
  b <- sapply(1:length(N2), function(x) rbinom(1,N2[x],p2[x]))
  Z1 = a
  Z2 = N1 - Z1
  Z3 = b
  Z4 = N2 - Z3
  dat = as.matrix(cbind(Z1,Z3,Z1+Z2,Z3+Z4))
  
  cellcountsmarginals <- getcellcountsandmarginals_DE(dat)
  simallcellcounts <- cellcountsmarginals[[1]]
  pvsupp <-
    lapply(
      simallcellcounts,
      FUN = function(x)
        pvalFETSupport(x)
    )
  rawpvalues <- unlist(lapply(
    pvsupp,
    FUN = function(x)
      x$rawpvalues
  ))
  rndpvalues <- unlist(lapply(
    pvsupp,
    FUN = function(x)
      x$rndpvalues
  ))
  plessvalues <-
    unlist(lapply(
      pvsupp,
      FUN = function(x)
        x$plessvalues
    ))
  pvsupp <- lapply(
    pvsupp,
    FUN = function(x)
      x$support
  )  
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
res_25 <- pblapply(1:N,MT_compare,pi_0 = 0.9,m=5000,alpha=seq(0.01,0.25,0.01),mu=25,cl = cl)
res_20 <- pblapply(1:N,MT_compare,pi_0 = 0.9,m=5000,alpha=seq(0.01,0.25,0.01),mu=20,cl = cl)

stopCluster(cl)
gc()


######### mu = 20
FDR_res_20 <- lapply(1:N, function(x) res_20[[x]]$FDR)
power_res_20 <- lapply(1:N, function(x) res_20[[x]]$power)
FDR_res_20 <- as.data.frame(Reduce("+", FDR_res_20) / length(FDR_res_20))
power_res_20 <- as.data.frame(Reduce("+", power_res_20) / length(power_res_20))
FDR_res_20$result <- "FDR"
power_res_20$result <- "power"
names(FDR_res_20) <- names(power_res_20)
total_20 <- rbind(FDR_res_20,power_res_20)
t <- gather(total_20,type,value,ends_with("power"))
t$type <- str_remove(t$type,"_.*")
######### plot together 
case1_20 <- ggplot(t, aes(x = alpha, y = value, colour = type, linetype=type)) +
  geom_line(size=1.05) + 
  labs(title = paste0("mu1 = mu2 = ",20)) +
  facet_grid(result~.,scales = "free")

#theme(plot.title = element_text(color="black", size=14,
#                                face="bold.italic", hjust = 0.5),
#      legend.text=element_text(size=10,face="bold.italic"))
case1_20 + ggsave("setup1_20.png")

############# mu = 25
FDR_res_25 <- lapply(1:N, function(x) res_25[[x]]$FDR)
power_res_25 <- lapply(1:N, function(x) res_25[[x]]$power)
FDR_res_25 <- as.data.frame(Reduce("+", FDR_res_25) / length(FDR_res_25))
power_res_25 <- as.data.frame(Reduce("+", power_res_25) / length(power_res_25))
FDR_res_25$result <- "FDR"
power_res_25$result <- "power"
names(FDR_res_25) <- names(power_res_25)
total_25 <- rbind(FDR_res_15,power_res_25)
t <- gather(total_25,type,value,ends_with("power"))
t$type <- str_remove(t$type,"_.*")
######### plot together 
case1_25 <- ggplot(t, aes(x = alpha, y = value, colour = type, linetype=type)) +
  geom_line(size=1.05) + 
  labs(title = paste0("mu = ",25)) +
  facet_grid(result~.,scales = "free")
case2_25 + ggsave("setup1_25.png")
save.image("mcf_simulation1.RData")



# FDR_res_25 <- lapply(1:N, function(x) res_25[[x]]$FDR)
# power_res_25 <- lapply(1:N, function(x) res_25[[x]]$power)
# FDR_res_25 <- as.data.frame(Reduce("+", FDR_res_25) / length(FDR_res_25))
# power_res_25 <- as.data.frame(Reduce("+", power_res_25) / length(power_res_25))
# 
# FDR_res_20 <- lapply(1:N, function(x) res_20[[x]]$FDR)
# power_res_20 <- lapply(1:N, function(x) res_20[[x]]$power)
# FDR_res_20 <- as.data.frame(Reduce("+", FDR_res_20) / length(FDR_res_20))
# power_res_20 <- as.data.frame(Reduce("+", power_res_20) / length(power_res_20))
# 
# library(tidyverse)
# t <- gather(FDR_res_25, type, value, ends_with("fdr"))
# ggplot(t, aes(x = alpha, y = value, colour = type)) + 
#   geom_line() + labs(title="Power results")
# 
# 
# t <- gather(power_res_25, type, value, ends_with("power"))
# 
# ggplot(t, aes(x = alpha, y = value, colour = type)) + 
#   geom_line()+ ggtitle("FDR results")
# save.image("mcf_simulation1.RData")
# 
# 
