#require(devtools)
#install_version("DiscreteFDR", version = "1.0", repos = "http://cran.us.r-project.org")
#install.packages("DiscreteFDR")
library(sgof)
library(discreteMTP)
library(parallel)
library(ggplot2)
library(DiscreteFDR)
library(pbapply)
library(qvalue)
library(knitr)
library(fdrDiscreteNull)
library("splines")
library(adaptMT)
library(mgcv)

get_prev <- function(pval_supp, p_org) {
  # Assuming the pval_supp is already sorted ascendingly
  idx <- which(pval_supp == p_org)[1]
  if (idx > 1) {
    return(pval_supp[idx - 1])
  } else {
    return(0)
  }
}

mcf_detect <- function(lambda_star, p_prev, p_org, randp_ecdf) {
  cdf <- sum(randp_ecdf < lambda_star)/length(randp_ecdf)
  r <- unlist(lapply(1:length(p_org), FUN = function(x) (lambda_star > p_org[x])*1 + (lambda_star < p_prev[x])*0 + 
                       (lambda_star >= p_prev[x] && lambda_star <= p_org[x])*(lambda_star - p_prev[x])/(p_org[x] - p_prev[x])))
  detection <- which(r >= quantile(r, probs = 1 - cdf))
  return(detection)
}

qvalue1 <- function (p, fdr.level = NULL, pfdr = FALSE, lfdr.out = TRUE, 
                     pi0 = NULL, ...) 
{
  p_in <- qvals_out <- lfdr_out <- p
  rm_na <- !is.na(p)
  p <- p[rm_na]
  if (min(p) < 0 || max(p) > 1) {
    stop("p-values not in valid range [0, 1].")
  }
  else if (!is.null(fdr.level) && (fdr.level <= 0 || fdr.level > 
                                   1)) {
    stop("'fdr.level' must be in (0, 1].")
  }
  if (is.null(pi0)) {
    pi0s <- pi0est(p, ...)
  }
  else {
    if (pi0 > 0 && pi0 <= 1) {
      pi0s = list()
      pi0s$pi0 = pi0
    }
    else {
      pi0s = list()
      pi0s$pi0 = pi0
    }
  }
  m <- length(p)
  i <- m:1L
  o <- order(p, decreasing = TRUE)
  ro <- order(o)
  if (pfdr) {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m/(i * (1 - 
                                                        (1 - p[o])^m))))[ro]
  }
  else {
    qvals <- pi0s$pi0 * pmin(1, cummin(p[o] * m/i))[ro]
  }
  qvals_out[rm_na] <- qvals
  if (lfdr.out) {
    lfdr <- lfdr(p = p, pi0 = pi0s$pi0, ...)
    lfdr_out[rm_na] <- lfdr
  }
  else {
    lfdr_out <- NULL
  }
  if (!is.null(fdr.level)) {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, 
                   qvalues = qvals_out, pvalues = p_in, lfdr = lfdr_out, 
                   fdr.level = fdr.level, significant = (qvals <= fdr.level), 
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda, 
                   pi0.smooth = pi0s$pi0.smooth)
  }
  else {
    retval <- list(call = match.call(), pi0 = pi0s$pi0, 
                   qvalues = qvals_out, pvalues = p_in, lfdr = lfdr_out, 
                   pi0.lambda = pi0s$pi0.lambda, lambda = pi0s$lambda, 
                   pi0.smooth = pi0s$pi0.smooth)
  }
  class(retval) <- "qvalue"
  return(retval)
}



generate_data <- function(m,portion,q,m1_portion){
  #m = 800
  #portion = 0.8
  #q = 0.4
  #m1_portion <- 0.2
  m3 <- m*portion
  m1 <- m1_portion*(m-m3)
  m2 <- m - m1-m3
  p1 <- 0.01
  p2 <- 0.1
  p3 <- 0.1
  group1_m1 <- sapply(1:m1,function(x) rbinom(25,1,p1))
  group1_m2 <- sapply(1:m2,function(x) rbinom(25,1,p2))
  group1_m3 <- sapply(1:m3,function(x) rbinom(25,1,p3))
  group1 <- cbind(group1_m1,group1_m2,group1_m3)
  #group1 <- c(rbinom(m1,25,p1),rbinom(m2,25,p2),rbinom(m3,25,p3))
  y1 = 25
  group2_m1 <- sapply(1:m1,function(x) rbinom(25,1,p1))
  group2_m2 <- sapply(1:m2,function(x) rbinom(25,1,p2))
  group2_m3 <- sapply(1:m3,function(x) rbinom(25,1,q))
  group2 <- cbind(group2_m1,group2_m2,group2_m3)
  #group2 <- c(rbinom(m1,25,p1),rbinom(m2,25,p2),rbinom(m3,25,q))
  y2 = 25
  #ds <- cbind(group1,y1,group2,y2)
  ds1 <- cbind(colSums(group1),colSums(group2),y1,y2)
  #res <- fisher.pvalues.support(ds,input = 'marginal',alternative = 'two.sided')
  res <- GeneralizedFDREstimators(data = ds1, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = 0.05)
  res$ture_NULL <- c(1:(m1+m2))
  res$ture_Alter <- c((m1+m2+1):m)
  res$m <- m
  res$m3 <- m3
  res$m1 <- m1
  total <- rbind(group1,group2)
  res$x <- sapply(1:m,function(i) var(total[,i]))
  #res$x <- group1 + group2
  #res_temp <- GeneralizedFDREstimators(data = ds1, Test = "Fisher's Exact Test", FET_via = "IndividualMarginals", FDRlevel = 0.05)
  #res$Threshold <- res_temp$SARP$Threshold
  
  return(res)
}

calculate_power <- function(rejected_index,ture_Alter){
  return(sum(rejected_index %in% ture_Alter)/length(ture_Alter))
}
calculate_FDP <- function(rejected_index,true_NULL){
  return(sum(rejected_index %in% true_NULL)/max(length(rejected_index),1))
}

overall <- function(iter, m,portion,q,m1_portion){
  options(warn = -1)
  m = 800
  portion = 0.8
  q = 0.4
  m1_portion <- 0.2
  ds <- generate_data(m,portion,q,m1_portion)
  p_org <- ds$pvalues
  p_prev <- unlist(lapply(1:m, FUN = function(x, pval, pval_supp) get_prev(pval_supp[[x]], pval[x]), pval = p_org, pval_supp = ds$pvalSupp))
  p_org[p_org > 1] <- 1
  pCDFlist <- ds$pvalSupp
  
  ########### MCF
  B <- 200
  replicated_randp <- replicate(B, runif(m, p_prev, p_org))
  randp_ecdf <- as.vector(replicated_randp)
  lambda_storey <- 0.5
  pi0 <- mean(apply(replicated_randp, 2, FUN = function(x) (1 + sum(x > lambda_storey))/((1 - lambda_storey)*m)))
  lambda_star <-ds$SARP$Threshold
  mcf_detection <- mcf_detect(lambda_star, p_prev, p_org, randp_ecdf)
  
  ############## Adapt 
  x <-  data.frame(x = ds$x)
  dist <- beta_family()
  #formulas <- paste0("ns(x, df = ", 6:10, ")")
  formula <- paste0("s(x)")
  #models <- lapply(formulas, function(formula){
  #  piargs <- muargs <- list(formula = formula)
  #  gen_adapt_model(name = "glm", piargs = piargs, muargs = muargs)
  #})
  # Run adapt
  #ADAP <- adapt(x = x, pvals = p_org, models = models,
  #             dist = dist, nfits = 10)
  ADAP <- adapt_gam(x = x, pvals = p_org, pi_formulas = formula,
                    mu_formulas = formula, dist = dist, nfits = 5, alphas = 0.05)
  ADAP
  
  #print(str(ADAP))
  if(length(ADAP$rejs) == 0){
    ADAP = 0
  }else{
    ADAP = ADAP$rejs[[length(ADAP$rejs)]]
  }
  
  ######################
  DBH.su <- DBH(p_org, pCDFlist, ret.crit.consts = FALSE)
  ADBH.su <- ADBH(p_org, pCDFlist, ret.crit.consts = FALSE)
  pi_0 <- pi0est(p_org,lambda=0.5)
  #Storey <- qvalue(raw.pvalues,fdr.level = 0.05,lambda = 0.5)
  Storey <- qvalue1(p_org,fdr.level = 0.05,lambda = 0.5,pi0 = pi_0$pi0.lambda)
  Heyse <- p.discrete.adjust(p_org, pCDFlist, method = 'DBH', cutoff = 1, n = length(p_org))
  BH <- p.discrete.adjust(p_org, NULL, method = 'BH')
  
  
  ###### power
  power_HSU <- calculate_power(DBH.su$Indices,ds$ture_Alter)
  power_AHSU <- calculate_power(ADBH.su$Indices,ds$ture_Alter)
  power_Storey <- calculate_power(which(Storey$significant == T),ds$ture_Alter)
  power_Heyse <- calculate_power(which(Heyse <= 0.05),ds$ture_Alter)
  power_BH <- calculate_power(which(BH <= 0.05),ds$ture_Alter)
  power_MCF <- calculate_power(mcf_detection,ds$ture_Alter)
  power_APA <- calculate_power(ADAP,ds$ture_Alter)
  res_power <- c(ds$m,ds$m3,ds$m1,q,power_BH,power_Storey,power_Heyse,power_HSU,power_AHSU,power_MCF,power_APA)
  names(res_power) <- c("m","m3","m1","q","BH","Storey","Heyse","HSU",'AHSU',"MCF","ADAPT")
  
  ###### FDR
  FDR_HSU <- calculate_FDP(DBH.su$Indices,ds$ture_NULL)
  FDR_AHSU <- calculate_FDP(ADBH.su$Indices,ds$ture_NULL)
  FDR_Storey <- calculate_FDP(which(Storey$significant == T),ds$ture_NULL)
  FDR_Heyse <- calculate_FDP(which(Heyse <= 0.05),ds$ture_NULL)
  FDR_BH <- calculate_FDP(which(BH <= 0.05),ds$ture_NULL)
  FDR_MCF <- calculate_FDP(mcf_detection,ds$ture_NULL)
  FDR_ADA <- calculate_FDP(ADAP,ds$ture_NULL)
  
  res_FDR <- c(ds$m,ds$m3,ds$m1,q,FDR_BH,FDR_Storey,FDR_Heyse,FDR_HSU,FDR_AHSU,FDR_MCF,FDR_ADA)
  names(res_FDR) <- c("m","m3","m1",'q',"BH","Storey","Heyse","HSU",'AHSU','MCF','ADAPT')
  
  options(warn = 0)
  res <- list(res_power,res_FDR)
  names(res) <- c('power',"FDR")
  #res <- cbind(res_power,res_FDR)
  return(res)
}

#overall(1,800,0.8,0.4,0.2) ->t
#pblapply(1:2,overall,m=800,portion = 0.8,q=0.4,m1_portion =0.2,cl=NULL) -> t


calculate_res <- function(iter,m,portion,q,m1_portion,cl=NULL){
  res <- pblapply(1:iter,overall,m=m,portion = portion,q=q,m1_portion = m1_portion,cl=cl)
  res_power <- sapply(1:iter,function(x) res[[x]]$power)
  res_FDR <- sapply(1:iter,function(x) res[[x]]$FDR)
  res_power <- rowMeans(res_power)
  res_FDR <- rowMeans(res_FDR)
  final_res <- list(res_power,res_FDR)
  names(final_res) <- c('power','FDR')
  return(final_res)
}
#calculate_res(2,800,0.8,0.4,0.2) -> t

trans_ds <- function(ds){
  #ds <- power_800
  ds <- ds[,5:11]
  row.names(ds) <- c("small","intermeidate","large prop.")
  new_ds <- data.frame(power = as.numeric(as.matrix(ds)),method =as.vector(sapply(1:length(colnames(ds)), function(x){rep(colnames(ds)[x],3)})),
                       scenerio = rep(row.names(ds),7))
  names(new_ds)[1] <- "Average Power"
  new_ds$method <- factor(new_ds$method,levels= c("BH","Storey","HSU","Heyse","AHSU","MCF","ADAPT"))
  new_ds$scenerio <- factor(new_ds$scenerio,levels = c("small","intermeidate","large prop."))
  return(new_ds)
}


##### q = 0.4, m =800

N <- 2
portion_of_m3 <- c(0.1,0.3,0.8)
#q_option <- c(0.15,0.25,0.4)
q_option <- 0.4
portion_of_m1 <- c(0.2,0.5,0.8)


#m_size <- c(800,2000)
m_size <- 800
#
#calculate_average_power(50,800,0.8,0.4,0.2,cl=cl)

#
cl <- makeCluster(detectCores())
clusterExport(cl,ls())
clusterEvalQ(cl,library(qvalue))
clusterEvalQ(cl,library(discreteMTP))
clusterEvalQ(cl,library(DiscreteFDR))
clusterEvalQ(cl,library(fdrDiscreteNull))
clusterEvalQ(cl,library(adaptMT))
#clusterEvalQ(cl,library(splines))
clusterEvalQ(cl,library(mgcv))

power_res <- NULL
FDR_res <- NULL
for (m in m_size) {
  for(portion in portion_of_m3){
    for (pi_m1 in portion_of_m1) {
      for (q in q_option) {
        print(paste0("m is: ",m,"portion is: ",portion,"pi_m1 is: ",pi_m1,"q is :", q))
        temp <- calculate_res(N,m,portion,q,pi_m1,cl=cl)
        power_res <- rbind(power_res,temp$power)
        FDR_res <- rbind(FDR_res,temp$FDR)
      }      
    }
  }
}

power_res <- as.data.frame(power_res)
FDR_res <- as.data.frame(FDR_res)
stopCluster(cl)
gc()

#kable(power_res,format = "latex")
#kable(FDR_res,format = "latex")

######### transform datafrom to plot 

kable(round(power_res,4),format = "latex")
kable(round(FDR_res,4),format = "latex")

power_800 <- power_res[power_res$m == 800 & power_res$q == 0.4 & 
                         power_res$m1/(power_res$m -power_res$m3) ==0.2,]

power_2000 <- power_res[power_res$m == 2000 & power_res$q == 0.4 & 
                          power_res$m1/(power_res$m -power_res$m3) ==0.2,]


power_800_new <- trans_ds(power_800)
power_2000_new <- trans_ds(power_2000)

power_800 <- power_res[power_res$m == 800 & power_res$q == 0.25 & 
                         power_res$m1/(power_res$m -power_res$m3) ==0.2,]

power_2000 <- power_res[power_res$m == 2000 & power_res$q == 0.25 & 
                          power_res$m1/(power_res$m -power_res$m3) ==0.2,]


power_800_new <- trans_ds(power_800)
power_2000_new <- trans_ds(power_2000)

ggplot(data=power_800_new, aes(x=scenerio, y=`Average Power`, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  ggtitle("m = 800") + ylab("Average power") +
  theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5))


ggplot(data=power_2000_new, aes(x=scenerio, y=`Average Power`, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  ggtitle("m = 2000") + ylab("Average power") +
  theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5))



fdr_2000 <- FDR_res[FDR_res$m == 2000 & FDR_res$q == 0.25 & 
                          FDR_res$m1/(FDR_res$m -FDR_res$m3) ==0.2,]


