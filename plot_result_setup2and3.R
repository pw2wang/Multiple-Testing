setwd("/Users/tedwang/Desktop/research_paper/simulation_MCF")
load('mcf_simulation1.RData')
load('mcf_simulation2.RData')

############ Simulation setup 1: EFT
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
t$type[t$type == 'storey'] = 'Storey'
t$type[t$type == 'mcf'] = 'MCF'

######### plot together 
case1_20 <-ggplot(t, aes(x = alpha, y = value, colour = type, linetype=type,shape=type)) +
  geom_line(size=0.6) + geom_point()+
#  labs(title = paste0("mu1 = mu2 = ",20)) +
  labs(y = "")+
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
t$type[t$type == 'storey'] = 'Storey'
t$type[t$type == 'mcf'] = 'MCF'

######### plot together 
case1_25 <- ggplot(t, aes(x = alpha, y = value, colour = type, linetype=type,shape=type)) +
  geom_line(size=0.6) + geom_point()+
#  labs(title = paste0("mu = ",25)) +
  labs(y = "")+
  facet_grid(result~.,scales = "free")
case2_25 + ggsave("setup1_25.png")
save.image("mcf_simulation1.RData")



##########################################################################
################ Simulation setup2: BT ###################################
##########################################################################



# FDR_res_10 <- lapply(1:N, function(x) res_10[[x]]$FDR)
# power_res_10 <- lapply(1:N, function(x) res_10[[x]]$power)
# FDR_res_10 <- as.data.frame(Reduce("+", FDR_res_10) / length(FDR_res_10))
# power_res_10 <- as.data.frame(Reduce("+", power_res_10) / length(power_res_10))
# 
# FDR_res_15 <- lapply(1:N, function(x) res_15[[x]]$FDR)
# power_res_15 <- lapply(1:N, function(x) res_15[[x]]$power)
# FDR_res_15 <- as.data.frame(Reduce("+", FDR_res_15) / length(FDR_res_15))
# power_res_15 <- as.data.frame(Reduce("+", power_res_15) / length(power_res_15))
# 
# t <- gather(FDR_res_25, type, value, ends_with("fdr"))
# 
# ggplot(t, aes(x = alpha, y = value, colour = type)) + 
#   geom_line() + labs(title="Power results")
# 
# 
# t <- gather(power_res_25, type, value, ends_with("power"))
# 
# ggplot(t, aes(x = alpha, y = value, colour = type)) + 
#   geom_line()+ ggtitle("FDR results")
# save.image("mcf_simulation1.RData")

######### mu = 10
FDR_res_10 <- lapply(1:N, function(x) res_10[[x]]$FDR)
power_res_10 <- lapply(1:N, function(x) res_10[[x]]$power)
FDR_res_10 <- as.data.frame(Reduce("+", FDR_res_10) / length(FDR_res_10))
power_res_10 <- as.data.frame(Reduce("+", power_res_10) / length(power_res_10))


FDR_res_10$result <- "FDR"
power_res_10$result <- "power"
names(FDR_res_10) <- names(power_res_10)
total_10 <- rbind(FDR_res_10,power_res_10)
t <- gather(total_10,type,value,ends_with("power"))
t$type <- str_remove(t$type,"_.*")
t$type[t$type == 'storey'] = 'Storey'
t$type[t$type == 'mcf'] = 'MCF'

######### plot together 
case2_10 <- ggplot(t, aes(x = alpha, y = value, colour = type, linetype=type,shape=type)) +
  geom_line(size=0.6) + geom_point()+
 # labs(title = paste0("mu = 10")) +
  labs(y = "")+
  facet_grid(result~.,scales = "free")
case2_10 + ggsave("setup2_10.png")



######### mu = 15
FDR_res_15 <- lapply(1:N, function(x) res_15[[x]]$FDR)
power_res_15 <- lapply(1:N, function(x) res_15[[x]]$power)
FDR_res_15 <- as.data.frame(Reduce("+", FDR_res_15) / length(FDR_res_15))
power_res_15 <- as.data.frame(Reduce("+", power_res_15) / length(power_res_15))


FDR_res_15$result <- "FDR"
power_res_15$result <- "power"
names(FDR_res_15) <- names(power_res_15)
total_15 <- rbind(FDR_res_15,power_res_15)
t <- gather(total_15,type,value,ends_with("power"))
t$type <- str_remove(t$type,"_.*")
t$type[t$type == 'storey'] = 'Storey'
t$type[t$type == 'mcf'] = 'MCF'

######### plot together 
case2_15 <- ggplot(t, aes(x = alpha, y = value, colour = type, linetype=type,shape=type)) +
  geom_line(size=0.6) + geom_point()+
  #labs(title = paste0("mu = 15")) +
  labs(y = "")+
  facet_grid(result~.,scales = "free")
case2_15
case2_15 + ggsave("setup2_15.png")


