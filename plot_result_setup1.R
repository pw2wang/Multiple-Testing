setwd("/Users/tedwang/Desktop/research_paper/code")
load("new_simulation_ns.RData")
library(knitr)
library(ggplot2)


#kable(round(power_res[,-5],4),format = "latex")
#kable(round(FDR_res[,-5],4),format = "latex")

######### transform datafrom to plot 
trans_ds <- function(ds){
  #ds <- power_800
  #ds <- ds[,6:12]
  ds <- ds[,6:11]
  row.names(ds) <- c("small","intermeidate","large prop.")
  new_ds <- data.frame(power = as.numeric(as.matrix(ds)),method =as.vector(sapply(1:length(colnames(ds)), function(x){rep(colnames(ds)[x],3)})),
                       scenerio = rep(row.names(ds),6))
  names(new_ds)[1] <- "Average Power"
  #new_ds$method <- factor(new_ds$method,levels= c("BH","Storey","HSU","Heyse","AHSU","MCF","ADAPT"))
  new_ds$method <- factor(new_ds$method,levels= c("BH","Storey","HSU","Heyse","AHSU","MCF"))
  new_ds$scenerio <- factor(new_ds$scenerio,levels = c("small","intermeidate","large prop."))
  return(new_ds)
}

power_res$ADAPT <- NULL
FDR_res$ADAPT <- NULL
power_res$alpha_level <- NULL
FDR_res$alpha_level <- NULL
kable(round(power_res,4),format = "latex")
kable(round(FDR_res,4),format = "latex")

power_800 <- power_res[power_res$m == 800 & power_res$q == 0.4 & 
                         power_res$m1/(power_res$m -power_res$m3) ==0.2,]
FDR_800 <- FDR_res[FDR_res$m == 800 & FDR_res$q == 0.4 & 
                     FDR_res$m1/(FDR_res$m -FDR_res$m3) ==0.2,]

power_2000 <- power_res[power_res$m == 2000 & power_res$q == 0.4 & 
                          power_res$m1/(power_res$m -power_res$m3) ==0.2,]
FDR_2000 <- FDR_res[FDR_res$m == 2000 & FDR_res$q == 0.4 & 
                     FDR_res$m1/(FDR_res$m -FDR_res$m3) ==0.2,]



power_800_new <- trans_ds(power_800)
power_2000_new <- trans_ds(power_2000)
FDR_800_new <- trans_ds(FDR_800)
FDR_2000_new <- trans_ds(FDR_2000)
FDR_800_new$m <- 800
FDR_2000_new$m <- 2000

power_800_new$m <- 800
power_2000_new$m <- 2000
power_tog <- rbind(power_800_new,power_2000_new)
FDR_tog <- rbind(FDR_800_new,FDR_2000_new)
names(FDR_tog)[1] <- 'FDR'

g_tog <- ggplot(data=power_tog, aes(x=scenerio, y=`Average Power`, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme(plot.title = element_text(color="black", size=10, face="bold.italic", hjust = 0.5),legend.text=element_text(size=10,face="bold.italic")) +
  facet_grid(.~m,scales = "free")
g_tog + ggsave("power_together.png")


g_tog <- ggplot(data=FDR_tog, aes(x=scenerio, y=FDR, fill=method)) +
  geom_bar(stat="identity", position=position_dodge())+
  scale_fill_brewer(palette="Paired")+
  theme(plot.title = element_text(color="black", size=10, face="bold.italic", hjust = 0.5),legend.text=element_text(size=10,face="bold.italic")) +
  facet_grid(.~m,scales = "free")
g_tog + ggsave("FDR_together.png")

# g_800 <- ggplot(data=power_800_new, aes(x=scenerio, y=`Average Power`, fill=method)) +
#   geom_bar(stat="identity", position=position_dodge())+
#   scale_fill_brewer(palette="Paired")+
#   #ggtitle(paste0("m : ",power_800$m[1],", pi_m1 : ",power_800$m1[1]/(power_800$m[1] - power_800$m3[1]),
#   #               ", q  :", power_800$q[1])) + ylab("Average power") +
#   ggtitle("m = 800") +
#   theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),legend.text=element_text(size=14,face="bold.italic"))
# 
# g_800 + ggsave("power_800.png")
# g_2000 <- ggplot(data=power_2000_new, aes(x=scenerio, y=`Average Power`, fill=method)) +
#   geom_bar(stat="identity", position=position_dodge())+
#   scale_fill_brewer(palette="Paired")+
#   #ggtitle(paste0("m : ",power_2000$m[1],", pi_m1 : ",power_2000$m1[1]/(power_2000$m[1] - power_2000$m3[1]),
#   #               ", q  :", power_2000$q[1])) + ylab("Average power") +
#   ggtitle("m = 2000") + 
#   theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),legend.text=element_text(size=14,face="bold.italic"))
# g_2000 + ggsave("power_2000.png")
# 
# ggplot(data=FDR_800_new, aes(x=scenerio, y=`Average Power`, fill=method)) +
#   geom_bar(stat="identity", position=position_dodge())+
#   scale_fill_brewer(palette="Paired")+
#   #ggtitle(paste0("m : ",FDR_800$m[1],", pi_m1 : ",FDR_800$m1[1]/(FDR_800$m[1] - FDR_800$m3[1]),
#   #               ", q  :", FDR_800$q[1])) + ylab("Average power") +
#   ggtitle("m = 800") +
#   theme(plot.title = element_text(color="black", size=14, face="bold.italic", hjust = 0.5),legend.text=element_text(size=14,face="bold.italic"))



fdr_2000 <- FDR_res[FDR_res$m == 2000 & FDR_res$q == 0.25 & 
                      FDR_res$m1/(FDR_res$m -FDR_res$m3) ==0.2,]


power_res[power_res$q == 0.4,]
