################################################################################
#
# SIMULATION OUTPUT
#
################################################################################
# Takes output files from the simulation and creates data tables and graphs
# Will automatically create tables/boxplot for the case of beta1=0
# Tables/plots for other scenarios commented out below

require(xtable)
require(reshape2)
require(ggplot2)
require(wesanderson)

## true beta2 (ie for temperature)
beta.truth.temp <- 0.017
beta.truth.vec.temp <- rep(beta.truth.temp,9)
beta.truth.mat.temp <- matrix(beta.truth.temp,nrow <- 2000,ncol=9)

## function to make results table
make_tab <- function(simest,simSE,truebeta.vec,truebeta.mat,cap="Simulation Results for Beta=0"){
  est.exp <- sprintf("%.4f",round(apply(simest,2,mean,na.rm=TRUE),4))
  bias.exp <-  sprintf("%.4f",round(apply(simest,2,mean,na.rm=TRUE)-truebeta.vec,4))
  per.bias.exp <-  sprintf("%.1f",round(100*(apply(simest,2,mean,na.rm=TRUE)-truebeta.vec)/truebeta.vec,1))
  SE.exp <-  sprintf("%.4f",round(apply(simSE,2,mean,na.rm=TRUE),4))
  SD.exp <- sprintf("%.4f",round(apply(simest,2,sd,na.rm=TRUE),4))
  RMSE.exp <-  sprintf("%.4f",round(sqrt(apply( (simest-truebeta.mat)^2,2,mean,na.rm=TRUE)),4))
  cvg.exp <-  sprintf("%.1f",round(100*apply( simest-1.96*simSE<truebeta.mat & simest+1.96*simSE>truebeta.mat,2,mean,na.rm=TRUE),1))
  
  op.char <- rbind(est.exp,bias.exp,per.bias.exp,SE.exp,SD.exp,RMSE.exp,cvg.exp)
  rownames(op.char) <- c('Est.','Bias','% Bias','SE','SD', 'RMSE','95% Coverage')
  colnames(op.char) <- c('True','Error','CS','RC','ML','B-Error','B-CS','B-RC','B-ML')
  xtable(op.char,caption=cap,align=c('c',rep('r',9)))
  return(op.char)
}
## function to make boxplot
make_boxplot <- function(simest,truebeta,limplus=0.05){
  boxdata <- melt(simest,na.rm=TRUE)
  boxdata$clr <- 0
  boxdata$clr[boxdata$variable=="Error" | boxdata$variable=="B-Error"] <- 1
  boxdata$clr[boxdata$variable=="CS" | boxdata$variable=="B-CS"] <- 2
  boxdata$clr[boxdata$variable=="RC" | boxdata$variable=="B-RC"] <- 3
  boxdata$clr[boxdata$variable=="ML" | boxdata$variable=="B-ML"] <- 4
  #boxdata.03$fade <- 1;  boxdata.03$fade[boxdata.03$variable=="B.Error" | boxdata.03$variable=="B.CS" | boxdata.03$variable=="B.RC"| boxdata.03$variable=="B.ML"] <-0.99 
  box <- ggplot(boxdata,aes(factor(variable),value))+
    geom_boxplot(aes(),alpha=0.8)+
    #geom_boxplot(aes(fill = factor(clr)),alpha=0.8)+
    #scale_fill_manual(values=wes_palette(n=5, name="Darjeeling")[c(1,3,4,5,2)],guide="none")+
    geom_hline(aes(yintercept=truebeta))+
    theme_bw() +
    theme(plot.background = element_blank(),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank()) +
    #ggtitle("Beta=0")+
    coord_cartesian(ylim=c(truebeta+limplus,truebeta-limplus))+
    labs(x="",y="Estimate")
  return(box)#box
}
##########################################
## Results for Simulation with beta=0 ##
##########################################

beta.truth <- 0
beta.truth.vec <- rep(beta.truth,9)
beta.truth.mat <- matrix(beta.truth,nrow <- 2000,ncol=9)

## LOAD FULL OUTPUT MATRIX
sim_full <- read.table("../Results/sim_par0_results_local.txt", header=TRUE, quote="\"")
## EXPOSURE
sim_est <- sim_full[,1:9]
sim_SE <- sim_full[,10:18]
## TEMPERATURE:
sim_temp_est <- sim_full[,19:27]
sim_temp_SE <- sim_full[,28:36]
## rename for plots
colnames(sim_est) <- colnames(sim_SE) <- colnames(sim_temp_est) <- colnames(sim_temp_SE) <- c('True','Error','CS','RC','ML','B-Error','B-CS','B-RC','B-ML')


## Table of Operating Characteristics ##
op.char.0 <- make_tab(sim_est,sim_SE,beta.truth.vec,beta.truth.mat,cap="Results for Beta=0")
xtable(op.char.0,caption="Simulation Results",align=c('c',rep('r',9)))
op.char.0.temp <- make_tab(sim_temp_est,sim_temp_SE,beta.truth.vec.temp,beta.truth.mat.temp,cap="Temperature Results for Beta=0")



## Boxplots ##
box.0 <- make_boxplot(sim_est,beta.truth,limplus=0.05)
box.0
ggsave("../Results/boxplot0.pdf",width=6,height=5)
box.0.temp <- make_boxplot(sim_temp_est,beta.truth.temp,limplus=0.05)
box.0.temp
ggsave("../Results/boxplot0_temp.pdf",width=6,height=5)




# ##########################################
# ## Results for Simulation with beta=0.015 ##
# ##########################################
# 
# beta.truth <- 0.015
# beta.truth.vec <- rep(beta.truth,9)
# beta.truth.mat <- matrix(beta.truth,nrow <- 2000,ncol=9)
# 
# ## LOAD FULL OUTPUT MATRIX
# sim_full <- read.table("../Results/sim_par15_results_local.txt", header=TRUE, quote="\"")
# 
# 
# ## EXPOSURE
# sim_est <- sim_full[,1:9]
# sim_SE <- sim_full[,10:18]
# ## TEMPERATURE:
# sim_temp_est <- sim_full[,19:27]
# sim_temp_SE <- sim_full[,28:36]
# ## rename for plots
# colnames(sim_est) <- colnames(sim_SE) <- colnames(sim_temp_est) <- colnames(sim_temp_SE) <- c('True','Error','CS','RC','ML','B-Error','B-CS','B-RC','B-ML')
# 
# 
# ## Table of Operating Characteristics ##
# op.char.015 <- make_tab(sim_est,sim_SE,beta.truth.vec,beta.truth.mat,cap="Results for Beta=0.015")
# op.char.015.temp <- make_tab(sim_temp_est,sim_temp_SE,beta.truth.vec.temp,beta.truth.mat.temp,cap="Temperature Results for Beta=0.015")
# 
# ## Boxplots ##
# box.015 <- make_boxplot(sim_est,beta.truth,limplus=0.05)
# box.015
# ggsave("../Results/boxplot015.pdf",width=6,height=5)
# box.015.temp <- make_boxplot(sim_temp_est,beta.truth.temp,limplus=0.05)
# box.015.temp
# ggsave("../Results/boxplot015_temp.pdf",width=6,height=5)
# 
# 
# ##########################################
# ## Results for Simulation with beta=0.03 ##
# ##########################################
# 
# beta.truth <- 0.03
# beta.truth.vec <- rep(beta.truth,9)
# beta.truth.mat <- matrix(beta.truth,nrow <- 2000,ncol=9)
# 
# ## LOAD FULL OUTPUT MATRIX
# sim_full <- read.table("../Results/sim_par3_results_local.txt", header=TRUE, quote="\"")
# 
# 
# ## EXPOSURE
# sim_est <- sim_full[,1:9]
# sim_SE <- sim_full[,10:18]
# ## TEMPERATURE:
# sim_temp_est <- sim_full[,19:27]
# sim_temp_SE <- sim_full[,28:36]
# ## rename for plots
# colnames(sim_est) <- colnames(sim_SE) <- colnames(sim_temp_est) <- colnames(sim_temp_SE) <- c('True','Error','CS','RC','ML','B-Error','B-CS','B-RC','B-ML')
# 
# 
# ## Table of Operating Characteristics ##
# op.char.03 <- make_tab(sim_est,sim_SE,beta.truth.vec,beta.truth.mat,cap="Results for Beta=0.03")
# op.char.03.temp <- make_tab(sim_temp_est,sim_temp_SE,beta.truth.vec.temp,beta.truth.mat.temp,cap="Temperature Results for Beta=0.03")
# 
# ## Boxplots ##
# box.03 <- make_boxplot(sim_est,beta.truth,limplus=0.05)
# box.03
# ggsave("../Results/boxplot03.pdf",width=6,height=5)
# box.03.temp <- make_boxplot(sim_temp_est,beta.truth.temp,limplus=0.05)
# box.03.temp
# ggsave("../Results/boxplot03_temp.pdf",width=6,height=5)
# 
# 
# ##########################################
# ## Results for Simulation with beta=0.045 ##
# ##########################################
# 
# beta.truth <- 0.045
# beta.truth.vec <- rep(beta.truth,9)
# beta.truth.mat <- matrix(beta.truth,nrow <- 2000,ncol=9)
# 
# ## LOAD FULL OUTPUT MATRIX
# sim_full <- read.table("../Results/sim_par45_results_local.txt", header=TRUE, quote="\"")
# 
# 
# ## EXPOSURE
# sim_est <- sim_full[,1:9]
# sim_SE <- sim_full[,10:18]
# ## TEMPERATURE:
# sim_temp_est <- sim_full[,19:27]
# sim_temp_SE <- sim_full[,28:36]
# ## rename for plots
# colnames(sim_est) <- colnames(sim_SE) <- colnames(sim_temp_est) <- colnames(sim_temp_SE) <- c('True','Error','CS','RC','ML','B-Error','B-CS','B-RC','B-ML')
# 
# ## Table of Operating Characteristics ##
# op.char.045 <- make_tab(sim_est,sim_SE,beta.truth.vec,beta.truth.mat,cap="Results for Beta=0.045")
# op.char.045.temp <- make_tab(sim_temp_est,sim_temp_SE,beta.truth.vec.temp,beta.truth.mat.temp,cap="Temperature Results for Beta=0.045")
# 
# 
# ## Boxplots ##
# box.045 <- make_boxplot(sim_est,beta.truth,limplus=0.05)
# box.045
# ggsave("../Results/boxplot045.pdf",width=6,height=5)
# box.045.temp <- make_boxplot(sim_temp_est,beta.truth.temp,limplus=0.05)
# box.045.temp
# ggsave("../Results/boxplot045_temp.pdf",width=6,height=5)
# 
# 
# 
# 
# 
# 
# ## MASTER TABLE
# 
# xtable(rbind(op.char.0,op.char.015,op.char.03,op.char.045),caption="Simulation Results",align=c('c',rep('r',9)))
# ## This is easier for rownames:
# xtable(op.char.0,caption="Simulation Results",align=c('c',rep('r',9)))
# xtable(op.char.015,caption="Simulation Results",align=c('c',rep('r',9)))
# xtable(op.char.03,caption="Simulation Results",align=c('c',rep('r',9)))
# xtable(op.char.045,caption="Simulation Results",align=c('c',rep('r',9)))
# 
# 
# ## same for temperature effects
# xtable(rbind(op.char.0.temp,op.char.015.temp,op.char.03.temp,op.char.045.temp),caption="Simulation Results for Temperature",align=c('c',rep('r',9)))
# ## This is easier for rownames:
# xtable(op.char.0.temp,caption="Simulation Results for Temperature",align=c('c',rep('r',9)))
# xtable(op.char.015.temp,caption="Simulation Results for Temperature",align=c('c',rep('r',9)))
# xtable(op.char.03.temp,caption="Simulation Results for Temperature",align=c('c',rep('r',9)))
# xtable(op.char.045.temp,caption="Simulation Results for Temperature",align=c('c',rep('r',9)))
# 
# 









