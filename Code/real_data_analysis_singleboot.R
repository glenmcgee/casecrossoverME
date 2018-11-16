################################################################################
#
# DATA ANALYSIS WITH ONLY PARAMETRIC BOOTSTRAP 
# PARALLELIZED
#
################################################################################
## Updated: 11/16/2018
  ## includes PM (avg) and AMBIENT TEMPERATURE (natural cubic splines) and pressure (continuous)
  ## some missing data on case times -- using admit time as error for now since it has much less missingness
  ## note the functions here and MLmatrix are slightly different from sim (to handle more covariates)
  ## Not creating MLmatrix here since it is slow
## NOTE: MLmatrix, along with pm.seq and controlz all begin at time 1
  ## X.mat begins at time 1002 (to generate data for times that can then be averaged back in time)
  ## sim.dat.new() then adds the necessary 1001 to the event times to get the proper covariates from MLmatrix
  ## raw times for true data and simulated data are all lined up with MLmatrix
  ## so no need to add 1001 to the times in any functions (unlike the sim)
    


## set TRUE to run locally
runLOCAL <- TRUE
## set FALSE to run job array on cluster via sbatch (requires changing file paths in "Cluster/Local" section below)

####################
##   Parameters   ##
####################
BB <- 200 ## BB parametric bootstrap resamples
SS <- 300 ## SS (validation sample size)
# JJ <- 200 ## number of jobs in array

#######################
##   Load packages   ##
#######################
require(survival)
require(foreign)
require(splines)
# library(foreach)
# library(doMPI)

#######################
##   Cluster/Local   ##
#######################
## set up parallel
if(runLOCAL==TRUE){
  path <- "../Results/" ## path for results
  datpath <- "../Data/" ## path for data
  loopvec <- 1:BB ## run entire loop
  suffix <- "_local" ## output file has a suffix to distinguish from other results files in dropbox
}else{
  path <- "../Results/" ## path for results
  datpath <- "../Data/" ## path for data
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  iter_no <- as.integer(args[1])  ## here there is only one argument, the iteration number
  loopvec <- ((BB/JJ)*(iter_no-1) +1):((BB/JJ)*iter_no) ## portion of the loop to run
  suffix <- paste0("_iter",iter_no) ## append output file names with the iteration number, to be combined later
}

# ##########################
# ##   Parallel cluster   ##
# ##########################
# ## make cluster for parallel bootstrap
# numcore=100
# #cl <- makeCluster(numcore-1)
# #registerDoParallel(cl)
# cl <- startMPIcluster(count=numcore-1)
# registerDoMPI(cl)

###################
##   Load Data   ##
###################
## load data
## Cases <- read.csv(paste0(datpath,"BrentCases_NF.csv")) ## DATA NOT PUBLICLY AVAILABLE
covs <- read.csv(paste0(datpath,"brentpmhourly.csv"))
MLmatrix <- read.table(paste0(datpath,"MLmatrix_realdata.txt"), header=TRUE, quote="\"")


####################
##  Process Data  ##
####################
pm.seq <- covs$PM25_pred_ma24_l0 ## exposure (PM2.5)
## replace missing values with average of preceding and following values 
## missing values occur in groups of 10-20ish
pm.seq[(1672-1):(1687+1)] <- seq(pm.seq[1671],pm.seq[1688],length.out=length((1672-1):(1687+1)))
pm.seq[(8519-1):(8540+1)] <- seq(pm.seq[8518],pm.seq[8541],length.out=length((8519-1):(8540+1)))
pm.seq[(16703-1):(16738+1)] <- seq(pm.seq[16702],pm.seq[16739],length.out=length((16703-1):(16738+1)))
pm.seq[(35388-1):(39657+1)] <- seq(pm.seq[35387],pm.seq[39658],length.out=length((35388-1):(39657+1)))
pm.seq[(40675-1):(40696+1)] <- seq(pm.seq[40674],pm.seq[40697],length.out=length((40675-1):(40696+1)))
pm.seq[(50636-1):(50658+1)] <- seq(pm.seq[50635],pm.seq[50659],length.out=length((50636-1):(50658+1)))
pm.seq[(51167-1):(51180+1)] <- seq(pm.seq[51166],pm.seq[51181],length.out=length((51167-1):(51180+1)))
pm.seq[(65995-1):(66000+1)] <- seq(pm.seq[65994],pm.seq[66001],length.out=length((65995-1):(66000+1)))
pm.seq[(78935-1):(78946+1)] <- seq(pm.seq[78934],pm.seq[78947],length.out=length((78935-1):(78946+1)))

temp <- covs$tempC_ma24_l0 ## covariate (temperature)
## replace missing values with average of preceding and following values 
temp[(12942-1):(12947+1)] <- seq(temp[12941],temp[12948],length.out=length((12942-1):(12947+1)))
temp[(19289-1):(19311+1)] <- seq(temp[19288],temp[19312],length.out=length((19289-1):(19311+1)))
temp[(23263-1):(23310+1)] <- seq(temp[23262],temp[23311],length.out=length((23263-1):(23310+1)))
temp[(45202-1):(45218+1)] <- seq(temp[45201],temp[45219],length.out=length((45202-1):(45218+1)))
temp[(91918-1):(91937+1)] <- seq(temp[91917],temp[91938],length.out=length((91918-1):(91937+1)))

dpt <- covs$dptC_ma24_l0
dpt[(12942-1):(12947+1)] <- seq(dpt[12941],dpt[12948],length.out=length((12942-1):(12947+1)))
dpt[(19289-1):(19311+1)] <- seq(dpt[19288],dpt[19312],length.out=length((19289-1):(19311+1)))
dpt[(23263-1):(23310+1)] <- seq(dpt[23262],dpt[23311],length.out=length((23263-1):(23310+1)))
dpt[(45202-1):(45218+1)] <- seq(dpt[45201],dpt[45219],length.out=length((45202-1):(45218+1)))
dpt[(91918-1):(91937+1)] <- seq(dpt[91917],dpt[91938],length.out=length((91918-1):(91937+1)))

press <- covs$pres_ma24_l0
# adding 2160 cuz were no longer taking covs <- covs[2161:nrow(covs)]
press[2160+(10782-1):(10787+1)] <- seq(press[2160+10781],press[2160+10788],length.out=length(2160+(10782-1):(10787+1)))
press[2160+(17131-1):(17151+1)] <- seq(press[2160+17130],press[2160+17152],length.out=length(2160+(17131-1):(17151+1)))
press[2160+(21103-1):(21150+1)] <- seq(press[2160+21102],press[2160+21151],length.out=length(2160+(21103-1):(21150+1)))
press[2160+(43042-1):(43058+1)] <- seq(press[2160+43041],press[2160+43059],length.out=length(2160+(43042-1):(43058+1)))
press[2160+(89758-1):(89777+1)] <- seq(press[2160+89757],press[2160+89778],length.out=length(2160+(89758-1):(89777+1)))
# replace starting NA's because they screw up matching
press[1:17] <- press[18]



### Case-time data:
Cases$OnsetDateTime_Best <- as.character(Cases$OnsetDateTime_Best)
## Impute date+time where missing (but date and time data exist)
Cases$OnsetDateTime_Best[1132] <- "07MAY05:11:00:00"
Cases$OnsetDateTime_Best[1172] <- "25MAR05:09:00:00"
Cases$OnsetDateTime_Best[1194] <- "08JAN06:08:00:00"
Cases$OnsetDateTime_Best[1265] <- "04APR05:18:00:00"

## Greg's paper writes that patients for whom the DATE is known but not the actual TIME, they imputed 9AM. 
## One patient still does not have a TIME but has a DATE:
Cases$OnsetDateTime_Best[1093] <- "14JUL05:09:00:00"

## Currently only consider times when true AND admit exist
Cases <- Cases[(as.character(Cases$OnsetDateTime_Best)!="") & (as.character(Cases$Admit_datetime)!=""),]

## get time index
time_true <- c()
time_admit <- c()
time_ED <- c()
for (jj in 1:nrow(Cases)){
  if (!is.na(as.character(Cases$OnsetDateTime_Best[jj]))){
    time_true <- c(time_true,which(as.character(covs$datetime)==as.character(Cases$OnsetDateTime_Best[jj])))
  }
  if (!is.na(as.character(Cases$Admit_datetime[jj]))){
    time_admit <- c(time_admit,which(as.character(covs$datetime)==as.character(Cases$Admit_datetime[jj])))
  }
  if (!is.na(as.character(Cases$ED_datetime[jj]))){
    time_ED <- c(time_ED,which(as.character(covs$datetime)==as.character(Cases$ED_datetime[jj])))
  }
}


## collect true delay times
true_delays <- time_admit-time_true



##########################
##  Create data matrix  ##
##########################
## time scale
hours.day <- 24
days.year <- 365
#years <- 5
years <- 5*2

days.month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
days.in.month <- c(seq(1,31,1),seq(1,28,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,30,1),seq(1,31,1))

days.month.leap <- c(31,29,31,30,31,30,31,31,30,31,30,31)
days.in.month.leap <- c(seq(1,31,1),seq(1,29,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,30,1),seq(1,31,1))

total.hours <- hours.day*days.year*years+hours.day*3  #3 leap years (00,04,08)
total.days <- days.year*years+3


month.seq <- c(rep(1,days.month[1]),rep(2,days.month[2]),rep(3,days.month[3]),rep(4,days.month[4]),rep(5,days.month[5]),rep(6,days.month[6]),
               rep(7,days.month[7]),rep(8,days.month[8]),rep(9,days.month[9]),rep(10,days.month[10]),rep(11,days.month[11]),rep(12,days.month[12]))
month.seq.leap <- c(rep(1,days.month[1]),rep(2,29),rep(3,days.month[3]),rep(4,days.month[4]),rep(5,days.month[5]),rep(6,days.month[6]),
                    rep(7,days.month[7]),rep(8,days.month[8]),rep(9,days.month[9]),rep(10,days.month[10]),rep(11,days.month[11]),rep(12,days.month[12]))
month.seq <- rep(month.seq, each=24)
month.seq.leap <- rep(month.seq.leap, each=24)


month.seq <- c(month.seq, month.seq.leap+12, month.seq+24, month.seq+36, 
               month.seq+48, month.seq.leap+60, month.seq+72, month.seq+84,
               month.seq+96, month.seq.leap+108, month.seq+120, month.seq+132)
hours.in.day <- rep(seq(1,24,1),total.days+365)
hours.seq <- seq(1,total.hours+(365*24))
day.seq <- (floor((hours.seq-1)/24))+1
week.seq <- floor((day.seq-1) / 7)+1
day.in.week <- rep(seq(1,7,1),each=24,length.out=(total.days+365)*24)

## delay times
true_delays <- true_delays[true_delays>0]
tab<-table(true_delays)    
frq<-tab; del<-as.numeric(names(tab)) ## del is only the unique values because well multiply by...
phi<-as.vector(frq/sum(frq)) ## delay distribution


## covariate time series
intercept <- rep(1,total.hours)
temp.spline <- matrix(ns(temp,df=3),ncol=3)
dpt.spline <- matrix(ns(dpt,df=3),ncol=3)

## covariate time series starting at time 1002 (for generating data)
press.post1001 <- press[1:(1001+total.hours+1000)]
temp.spline.post1001 <- temp.spline[1:(1001+total.hours+1000),] 
dpt.spline.post1001 <- dpt.spline[1:(1001+total.hours+1000),] 
X.mat <- cbind(intercept, pm.seq[1002:(1001+total.hours)],press.post1001[1002:(1001+total.hours)],temp.spline.post1001[1002:(1001+total.hours),],dpt.spline.post1001[1002:(1001+total.hours),])




## import MLmatrix just for now (quicker)
## this is a file that includes PM and temperature, and its faster to load than to wait 30 minutes for it to populate
MLmatrix <- data.matrix(MLmatrix)

## covariate series used to create MLmatrix
control1_0 <- press[1:(total.hours+2001)]
control2_0<- temp.spline[1:(total.hours+2001),1];control3_0<- temp.spline[1:(total.hours+2001),2];control4_0<- temp.spline[1:(total.hours+2001),3]
control5_0<- dpt.spline[1:(total.hours+2001),1];control6_0<- dpt.spline[1:(total.hours+2001),2];control7_0<- dpt.spline[1:(total.hours+2001),3]
controlz <- data.matrix(cbind(control1_0,control2_0,control3_0,control4_0,control5_0,control6_0,control7_0))





# ##############
# ### MATRIX ###
# ##############
# ## Code to create the file: "MLmatrix_realdata.txt"
# # total.hours indexes row by hour (ie same as the row index)
# hours.5years <- c(hours.seq[1:(total.hours+2001)])   
# 
# # n.matches dictates how many (4 or 5) observations are in a matched set for a specific event time (including the actual event/case)
# n.days.in.month <- c(rep(31,31),rep(28,28),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31))
# which.day <- (days.in.month-1)%%7+1     # indexes days 1:7 only within a month (unlike days.in.week)
# n.matches <- 4+(which.day+7*4<=n.days.in.month)
# 
# n.days.in.month.leap <- c(rep(31,31),rep(29,29),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31))
# which.day.leap <- (days.in.month.leap-1)%%7+1     # indexes days 1:7 only within a month (unlike days.in.week)
# n.matches.leap <- 4+(which.day.leap+7*4<=n.days.in.month.leap)
# 
# n.matches <- c(rep(n.matches,each=24),rep(n.matches.leap,each=24),rep(n.matches,each=24),rep(n.matches,each=24),
#                rep(n.matches,each=24),rep(n.matches.leap,each=24),rep(n.matches,each=24),rep(n.matches,each=24),
#                rep(n.matches,each=24),rep(n.matches.leap,each=24))
# 
# n.matches <- c(n.matches,head(n.matches,2001))    # add 1000 to end so we can delay to after total.hours
# 
# 
# 
# # which week of the month is it
# which.week <- 1+(days.in.month>7)+(days.in.month>14)+(days.in.month>21)+(days.in.month>28)
# which.week.leap <- 1+(days.in.month.leap>7)+(days.in.month.leap>14)+(days.in.month.leap>21)+(days.in.month.leap>28)
# 
# which.week <- c(rep(which.week,each=24),rep(which.week.leap,each=24),rep(which.week,each=24),rep(which.week,each=24),
#                 rep(which.week,each=24),rep(which.week.leap,each=24),rep(which.week,each=24),rep(which.week,each=24),
#                 rep(which.week,each=24),rep(which.week.leap,each=24))
# 
# which.week <- c(which.week,head(which.week,2001))    # add 1001 to beginning and 1000 to end
# 
# 
# 
# # exposure for case
# exp_0 <- pm.seq[1:(total.hours+2001)]
# # control variables/confounders for case
# #control1_0 <- control1.full    # full because we have added on the 1001 rows at the beginning
# control1_0 <- press[1:(total.hours+2001)]
# control2_0<- temp.spline[1:(total.hours+2001),1]; control3_0<- temp.spline[1:(total.hours+2001),2];control4_0<- temp.spline[1:(total.hours+2001),3];
# control5_0<- dpt.spline[1:(total.hours+2001),1]; control6_0<- dpt.spline[1:(total.hours+2001),2];control7_0<- dpt.spline[1:(total.hours+2001),3];
# 
# 
# controlz <- data.matrix(cbind(control1_0,control2_0,control3_0,control4_0,control5_0,control6_0,control7_0))
# 
# # hours (to be used for sim data)
# hours_0 <- hours.5years
# 
# 
# # initialize 
# exp_1 <- exp_2 <- exp_3 <- exp_4 <- rep(-99,total.hours+2001)                        # use -99 to avoid loops later on?
# control1_1 <- control1_2 <- control1_3 <- control1_4 <- rep(-99,total.hours+2001)
# control2_1 <- control2_2 <- control2_3 <- control2_4 <- rep(-99,total.hours+2001)
# control3_1 <- control3_2 <- control3_3 <- control3_4 <- rep(-99,total.hours+2001)
# control4_1 <- control4_2 <- control4_3 <- control4_4 <- rep(-99,total.hours+2001)
# control5_1 <- control5_2 <- control5_3 <- control5_4 <- rep(-99,total.hours+2001)
# control6_1 <- control6_2 <- control6_3 <- control6_4 <- rep(-99,total.hours+2001)
# control7_1 <- control7_2 <- control7_3 <- control7_4 <- rep(-99,total.hours+2001)
# hours_1 <- hours_2 <- hours_3 <- hours_4 <- rep(-99,total.hours+2001)
# 
# # initialize matrix
# MLmatrix <- data.frame(cbind(hours.5years,n.matches,which.week,
#                              exp_0,control1_0,control2_0,control3_0,control4_0,control5_0,control6_0,control7_0,
#                              exp_1,control1_1,control2_1,control3_1,control4_1,control5_1,control6_1,control7_1,
#                              exp_2,control1_2,control2_2,control3_2,control4_2,control5_2,control6_2,control7_2,
#                              exp_3,control1_3,control2_3,control3_3,control4_3,control5_3,control6_3,control7_3,
#                              exp_4,control1_4,control2_4,control3_4,control4_4,control5_4,control6_4,control7_4,
#                              hours_0,hours_1,hours_2,hours_3,hours_4))
# 
# # restrict everything to the 5 years
# month.seq.5years <- month.seq[1:total.hours]  
# month.seq.5years <- c(month.seq.5years,head(month.seq.5years,2001)+120)   # add 1000 to end
# 
# day.in.week.5years <- day.in.week[1:total.hours]
# day.in.week.5years <- c(day.in.week.5years,head(rep(seq(1,7),300),2001))       # add 1001 to beginning and 1000 to end
# 
# hours.in.day.5years <- hours.in.day[1:total.hours]
# hours.in.day.5years <- c(hours.in.day.5years,head(rep(seq(1,24),100),2001))       # add 1001 to beginning and 1000 to end
# 
# 
# for (jj in 1:(total.hours+2001)){
#   # subset matched data; exclude case since its already in there
#   matches <- MLmatrix[month.seq.5years==month.seq.5years[jj] & day.in.week.5years==day.in.week.5years[jj] & hours.in.day.5years==hours.in.day.5years[jj] & hours.5years!=jj,]
#   #first match
#   MLmatrix[jj,11+(1:8)] <- matches[1,4:11]
#   MLmatrix[jj,45] <- matches[1,1]      # to have a time column for sim.data
#   # second match
#   MLmatrix[jj,19+(1:8)] <- matches[2,4:11]
#   MLmatrix[jj,46] <- matches[2,1]
#   # third match
#   MLmatrix[jj,27+(1:8)] <- matches[3,4:11]
#   MLmatrix[jj,47] <- matches[3,1]
#   # possible 4th match
#   if ( nrow(matches)==4){
#     MLmatrix[jj,35+(1:8)] <- matches[4,4:11]
#     MLmatrix[jj,48] <- matches[4,1]
#     
#   }  
# }
# 
# match5 <- as.numeric(MLmatrix$n.matches==5)     ### indicator of whether or not to include last match in denominator calculation
# MLmatrix <- cbind(MLmatrix,match5)              ### add indicator to make code run faster
# # write.table(x=MLmatrix, file='~/MeasErr/Results/MLmatrix_realdata.txt', row.names=F, col.names=T, quote=F)
# MLmatrix <- data.matrix(MLmatrix)







############################
### FUNCTION DEFINITIONS ###
############################
## pretty sure we dont have to shift here since the numbering counts from Jan1st so in the data we have like 2000 etc.
### negative logLikelihood (to be maximized)
neg.logL.MLnew<-function(b,dataframe,del,phi){
  subj.times <- dataframe$time.error[dataframe$case.vece==1]     ## already shifted by 1001
  res<-0
  #ct.nan<-0
  for(a in subj.times){
    
    subj.data <- MLmatrix[a-del,]
    
    case.contrib <- exp(subj.data[,3+(1:8) ]%*%b )                            # numerator
    match.contrib <- case.contrib+
      exp(subj.data[,11+(1:8) ]%*%b )+
      exp(subj.data[,19+(1:8) ]%*%b )+
      exp(subj.data[,27+(1:8) ]%*%b )+
      exp(subj.data[,35+(1:8) ]%*%b )*subj.data[,ncol(subj.data)]     # multiply by indicator
    # (ie only include if it exists)
    vals <- case.contrib/match.contrib
    
    res <- res + log(sum(vals*phi))
    
    #if ( res.1=='NaN' | is.na(res.1)==TRUE ) ct.nan<-ct.nan+1
    #if ( res.1!='NaN' & is.na(res.1)!=TRUE ) res<-res+res.1
  }
  return ( -res )
}


### gradient (vectorized)
grad.MLnew<-function(b,dataframe,del,phi){
  subj.times <- dataframe$time.error[dataframe$case.vece==1]     ## already shifted by 1001
  res <- rep(0,length(b))
  #ct.nan<-0
  for(a in subj.times){
    subj.data <- MLmatrix[a-del,]
    
    tmp.1.1 <- c(exp(subj.data[,4+(0:7) ]%*%b ))
    tmp.1.2 <- subj.data[,4]*tmp.1.1
    tmp.g.1 <- subj.data[,4+(1:7)]*tmp.1.1
    # initialize (including first term from above)
    tmp.2.1 <- tmp.1.1
    tmp.2.2 <- tmp.1.2
    tmp.g.2 <- tmp.g.1
    #for (nn in 2:4){                               # loop over matches (depending on number of matches)
    for (nn in 1:3){
      tmp.nn.1 <- c(exp(subj.data[,4+8*nn+(0:7) ]%*%b ) )
      tmp.2.1 <- tmp.2.1+tmp.nn.1
      tmp.2.2 <- tmp.2.2+subj.data[,4+8*nn]*tmp.nn.1
      tmp.g.2 <- tmp.g.2+subj.data[,4+8*nn+(1:7)]*tmp.nn.1
    }
    
    # potential match
    tmp.5.1 <- c(exp(subj.data[,4+8*4+(0:7) ]%*%b )*subj.data[,ncol(subj.data)] )     # multiply by indicator
    tmp.2.1 <- tmp.2.1+tmp.5.1
    tmp.2.2 <- tmp.2.2+subj.data[,4+8*4]*tmp.5.1
    tmp.g.2 <- tmp.g.2+subj.data[,4+8*4+(1:7)]*tmp.5.1
    
    A.b <- sum((tmp.1.2/tmp.2.1)*phi)
    B.b <- sum((tmp.1.1/tmp.2.1*tmp.2.2/tmp.2.1)*phi)
    C <- sum((tmp.1.1/tmp.2.1)*phi)
    
    A.g<-colSums(tmp.g.1/tmp.2.1 * phi)
    B.g<-colSums(tmp.g.2/tmp.2.1*tmp.1.1/tmp.2.1 * phi)
    res1<-(A.b-B.b)/C
    res2<-(A.g-B.g)/C
    
    res <- res+c(res1,res2)
    #if ( res.1=='NaN' | is.na(res.1)==TRUE ) ct.nan<-ct.nan+1
    #if ( res.1!='NaN' & is.na(res.1)!=TRUE ) res<-res+res.1
  }
  return ( -res )
}

### Marginal likelihood approach  
fit_MargLik <- function(ddat,beta.init){
  
  ## beta.init<-cc.fit.error$coefficients ## starting values
  
  fit<-tryCatch(optim(beta.init,neg.logL.MLnew,grad.MLnew,dataframe=ddat$error,del=ddat$del,phi=ddat$phi,hessian=TRUE,control=list(maxit=2000)),error=function(a){NA}) # if hessian=TRUE, then hessian can be computed numerically
  if(sum(is.na(fit))>0){ ## in case of convergence error
    return(list(est=rep(NA,length(beta.init)),SE=rep(NA,length(beta.init))))
  }else{
    est <- fit$par ## extract ML estimates
    SE <- sqrt(diag(solve(fit$hessian))) ## extract SEs
    
    return(list(est=est,SE=SE))
  }  
  ## Note: using the numerical information matrix is only slightly less accurate (only noticeable at 6th decimal point)
  ## To compute actual information matrix (very slowly): ## SE <- infomat.MLnew(est,dataframe=ddat$error,del=ddat$del,phi=ddat$phi)   
  
  
}

### Fit regression calibration 
fit_regcal <- function(ddat,beta.init){
  
  data.cc.error <- ddat$error
  ref.time.error<-unique(data.cc.error$time.error)                    # unique error times
  
  pred.pm<-rep(0,length(ref.time.error))
  pred.cont<-matrix(0,nrow=length(ref.time.error),ncol=7)
  for(a in 1:length(ref.time.error)){                                 # loop over unique error times
    idx<-ref.time.error[a]-ddat$del                                   # index for possible exposure times (to be summed over)
    tmp.x<-pm.seq[idx]                                                # each possible exposure to sum over
    tmp.cont<-controlz[idx,] 
    pred.pm[a]<-sum(ddat$phi*tmp.x) #pred.pm[a]<-sum(ddat$phi*tmp.x)  # expected exposure 
    pred.cont[a,]<-ddat$phi%*%tmp.cont
    
  }
  exp.vecr<-rep(0,nrow(data.cc.error))          
  cont.vecr<-matrix(0,nrow=nrow(data.cc.error),ncol=7)      
  for(a in 1:nrow(data.cc.error)){
    exp.vecr[a]<-pred.pm[which(ref.time.error==data.cc.error$time.error[a])]    # assign predicted vals (do it this way because we found predicteds for unique times--> must distribute)
    cont.vecr[a,]<-pred.cont[which(ref.time.error==data.cc.error$time.error[a]),] #cont.vecr[a]<-pred.cont[which(ref.time.error==data.cc.error$time.error[a])]
    
  }
  data.cc.RC<-data.cc.error 
  data.cc.RC[,3]<-exp.vecr          # replace exp.vece with expected value exp.vecr -- it will still be called exp.vece in the data
  data.cc.RC[,5:11]<-cont.vecr
  fit <- clogit(case.vece ~ exp.vece+cont.vece1+cont.vece2+cont.vece3+cont.vece4+cont.vece5+cont.vece6+cont.vece7+strata(subj.vece), data=data.cc.RC)
  est <- fit$coef
  SE <- sqrt(diag(fit$var))
  
  return(list(est=est,SE=SE))
  
}



## get SIGMA.u for conditional score approach
get_sigma <- function(ddat,n.reps=50){
  sd.u <- NA
  while(sum(is.na(sd.u))>0){
    # estimate measurement error variance
    U <- U.cont1 <- U.cont2 <- U.cont3 <- U.cont4 <- U.cont5 <- U.cont6 <- U.cont7 <- NULL
    for(a in 1:n.reps){  
      
      ## sampling from validation sample of delay times (300 of them)
      delay2 <- sample(ddat$delays,length(ddat$hours.events),replace=T) 
      
      ## compute U 
      u<-pm.seq[hours.events.error-delay2] - pm.seq[hours.events.error]            ## minus since it should be U not U
      U<-cbind(U,u)
      
      u.cont1<-controlz[hours.events.error-delay2,1]-controlz[hours.events.error,1]
      U.cont1<-cbind(U.cont1,u.cont1)
      
      u.cont2<-controlz[hours.events.error-delay2,2]-controlz[hours.events.error,2]
      U.cont2<-cbind(U.cont2,u.cont2)
      
      u.cont3<-controlz[hours.events.error-delay2,3]-controlz[hours.events.error,3]
      U.cont3<-cbind(U.cont3,u.cont3)
      
      u.cont4<-controlz[hours.events.error-delay2,4]-controlz[hours.events.error,4]
      U.cont4<-cbind(U.cont4,u.cont4)
      
      u.cont5<-controlz[hours.events.error-delay2,5]-controlz[hours.events.error,5]
      U.cont5<-cbind(U.cont5,u.cont5)
      
      u.cont6<-controlz[hours.events.error-delay2,6]-controlz[hours.events.error,6]
      U.cont6<-cbind(U.cont6,u.cont6)
      
      u.cont7<-controlz[hours.events.error-delay2,7]-controlz[hours.events.error,7]
      U.cont7<-cbind(U.cont7,u.cont7)
      
      ## collect data
      U<-cbind(U,u)
      U.cont1<-cbind(U.cont1,u.cont1)
      U.cont2<-cbind(U.cont2,u.cont2)
      U.cont3<-cbind(U.cont3,u.cont3)
      U.cont4<-cbind(U.cont4,u.cont4)
      U.cont5<-cbind(U.cont5,u.cont5)
      U.cont6<-cbind(U.cont6,u.cont6)
      U.cont7<-cbind(U.cont7,u.cont7)
    }
    
    s2.u11<-mean(apply(U,1,var))
    s2.u22<-mean(apply(U.cont1,1,var))
    s2.u33<-mean(apply(U.cont2,1,var))
    s2.u44<-mean(apply(U.cont3,1,var))
    s2.u55<-mean(apply(U.cont4,1,var))
    s2.u66<-mean(apply(U.cont5,1,var))
    s2.u77<-mean(apply(U.cont6,1,var))
    s2.u88<-mean(apply(U.cont7,1,var))
    
    s2.u12<-mean(diag((U-apply(U,1,mean))%*%t((U.cont1-apply(U.cont1,1,mean)))/(n.reps-1)))
    s2.u13<-mean(diag((U-apply(U,1,mean))%*%t((U.cont2-apply(U.cont2,1,mean)))/(n.reps-1)))
    s2.u14<-mean(diag((U-apply(U,1,mean))%*%t((U.cont3-apply(U.cont3,1,mean)))/(n.reps-1)))
    s2.u15<-mean(diag((U-apply(U,1,mean))%*%t((U.cont4-apply(U.cont4,1,mean)))/(n.reps-1)))
    s2.u16<-mean(diag((U-apply(U,1,mean))%*%t((U.cont5-apply(U.cont5,1,mean)))/(n.reps-1)))
    s2.u17<-mean(diag((U-apply(U,1,mean))%*%t((U.cont6-apply(U.cont6,1,mean)))/(n.reps-1)))
    s2.u18<-mean(diag((U-apply(U,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))#
    
    s2.u23<-mean(diag((U.cont1-apply(U.cont1,1,mean))%*%t((U.cont2-apply(U.cont2,1,mean)))/(n.reps-1)))
    s2.u24<-mean(diag((U.cont1-apply(U.cont1,1,mean))%*%t((U.cont3-apply(U.cont3,1,mean)))/(n.reps-1)))
    s2.u25<-mean(diag((U.cont1-apply(U.cont1,1,mean))%*%t((U.cont4-apply(U.cont4,1,mean)))/(n.reps-1)))
    s2.u26<-mean(diag((U.cont1-apply(U.cont1,1,mean))%*%t((U.cont5-apply(U.cont5,1,mean)))/(n.reps-1)))
    s2.u27<-mean(diag((U.cont1-apply(U.cont1,1,mean))%*%t((U.cont6-apply(U.cont6,1,mean)))/(n.reps-1)))
    s2.u28<-mean(diag((U.cont1-apply(U.cont1,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))
    
    s2.u34<-mean(diag((U.cont2-apply(U.cont2,1,mean))%*%t((U.cont3-apply(U.cont3,1,mean)))/(n.reps-1)))
    s2.u35<-mean(diag((U.cont2-apply(U.cont2,1,mean))%*%t((U.cont4-apply(U.cont4,1,mean)))/(n.reps-1)))
    s2.u36<-mean(diag((U.cont2-apply(U.cont2,1,mean))%*%t((U.cont5-apply(U.cont5,1,mean)))/(n.reps-1)))
    s2.u37<-mean(diag((U.cont2-apply(U.cont2,1,mean))%*%t((U.cont6-apply(U.cont6,1,mean)))/(n.reps-1)))
    s2.u38<-mean(diag((U.cont2-apply(U.cont2,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))
    
    s2.u45<-mean(diag((U.cont3-apply(U.cont3,1,mean))%*%t((U.cont4-apply(U.cont4,1,mean)))/(n.reps-1)))
    s2.u46<-mean(diag((U.cont3-apply(U.cont3,1,mean))%*%t((U.cont5-apply(U.cont5,1,mean)))/(n.reps-1)))
    s2.u47<-mean(diag((U.cont3-apply(U.cont3,1,mean))%*%t((U.cont6-apply(U.cont6,1,mean)))/(n.reps-1)))
    s2.u48<-mean(diag((U.cont3-apply(U.cont3,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))
    
    s2.u56<-mean(diag((U.cont4-apply(U.cont4,1,mean))%*%t((U.cont5-apply(U.cont5,1,mean)))/(n.reps-1)))
    s2.u57<-mean(diag((U.cont4-apply(U.cont4,1,mean))%*%t((U.cont6-apply(U.cont6,1,mean)))/(n.reps-1)))
    s2.u58<-mean(diag((U.cont4-apply(U.cont4,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))
    
    s2.u67<-mean(diag((U.cont5-apply(U.cont5,1,mean))%*%t((U.cont6-apply(U.cont6,1,mean)))/(n.reps-1)))
    s2.u68<-mean(diag((U.cont5-apply(U.cont5,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))
    
    s2.u78<-mean(diag((U.cont6-apply(U.cont6,1,mean))%*%t((U.cont7-apply(U.cont7,1,mean)))/(n.reps-1)))
    
    sd.u <- rbind(c(s2.u11,s2.u12,s2.u13,s2.u14,s2.u15,s2.u16,s2.u17,s2.u18),
                  c(s2.u12,s2.u22,s2.u23,s2.u24,s2.u25,s2.u26,s2.u27,s2.u28),
                  c(s2.u13,s2.u23,s2.u33,s2.u34,s2.u35,s2.u36,s2.u37,s2.u38),
                  c(s2.u14,s2.u24,s2.u34,s2.u44,s2.u45,s2.u46,s2.u47,s2.u48),
                  c(s2.u15,s2.u25,s2.u35,s2.u45,s2.u55,s2.u56,s2.u57,s2.u58),
                  c(s2.u16,s2.u26,s2.u36,s2.u46,s2.u56,s2.u66,s2.u67,s2.u68),
                  c(s2.u17,s2.u27,s2.u37,s2.u47,s2.u57,s2.u67,s2.u77,s2.u78),
                  c(s2.u18,s2.u28,s2.u38,s2.u48,s2.u58,s2.u68,s2.u78,s2.u88))
    
    
  }
  
  return(sd.u)
}


## fit conditional score approach
fit_condscore <- function(ddat,beta.init,sig,eps=1e-4){
  
  ## get data
  dat.err <- ddat$error
  
  ## subject index
  subj<-unique(dat.err$subj.vece)
  
  ## starting vector for beta
  b0<-beta.init 
  s2.u <- sig ## SIGMA
  
  iter<-1; dif<-1
  repeat{ ## stop if convergence criterion reached of after 1000 iterations
    
    exp.correct<-NULL
    cont1.correct<-NULL
    cont2.correct<-NULL
    cont3.correct<-NULL
    cont4.correct<-NULL
    cont5.correct<-NULL
    cont6.correct<-NULL
    cont7.correct<-NULL
    for(j in subj){
      idx <- which(dat.err[,1]==j) ## index for subject
      idx.1 <- which(dat.err[idx,2]==1) ## index for case (for the subject)
      
      ## correct exposure
      exp.vece.corrected <- (dat.err[idx,3]-dat.err[idx[idx.1],3])-(s2.u%*%b0)[1]
      exp.vece.corrected[idx.1] <- 0 ## assign 0 to cases (since the difference between case and case is 0)
      exp.correct<-append(exp.correct,exp.vece.corrected) ## collect results
      
      ## correct covariates
      cont1.vece.corrected <- (dat.err[idx,5]-dat.err[idx[idx.1],5])-(s2.u%*%b0)[2]
      cont1.vece.corrected[idx.1] <- 0
      cont1.correct<-append(cont1.correct,cont1.vece.corrected)
      
      cont2.vece.corrected <- (dat.err[idx,6]-dat.err[idx[idx.1],6])-(s2.u%*%b0)[3]
      cont2.vece.corrected[idx.1] <- 0
      cont2.correct<-append(cont2.correct,cont2.vece.corrected)
      
      cont3.vece.corrected <- (dat.err[idx,7]-dat.err[idx[idx.1],7])-(s2.u%*%b0)[4]
      cont3.vece.corrected[idx.1] <- 0
      cont3.correct<-append(cont3.correct,cont3.vece.corrected)
      
      cont4.vece.corrected <- (dat.err[idx,8]-dat.err[idx[idx.1],8])-(s2.u%*%b0)[5]
      cont4.vece.corrected[idx.1] <- 0
      cont4.correct<-append(cont4.correct,cont4.vece.corrected)
      
      cont5.vece.corrected <- (dat.err[idx,9]-dat.err[idx[idx.1],9])-(s2.u%*%b0)[6]
      cont5.vece.corrected[idx.1] <- 0
      cont5.correct<-append(cont5.correct,cont5.vece.corrected)
      
      cont6.vece.corrected <- (dat.err[idx,10]-dat.err[idx[idx.1],10])-(s2.u%*%b0)[7]
      cont6.vece.corrected[idx.1] <- 0
      cont6.correct<-append(cont6.correct,cont6.vece.corrected)
      
      cont7.vece.corrected <- (dat.err[idx,11]-dat.err[idx[idx.1],11])-(s2.u%*%b0)[8]
      cont7.vece.corrected[idx.1] <- 0
      cont7.correct<-append(cont7.correct,cont7.vece.corrected)
    }
    
    ## corrected data
    data.cc.correct <- cbind(dat.err,exp.correct,cont1.correct,cont2.correct,cont3.correct,cont4.correct,cont5.correct,cont6.correct,cont7.correct)
    
    ## fit model to corrected data
    cc.fit.correct <- tryCatch(clogit(case.vece ~ exp.correct+cont1.correct+cont2.correct+cont3.correct+cont4.correct+cont5.correct+cont6.correct+cont7.correct+strata(subj.vece), data=data.cc.correct),
                               error=function(a){list(coefficients=rep(NA,8))}) ## report NAs if convergence error
    b<-cc.fit.correct$coefficients; dif<-sum(abs(b0-b)) ## get new estimates
    if (sum(is.na(b))>0) break ## stop if error
    if ( dif<eps ) break  ## stop if reached convergence
    if ( iter>1000 ) break ## stop after 1000 iterations 
    
    b0<-b ## get b0 to compute corrected exposures
    iter<-iter+1 ## next iteration
  }
  
  return(b)
  
}



## get lambda
get_lambda <- function(bbeta,yy){
  
  ## get intercept
  bb0 <- glm(yy~offset(X.mat[,2:ncol(X.mat)]%*%bbeta),family=poisson(link='log'))$coefficients[1]
  ## compute lambda based on beta estimate and intercept
  new_lambda <-  exp(X.mat%*%c(bb0,bbeta))
  
  return(new_lambda)
}



### simulate data
sim.data.new<-function(i,lambda,vsize=SS,delays=true_delays){  
  set.seed(1234+4343*i)
  
  y.sim <- rpois(total.hours, lambda)                 # draw outcomes for every hour, (except the first 1000 hrs for some reason)
  sumy <- sum(y.sim)                                  # total outcomes
  
  hours.5years <- hours.seq[1:total.hours]            # only taking relevant 5 years 
  hours.events <- hours.5years[y.sim>0]               # hours at which an event occurs
  hours.events2 <- hours.events                       # duplicate
  y.sim.events <- y.sim[y.sim>0]                      # non-zero events only
  hours.events <- rep(hours.events, y.sim.events)     # vector of hours at which events occur, where each hour is repeated by # of events at that hour
  
  delay <- sample(delays,sumy,replace=T)              # random sample of delay time, one for each event
  hours.events.error <- hours.events+delay            # time of mismeasurement
  
  
  hours.events <- 1001+hours.events                   # +1001 since X.mat actually starts at 1002
  hours.events.error <- 1001+hours.events.error       # +1001 since X.mat actually starts at 1002
  data.cc.true <- MLmatrix[hours.events,4:48]         
  data.cc.error <- MLmatrix[hours.events.error,4:48]
  
  subj.vect <- seq(1,sumy)
  subj.vece <- seq(1,sumy)
  
  case.vect_0 <- rep(1,sumy)
  case.vect_1 <- case.vect_2 <- case.vect_3 <- case.vect_4 <- rep(0,sumy)
  
  case.vece_0 <- rep(1,sumy)
  case.vece_1 <- case.vece_2 <- case.vece_3 <- case.vece_4 <- rep(0,sumy)
  
  data.cc.true <- cbind(subj.vect,data.cc.true,case.vect_0,case.vect_1,case.vect_2,case.vect_3,case.vect_4)
  data.cc.error <- cbind(subj.vece,data.cc.error,case.vece_0,case.vece_1,case.vece_2,case.vece_3,case.vece_4)
  
  data.cc.true <- data.frame(data.cc.true)
  data.cc.error <- data.frame(data.cc.error)
  
  data.cc.true <- reshape(data.cc.true, varying = 2:51, sep = "_", direction = 'long')
  data.cc.error <- reshape(data.cc.error, varying = 2:51, sep = "_", direction = 'long')
  
  data.cc.true <- data.cc.true[order(data.cc.true[,1],data.cc.true[,2],decreasing=FALSE),]
  data.cc.error <- data.cc.error[order(data.cc.error[,1],data.cc.error[,2],decreasing=FALSE),]
  
  data.cc.true <- data.cc.true[,c(1,12,3,11,4,5,6,7,8,9,10)]
  data.cc.error <- data.cc.error[,c(1,12,3,11,4,5,6,7,8,9,10)]
  
  colnames(data.cc.true) <- c("subj.vect", "case.vect", "exp.vect", "time.true", "cont.vect1", "cont.vect2", "cont.vect3", "cont.vect4", "cont.vect5", "cont.vect6", "cont.vect7")
  colnames(data.cc.error) <- c("subj.vece", "case.vece", "exp.vece", "time.error", "cont.vece1", "cont.vece2", "cont.vece3", "cont.vece4", "cont.vece5", "cont.vece6", "cont.vece7")
  
  data.cc.true <- data.cc.true[data.cc.true$exp.vect!=-99 & data.cc.true$cont.vect1!=-99 & data.cc.true$cont.vect2!=-99 & data.cc.true$cont.vect3!=-99 & data.cc.true$cont.vect4!=-99 & data.cc.true$cont.vect5!=-99 & data.cc.true$cont.vect6!=-99 & data.cc.true$cont.vect7!=-99,]
  data.cc.error <- data.cc.error[data.cc.error$exp.vece!=-99 & data.cc.error$cont.vece1!=-99 & data.cc.error$cont.vece2!=-99 & data.cc.error$cont.vece3!=-99 & data.cc.error$cont.vece4!=-99 & data.cc.error$cont.vece5!=-99 & data.cc.error$cont.vece6!=-99 & data.cc.error$cont.vece7!=-99,]
  
  
  ## validation sample
  valid_id <- sample(1:sumy,vsize,replace = FALSE) ## vsize is validation size
  # valid_id <- 1:sumy ## complete validation sample
  
  ## validation distribution of delays
  valid_delays <- delay[valid_id]
  valid_dist <- table(valid_delays)
  valid_del <- as.numeric(names(valid_dist))
  valid_phi <-as.vector(valid_dist/sum(valid_dist)) 
  
  
  return (list(true=data.cc.true, error=data.cc.error, hours.events=hours.events.error, sumy=sumy,
               del=valid_del, phi=valid_phi, delays=valid_delays))     
  
  ## true (data frame): data at true onset time
  ## error (data.frame) : data at delayed time
  ## hours.events (numeric vector) : onset time index
  ## sumy (numeric scalar) : total number of events
  ## del (numeric vector) : validation sample of observed delay times
  ## phi (numeric vector) : relative frequencies of delay times based on validation delays
  
  ## data.cc.xxx consists of :
  ## (1) subject id (numeric), 
  ## (2) case variable (numeric) : '1' for case(event), '0' for control(matched set), 
  ## (3) PM exposure (numeric) : PM level at time index,
  ## (4) time index
  ## and add control variables (temperature or pressure ...)
  
}




##############################
##  Format True/Error Data  ##
##############################
hours.events <- time_true
hours.events.error <- time_admit    ### WITH ADMISSION
sumy <- length(hours.events)
sumy.error <- length(hours.events.error)

data.cc.true <- MLmatrix[hours.events,4:48]           
data.cc.error <- MLmatrix[hours.events.error,4:48]

subj.vect <- seq(1,sumy)
subj.vece <- seq(1,sumy.error)

case.vect_0 <- rep(1,sumy)
case.vect_1 <- case.vect_2 <- case.vect_3 <- case.vect_4 <- rep(0,sumy)

case.vece_0 <- rep(1,sumy.error)
case.vece_1 <- case.vece_2 <- case.vece_3 <- case.vece_4 <- rep(0,sumy.error)

data.cc.true <- cbind(subj.vect,data.cc.true,case.vect_0,case.vect_1,case.vect_2,case.vect_3,case.vect_4)
data.cc.error <- cbind(subj.vece,data.cc.error,case.vece_0,case.vece_1,case.vece_2,case.vece_3,case.vece_4)

data.cc.true <- data.frame(data.cc.true)
data.cc.error <- data.frame(data.cc.error)

data.cc.true <- reshape(data.cc.true, varying = 2:51, sep = "_", direction = 'long')
data.cc.error <- reshape(data.cc.error, varying = 2:51, sep = "_", direction = 'long')

data.cc.true <- data.cc.true[order(data.cc.true[,1],data.cc.true[,2],decreasing=FALSE),]
data.cc.error <- data.cc.error[order(data.cc.error[,1],data.cc.error[,2],decreasing=FALSE),]

data.cc.true <- data.cc.true[,c(1,12,3,11,4,5,6,7,8,9,10)]
data.cc.error <- data.cc.error[,c(1,12,3,11,4,5,6,7,8,9,10)]

colnames(data.cc.true) <- c("subj.vect", "case.vect", "exp.vect", "time.true", "cont.vect1", "cont.vect2", "cont.vect3", "cont.vect4", "cont.vect5", "cont.vect6", "cont.vect7")
colnames(data.cc.error) <- c("subj.vece", "case.vece", "exp.vece", "time.error", "cont.vece1", "cont.vece2", "cont.vece3", "cont.vece4", "cont.vece5", "cont.vece6", "cont.vece7")

data.cc.true <- data.cc.true[data.cc.true$exp.vect!=-99 & data.cc.true$cont.vect1!=-99 & data.cc.true$cont.vect2!=-99 & data.cc.true$cont.vect3!=-99 & data.cc.true$cont.vect4!=-99 & data.cc.true$cont.vect5!=-99 & data.cc.true$cont.vect6!=-99 & data.cc.true$cont.vect7!=-99,]
data.cc.error <- data.cc.error[data.cc.error$exp.vece!=-99 & data.cc.error$cont.vece1!=-99 & data.cc.error$cont.vece2!=-99 & data.cc.error$cont.vece3!=-99 & data.cc.error$cont.vece4!=-99 & data.cc.error$cont.vece5!=-99 & data.cc.error$cont.vece6!=-99 & data.cc.error$cont.vece7!=-99,]


#######################
## Validation Sample ##
#######################
## validation sample of size SS=300
set.seed(1234*1000+99)
valid_id <- sample(1:sumy,SS,replace = FALSE) ## vsize is validation size
## validation distribution of delays
valid_delays <- true_delays[valid_id]
valid_dist <- table(valid_delays)
valid_del <- as.numeric(names(valid_dist))
valid_phi <-as.vector(valid_dist/sum(valid_dist)) 



#####################
##  Data Analysis  ##
#####################

### create dataset
dat <- list(true=data.cc.true,error=data.cc.error,hours.events=hours.events.error,
            del=valid_del,phi=valid_phi,delays=valid_delays)


## true model (no mismeasurement)
cc.fit.true <- clogit(case.vect ~ exp.vect+cont.vect1+cont.vect2+cont.vect3+cont.vect4+cont.vect5+cont.vect6+cont.vect7+strata(subj.vect), data=dat$true)
cc.fit.true.coef <- cc.fit.true$coef
cc.fit.true.SE <- sqrt(diag(cc.fit.true$var))

## naive estimator
cc.fit.error <- clogit(case.vece ~exp.vece+cont.vece1+cont.vece2+cont.vece3+cont.vece4+cont.vece5+cont.vece6+cont.vece7+strata(subj.vece), data=dat$error)
cc.fit.error.coef <- cc.fit.error$coef
cc.fit.error.SE <- sqrt(diag(cc.fit.error$var))

## conditional score method
sigMat <- get_sigma(dat)
cc.fit.cs <- fit_condscore(ddat=dat,beta.init=cc.fit.error.coef,sig=sigMat)
cc.fit.cs.coef <- cc.fit.cs
#cc.fit.cs.SE <- NA ## get SEs from bootstrap

## regression calibration method
cc.fit.rc <- fit_regcal(ddat=dat,beta.init=cc.fit.error.coef)
cc.fit.rc.coef <- cc.fit.rc$est
cc.fit.rc.SE <- cc.fit.rc$SE

## marginalized likelihood method
cc.fit.ml <- fit_MargLik(ddat=dat,beta.init=cc.fit.error.coef)
cc.fit.ml.coef <- cc.fit.ml$est
cc.fit.ml.SE <- cc.fit.ml$SE






##########################
## Parametric Bootstrap ##
##########################
### Bootstrap to obtain bias and standard error for CS, RC and ML

## generate results matrices
boot.error <- boot.cs <- boot.rc <- boot.ml <- matrix(NA,nr=BB,nc=length(cc.fit.true.coef), dimnames=list(NULL, names(cc.fit.true.coef)))

## Get lambdas for parametric bootstrap
b.idx<-which(dat$error$case.vece==1) ## case index
b.case<-dat$error$time.error[b.idx] ## case times
b.Y<-rep(0,total.hours)  ## EDITED: +1000
for (a in 1:length(b.case)) {b.Y[b.case[a]]<-b.Y[b.case[a]]+1}
## compute lambdas
b.lambda.error  <- get_lambda(bbeta=cc.fit.error.coef,yy=b.Y)  ## error (naive estimator)
# b.lambda.cs   <- get_lambda(bbeta=cc.fit.cs.coef,yy=b.Y)     ## conditional score 
b.lambda.rc     <- get_lambda(bbeta=cc.fit.rc.coef,yy=b.Y)     ## regression calibration
b.lambda.ml     <- get_lambda(bbeta=cc.fit.ml.coef,yy=b.Y)     ## marginalized likelihood



## bootstrap loop
for(bb in loopvec){
  k <- 2000+bb ## set seed
  
  ## error estimator (naive)
  b.dat<-sim.data.new(k,b.lambda.error,delays=dat$delays) ## sim data
  boot.fit.error <- clogit(case.vece ~exp.vece+cont.vece1+cont.vece2+cont.vece3+cont.vece4+cont.vece5+cont.vece6+cont.vece7+strata(subj.vece), data=b.dat$error) ## fit data
  if(length(boot.fit.error$coef)==length(cc.fit.true.coef)){ boot.error[bb,] <- boot.fit.error$coef } ## collect ests
  
  # ## conditional score
  # b.dat<-sim.data.new(k,b.lambda.cs,delays=dat$delays) ## sim data
  # sigMat <- get_sigma(b.dat) ## get new sigma
  # boot.fit.cs <- fit_condscore(ddat=b.dat,beta.init=cc.fit.cs.coef,sig=sigMat) ## fit data
  # if(length(boot.fit.cs)==length(cc.fit.true.coef)){ boot.cs[bb,] <- boot.fit.cs } ## collect ests
  
  ## regression calibration
  b.dat<-sim.data.new(k,b.lambda.rc,delays=dat$delays) ## sim data
  boot.fit.rc <- fit_regcal(ddat=b.dat,beta.init=cc.fit.rc.coef) ## fit data
  if(length(boot.fit.rc$est)==length(cc.fit.true.coef)){ boot.rc[bb,] <- boot.fit.rc$est } ## collect ests
  
  ## marginalized likelihood
  b.dat<-sim.data.new(k,b.lambda.ml,delays=dat$delays) ## sim data
  boot.fit.ml <- fit_MargLik(ddat=b.dat,beta.init=cc.fit.ml.coef) ## fit data
  if(length(boot.fit.ml$est)==length(cc.fit.true.coef)){ boot.ml[bb,] <- boot.fit.ml$est } ## collect ests
  
  print(bb)
}

#######################
##  Compile Results  ##
#######################

## apply bias correction and compile results
ests <- cbind(cc.fit.true.coef,cc.fit.error.coef,cc.fit.cs.coef,cc.fit.rc.coef,cc.fit.ml.coef, ## raw versions
              2*cc.fit.error.coef - apply(boot.error,2,mean,na.rm=T), ## bootstrap-corrected versions
              2*cc.fit.cs.coef - apply(boot.cs,2,mean,na.rm=T),  
              2*cc.fit.rc.coef - apply(boot.rc,2,mean,na.rm=T),     
              2*cc.fit.ml.coef - apply(boot.ml,2,mean,na.rm=T))

# sim.est <- ests[,1]       ## first column of ests is for PM2.5
# sim.temp.est <- ests[,2]  ## second column for temperature

SEs <- cbind(cc.fit.true.SE,cc.fit.error.SE,apply(boot.cs,2,sd,na.rm=T),cc.fit.rc.SE,cc.fit.ml.SE, ## raw versions
             apply(boot.error,2,sd,na.rm=T), ## bootstrap-corrected versions
             apply(boot.cs,2,sd,na.rm=T),
             apply(boot.rc,2,sd,na.rm=T),
             apply(boot.ml,2,sd,na.rm=T))

# sim.SE <- SEs[,1]       ## first column of ests is for PM2.5
# sim.temp.SE <- SEs[,2]  ## second column for temperature

## confidence intervals
loCI <- ests-1.96*SEs
hiCI <- ests+1.96*SEs

## combine estimates and intervals
df_res <- c()
for (cc in 1:ncol(ests)){
  df_res <- cbind(df_res,ests[,cc],loCI[,cc],hiCI[,cc])
}
## rename columns
colnames(df_res) <- paste0(rep(c("True","Error","CS","RC","ML","B.Error","B.CS","B.RC","B.ML"),each=3) ,
                           rep(c(".Est",".lowCI",".hiCI"),times=9) )

## save results
write.table(x=df_res, file=paste0(path,'realdata_results_300.txt'), row.names=T, col.names=T, quote=F)
