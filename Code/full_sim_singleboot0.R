#################################################
##                                             ##
## Simulation Study, with Parametric Bootstrap ##
##                                             ##
#################################################
## Updated: 11/16/2018
  ## includes PM and AMBIENT TEMPERATURE
  ## SEs and bias correction both from parametric bootstrap
  ## Embarrassingly Parallel

## set TRUE to run locally
runLOCAL <- TRUE
## set FALSE to run job array on cluster via sbatch (requires changing file paths in "Cluster/Local" section below)

####################
##   Parameters   ##
####################
RR <- 2000 ## RR iterations
BB <- 200 ## BB parametric bootstrap resamples
JJ <- 200 ## number of jobs in array ## only used if using job array via sbatch (runLocal=FALSE below)

## exposure effect (make sure to change b0 just below as we change b1)
## uncomment line for different exposure effects used in paper
beta1<- 0.00; beta0<- -3.811837; version_suffix <- 0  ##-3.8
# beta1<-0.015; beta0<- -3.955918; version_suffix <- 15  ##-3.95
# beta1<-0.03; beta0<- -4.1; version_suffix <- 3
# beta1<-0.045; beta0<- -4.244081; version_suffix <- 45  ##-4.24

## all coefficients
beta <- rbind(beta0,beta1,0.017) 


#######################
##   Load packages   ##
#######################
library(survival)
library(foreign)
# library(foreach)
# library(doMPI)



#######################
##   Cluster/Local   ##
#######################
## set up parallel
if(runLOCAL==TRUE){
  path <- "../Results/" ## path for results
  datpath <- "../Data/" ## path for data
  loopvec <- 1:RR ## run entire loop
  suffix <- "_local" ## output file has a suffix to distinguish from other results files
}else{ ## paths need to be changed to correspond with cluster
  path <- "~/Results/" ## path for results
  datpath <- "~/Data/" ## path for data
  args <- commandArgs(trailingOnly = TRUE) ## collect arguments from shell script
  iter_no <- as.integer(args[1])  ## here there is only one argument, the iteration number
  loopvec <- ((RR/JJ)*(iter_no-1) +1):((RR/JJ)*iter_no) ## portion of the loop to run
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
true_delays <- read.table(paste0(datpath,"delays.csv"), header=TRUE, quote="\"")
covs <- read.csv(paste0(datpath,"brentpmhourly.csv"))
MLmatrix <- read.table(paste0(datpath,"MLmatrix.txt"), header=TRUE, quote="\"")


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


##########################
##  Create data matrix  ##
##########################
## time scale
hours.day <- 24
days.year <- 365
years <- 5

days.month <- c(31,28,31,30,31,30,31,31,30,31,30,31)
days.in.month <- c(seq(1,31,1),seq(1,28,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,31,1),seq(1,30,1),seq(1,31,1),seq(1,30,1),seq(1,31,1))

total.hours <- hours.day*days.year*years
total.days <- days.year*years

month.seq <- c(rep(1,days.month[1]),rep(2,days.month[2]),rep(3,days.month[3]),rep(4,days.month[4]),rep(5,days.month[5]),rep(6,days.month[6]),
               rep(7,days.month[7]),rep(8,days.month[8]),rep(9,days.month[9]),rep(10,days.month[10]),rep(11,days.month[11]),rep(12,days.month[12]))
month.seq <- rep(month.seq, rep(24,365))
month.seq <- c(month.seq, month.seq+12, month.seq+24, month.seq+36, month.seq+48, month.seq+60)

hours.in.day <- rep(seq(1,24,1),total.days+365)
hours.seq <- seq(1,total.hours+(365*24))
day.seq <- (floor((hours.seq-1)/24))+1
week.seq <- floor((day.seq-1) / 7)+1
day.in.week <- rep(seq(1,7,1),each=24,length.out=(total.days+365)*24)


## delay times
true_delays <- true_delays[true_delays>0]
tab<-table(true_delays)    
frq<-tab; del<-as.numeric(names(tab)) # del is only the unique values because well multiply by...
phi<-as.vector(frq/sum(frq))  # delay distribution


## covariate time series
intercept <- rep(1,total.hours+1000)
X.mat <- cbind(intercept, pm.seq[1002:(1001+total.hours+1000)]) 
control1.full <- temp[1:(1001+total.hours+1000)] ## select temperature (dry bulb) as confounder/covariate/control variable
control1.full <- control1.full-mean(control1.full[1002:(1001+total.hours)]) ## center data 
X.mat <- cbind(X.mat,control1.full[1002:(1001+total.hours+1000)])
lambda <-  exp(X.mat%*%beta) ## calculate mean for counts


## import MLmatrix 
## includes PM and temperature
## faster just to load in than to wait for it to populate
## code to create MLmatrix commented out below
MLmatrix <- data.matrix(MLmatrix)


# ##############
# ### MATRIX ###
# ##############
# 
# # total.hours indexes row by hour (ie same as the row index)
# hours.5years <- c(-1000:0,hours.seq[1:(total.hours+1000)])   ## added 1001 hours to beginning and 1000 to end just to be sure
# 
# # n.matches dictates how many (4 or 5) observations are in a matched set for a specific event time (including the actual event/case)
# n.days.in.month <- c(rep(31,31),rep(28,28),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31),rep(31,31),rep(30,30),rep(31,31),rep(30,30),rep(31,31))
# which.day <- (days.in.month-1)%%7+1     # indexes days 1:7 only within a month (unlike days.in.week)
# n.matches <- 4+(which.day+7*4<=n.days.in.month)
# n.matches <- rep(rep(n.matches,each=24),5)
# n.matches <- c(tail(n.matches,1001),n.matches,head(n.matches,1000))    # add 1001 to beginning and 1000 to end
# 
# # which week of the month is it
# which.week <- 1+(days.in.month>7)+(days.in.month>14)+(days.in.month>21)+(days.in.month>28)
# which.week <- rep(rep(which.week,each=24),5)
# which.week <- c(tail(which.week,1001),which.week,head(which.week,1000))    # add 1001 to beginning and 1000 to end
# 
# # exposure for case
# exp_0 <- pm.seq[1:(1001+total.hours+1000)]      ## added 1001 hours to beginning and 1000 to end just to be sure
# # control variables/confounders for case
# control1_0 <- control1.full    # full because we have added on the 1001 rows at the beginning
# controlz <- data.matrix(control1.full)
# 
# # hours (to be used for sim data)
# hours_0 <- hours.5years
# 
# 
# # initialize 
# exp_1 <- exp_2 <- exp_3 <- exp_4 <- rep(-99,1001+total.hours+1000)                        # use -99 to avoid loops later on?
# control1_1 <- control1_2 <- control1_3 <- control1_4 <- rep(-99,1001+total.hours+1000)
# hours_1 <- hours_2 <- hours_3 <- hours_4 <- rep(-99,1001+total.hours+1000)
# 
# # initialize matrix
# MLmatrix <- data.frame(cbind(hours.5years,n.matches,which.week,
#                              exp_0,control1_0,
#                              exp_1,control1_1,
#                              exp_2,control1_2,
#                              exp_3,control1_3,
#                              exp_4,control1_4,
#                              hours_0,hours_1,hours_2,hours_3,hours_4))
# 
# # restrict everything to the 5 years
# month.seq.5years <- month.seq[1:43800]   
# month.seq.5years <- c(tail(month.seq.5years,1001)-60,month.seq.5years,head(month.seq.5years,1000)+60)   # add 1001 to beginning and 1000 to end
# day.in.week.5years <- day.in.week[1:43800]
# day.in.week.5years <- c(tail(rep(seq(1,7),150),1001),day.in.week.5years,head(rep(seq(1,7),150),1000))       # add 1001 to beginning and 1000 to end
# hours.in.day.5years <- hours.in.day[1:43800]
# hours.in.day.5years <- c(tail(rep(seq(1,24),50),1001),hours.in.day.5years,head(rep(seq(1,24),50),1000))       # add 1001 to beginning and 1000 to end
# 
# for (jj in 1:(1001+total.hours+1000)){
#   # subset matched data; exclude case since its already in there
#   matches <- MLmatrix[month.seq.5years==month.seq.5years[jj] & day.in.week.5years==day.in.week.5years[jj] & hours.in.day.5years==hours.in.day.5years[jj] & hours.5years!=jj-1001,]
#   #first match
#   MLmatrix[jj,5+(1:2)] <- matches[1,4:5]
#   MLmatrix[jj,15] <- matches[1,14]      # to have a time column for sim.data
#   # second match
#   MLmatrix[jj,7+(1:2)] <- matches[2,4:5]
#   MLmatrix[jj,16] <- matches[2,14]
#   # third match
#   MLmatrix[jj,9+(1:2)] <- matches[3,4:5]
#   MLmatrix[jj,17] <- matches[3,14]
#   # possible 4th match
#   if ( nrow(matches)==4){
#     MLmatrix[jj,11+(1:2)] <- matches[4,4:5]
#     MLmatrix[jj,18] <- matches[4,14]
#     
#   }  
# }
# 
# match5 <- as.numeric(MLmatrix$n.matches==5)     ### indicator of whether or not to include last match in denominator calculation
# MLmatrix <- cbind(MLmatrix,match5)              ### add indicator to make code run faster
# # write.table(x=MLmatrix, file='~/MeasErr/Results/MLmatrix_temp.txt', row.names=F, col.names=T, quote=F)
# MLmatrix <- data.matrix(MLmatrix)





############################
##  Function Definitions  ##
############################

### negative logLikelihood (to be maximized)
neg.logL.MLnew<-function(b,dataframe,del,phi){
  subj.times <- 1001+dataframe$time.error[dataframe$case.vece==1]     # shift by 1001
  res<-0
  #ct.nan<-0 ## for troubleshooting NaNs
  for(a in subj.times){
    
    subj.data <- MLmatrix[a-del,]
    
    case.contrib <- exp(subj.data[,2+2*1+(0:1) ]%*%b )                            # numerator
    match.contrib <- case.contrib+
      exp(subj.data[,2+2*2+(0:1) ]%*%b )+
      exp(subj.data[,2+2*3+(0:1) ]%*%b )+
      exp(subj.data[,2+2*4+(0:1) ]%*%b )+
      exp(subj.data[,2+2*5+(0:1) ]%*%b )*subj.data[,ncol(subj.data)]     # multiply by indicator 
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
  subj.times <- 1001+dataframe$time.error[dataframe$case.vece==1]     # shift by 1001
  res <- rep(0,2)
  #ct.nan<-0
  for(a in subj.times){
    subj.data <- MLmatrix[a-del,]
    
    tmp.1.1 <- c(exp(subj.data[,2+2*1+(0:1) ]%*%b )) 
    tmp.1.2 <- subj.data[,4]*tmp.1.1
    tmp.g.1 <- subj.data[,5]*tmp.1.1  
    # initialize (including first term from above)
    tmp.2.1 <- tmp.1.1
    tmp.2.2 <- tmp.1.2
    tmp.g.2 <- tmp.g.1                                                
    for (nn in 2:4){                               # loop over matches (depending on number of matches)
      tmp.nn.1 <- c(exp(subj.data[,2+2*nn+(0:1) ]%*%b ) )
      tmp.2.1 <- tmp.2.1+tmp.nn.1
      tmp.2.2 <- tmp.2.2+subj.data[,2+2*nn]*tmp.nn.1
      tmp.g.2 <- tmp.g.2+subj.data[,2+2*nn+(1)]*tmp.nn.1
    }
    
    # potential match
    tmp.5.1 <- c(exp(subj.data[,2+2*5+(0:1) ]%*%b )*subj.data[,ncol(subj.data)] )     # multiply by indicator 
    tmp.2.1 <- tmp.2.1+tmp.5.1
    tmp.2.2 <- tmp.2.2+subj.data[,2+2*5]*tmp.5.1
    tmp.g.2 <- tmp.g.2+subj.data[,2+2*5+(1)]*tmp.5.1
    
    A.b <- sum((tmp.1.2/tmp.2.1)*phi)
    B.b <- sum((tmp.1.1/tmp.2.1*tmp.2.2/tmp.2.1)*phi)
    C <- sum((tmp.1.1/tmp.2.1)*phi)
    
    A.g<-sum(tmp.g.1/tmp.2.1 * phi)
    B.g<-sum(tmp.g.2/tmp.2.1*tmp.1.1/tmp.2.1 * phi)
    res1<-(A.b-B.b)/C
    res2<-(A.g-B.g)/C
    
    res <- res+c(res1,res2)
    #if ( res.1=='NaN' | is.na(res.1)==TRUE ) ct.nan<-ct.nan+1
    #if ( res.1!='NaN' & is.na(res.1)!=TRUE ) res<-res+res.1
  }
  return ( -res ) 
}

### Information Matrix
infomat.MLnew<-function(b,dataframe,del,phi){
  # infomat.ml : information which considers 3 additional control variables, resulting in a matrix of size 4x4
  subj.times <- 1001+dataframe$time.error[dataframe$case.vece==1]     # shift by 1001
  I.11<-0; I.21<-0; I.22<-0; ct.nan<-0
  for(a in subj.times){
    subj.data <- MLmatrix[a-del,]
    
    tmp.1.1 <- c(exp(subj.data[,2+2*1+(0:1) ]%*%b )) 
    tmp.1.2 <- subj.data[,4]*tmp.1.1
    tmp.1.3 <- ((subj.data[,4])^2)*tmp.1.1
    tmp.g.1 <- subj.data[,5]*tmp.1.1  
    tmp.gg.1 <- (subj.data[,5]^2)*tmp.1.1   
    
    # initialize (including first term from above)
    tmp.2.1 <- tmp.1.1
    tmp.2.2 <- tmp.1.2
    tmp.2.3 <- tmp.1.3
    tmp.bg.1 <- subj.data[,5]*tmp.1.1
    tmp.bg.2 <- subj.data[,5]*subj.data[,4]*tmp.1.1
    tmp.gg.2 <- (3*(subj.data[,5]^2) )*tmp.1.1
    
    ## loop over matches (depending on number of matches)
    for (nn in 2:4){                               
      tmp.nn.1 <- c(exp(subj.data[,2+2*nn+(0:1) ]%*%b ) )
      tmp.2.1 <- tmp.2.1+tmp.nn.1
      tmp.2.2 <- tmp.2.2+subj.data[,2+2*nn]*tmp.nn.1
      tmp.2.3 <- tmp.2.3+((subj.data[,2+2*nn])^2)*tmp.nn.1
      tmp.bg.1 <- tmp.bg.1+subj.data[,2+2*nn+(1)]*tmp.nn.1
      tmp.bg.2 <- tmp.bg.2+subj.data[,2+2*nn+(1)]*(subj.data[,2+2*nn]*tmp.nn.1)
      tmp.gg.2 <- tmp.gg.2+( 2*subj.data[,2+2*nn+(1)]*subj.data[,2+2*1+(1)]+subj.data[,2+2*nn+(1)]^2 )*tmp.nn.1  
    }
    
    ## potential 5th match
    tmp.5.1 <- c(exp(subj.data[,2+2*5+(0:1) ]%*%b )*subj.data[,ncol(subj.data)])     # multiply by indicator 
    tmp.2.1 <- tmp.2.1+tmp.5.1
    tmp.2.2 <- tmp.2.2+subj.data[,2+2*5]*tmp.5.1
    tmp.2.3 <- tmp.2.3+((subj.data[,2+2*5])^2)*tmp.5.1 
    tmp.bg.1 <- tmp.bg.1+subj.data[,2+2*5+(1)]*tmp.5.1
    tmp.bg.2 <- tmp.bg.2+subj.data[,2+2*5+(1)]*(subj.data[,2+2*5]*tmp.5.1)
    tmp.gg.2 <- tmp.gg.2+( 2*subj.data[,2+2*5+(1)]*subj.data[,2+2*1+(1)]+subj.data[,2+2*5+(1)]^2 )*tmp.5.1  
    
    ## putting it all together
    A.b <- sum((tmp.1.2/tmp.2.1)*phi)
    B.b <- sum((tmp.1.1/tmp.2.1*tmp.2.2/tmp.2.1)*phi)
    A.g <- sum(tmp.g.1/tmp.2.1 * phi)
    B.g <- sum(tmp.bg.1/tmp.2.1*tmp.1.1/tmp.2.1 * phi)
    C <- sum((tmp.1.1/tmp.2.1)*phi)
    A.bb <- sum(( tmp.1.3/tmp.2.1 - tmp.1.2/tmp.2.1*tmp.2.2/tmp.2.1 ) * phi)
    B.bb <- sum(( tmp.1.2/tmp.2.1*tmp.2.2/tmp.2.1 + tmp.1.1/tmp.2.1*tmp.2.3/tmp.2.1 - 2*tmp.1.1/tmp.2.1*(tmp.2.2/tmp.2.1)^2 ) * phi)
    i.11 <- (A.bb-B.bb)/C - (A.b-B.b)^2/C^2
    
    A.bg <- sum((tmp.1.2/tmp.2.1) * subj.data[,5] * phi)- sum( (tmp.1.1/tmp.2.1*tmp.2.2/tmp.2.1) * subj.data[,5] * phi)
    B.bg <- sum((tmp.1.2/tmp.2.1*tmp.bg.1/tmp.2.1) * phi) + sum(tmp.1.1/tmp.2.1*tmp.bg.2/tmp.2.1 * phi) - 2* sum(tmp.1.1/tmp.2.1*tmp.2.2/tmp.2.1*tmp.bg.1/tmp.2.1 * phi)
    i.21 <- (A.bg-B.bg)/C - (A.g-B.g)*(A.b-B.b)/C^2
    
    
    ## all possible cross terms (theres probably a faster way to do it but its not that important)
    tmp.gg.3 <- 0
    for (nn in 1:5){                               
      
      if (nn==5){
        tmp.nn.1 <- c(exp(subj.data[,2+2*nn+(0:1) ]%*%b )*subj.data[,ncol(subj.data)])     # multiply by indicator 
      } else {
        tmp.nn.1 <- c(exp(subj.data[,2+2*nn+(0:1) ]%*%b ) )
      }
      for (mm in 1:5){    
        if (mm==5){
          tmp.mm.1 <- c(exp(subj.data[,2+2*mm+(0:1) ]%*%b )*subj.data[,ncol(subj.data)])     # multiply by indicator 
        } else {
          tmp.mm.1 <- c(exp(subj.data[,2+2*mm+(0:1) ]%*%b ) )
        }
        
        tmp.gg.3 <- tmp.gg.3+( subj.data[,2+2*nn+(1)]*subj.data[,2+2*mm+(1)] )*tmp.nn.1*tmp.mm.1
      }      
    }
    
    chunk1 <- sum((tmp.gg.1/tmp.2.1) * phi)
    chunk2 <- sum((tmp.1.1/tmp.2.1*tmp.gg.2/tmp.2.1) * phi)
    chunk3 <- sum((tmp.1.1/tmp.2.1*( tmp.gg.3   )/(tmp.2.1)^2) * phi)
    A.gg_B.gg <- chunk1 - chunk2 + 2*chunk3    
    i.22 <- A.gg_B.gg/C - (A.g-B.g)%*%t(A.g-B.g)/C^2
    
    
    if ( i.11=='NaN' | sum(i.21=='NaN')>0 | sum(i.22=='NaN')>0 | is.na(i.11)==TRUE | sum(is.na(i.21)==TRUE)>0 | sum(is.na(i.22)==TRUE)>0 ) ct.nan<-ct.nan+1
    if ( i.11!='NaN' & sum(i.21=='NaN')==0 & sum(i.22=='NaN')==0 & is.na(i.11)!=TRUE & sum(is.na(i.21)==TRUE)==0 & sum(is.na(i.22)==TRUE)==0 ) { I.11<-I.11+i.11; I.21<-I.21+i.21; I.22<-I.22+i.22 }
    
  }
  res <- cbind( c(I.11,I.21), rbind(t(I.21),I.22) )
  return( -res )
}

### Marginal likelihood approach  
fit_MargLik <- function(ddat,beta.init){
  
  ## beta.init<-cc.fit.error$coefficients ## starting values
  
  fit<-tryCatch(optim(beta.init,neg.logL.MLnew,grad.MLnew,dataframe=ddat$error,del=ddat$del,phi=ddat$phi,method='L-BFGS-B',hessian=TRUE),error=function(a){NA}) # if hessian=TRUE, then hessian can be computed numerically
  if(sum(is.na(fit))>0){ ## in case of convergence error
    return(list(est=c(NA,NA),SE=c(NA,NA)))
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
  ref.time.error<-unique(data.cc.error$time.error)          # unique error times
  
  pred.pm<-rep(0,length(ref.time.error))
  pred.cont<-rep(0,length(ref.time.error))
  for(a in 1:length(ref.time.error)){                       # loop over unique error times
    idx<-ref.time.error[a]-ddat$del                         # index for possible exposure times (to be summed over)
    tmp.x<-pm.seq[1001+idx]                                 # each possible exposure to sum over
    tmp.cont<-control1.full[1001+idx] 
    pred.pm[a]<-sum(ddat$phi*tmp.x)                         # expected exposure 
    pred.cont[a]<-ddat$phi%*%tmp.cont
    
  }
  exp.vecr<-rep(0,nrow(data.cc.error))          
  cont.vecr<-rep(0,nrow(data.cc.error))     
  for(a in 1:nrow(data.cc.error)){
    exp.vecr[a]<-pred.pm[which(ref.time.error==data.cc.error$time.error[a])]    # assign predicted vals (do it this way because we found predicteds for unique times--> must distribute)
    cont.vecr[a]<-pred.cont[which(ref.time.error==data.cc.error$time.error[a])]
    
  }
  data.cc.RC<-data.cc.error 
  data.cc.RC[,3]<-exp.vecr          # replace exp.vece with expected value exp.vecr -- it will still be called exp.vece in the data
  data.cc.RC[,5]<-cont.vecr
  fit <- clogit(case.vece ~ exp.vece+cont.vece1+strata(subj.vece), data=data.cc.RC)
  est <- fit$coef
  SE <- sqrt(diag(fit$var))
  
  return(list(est=est,SE=SE))
  
}


## get SIGMA.u for conditional score approach
get_sigma <- function(ddat,n.reps=50){
  sd.u <- NA
  while(sum(is.na(sd.u))>0){
    # estimate measurement error variance
    U<-NULL
    U.cont1<-NULL
    for(a in 1:n.reps){  
      
      ## sampling from validation sample of delay times (300 of them)
      delay2 <- sample(ddat$delays,length(ddat$hours.events),replace=T) 
      
      ## compute U 
      u<-pm.seq[1001+ddat$hours.events+delay2] - pm.seq[1001+ddat$hours.events]  
      u.cont1<-control1.full[1001+ddat$hours.events+delay2]-control1.full[1001+ddat$hours.events]
      
      ## collect data
      U<-cbind(U,u)
      U.cont1<-cbind(U.cont1,u.cont1)
    }
    
    s2.u11<-mean(apply(U,1,var))
    s2.u22<-mean(apply(U.cont1,1,var))
    s2.u12<-mean(diag((U-apply(U,1,mean))%*%t((U.cont1-apply(U.cont1,1,mean)))/(n.reps-1)))
    sd.u <- cbind(c(s2.u11,s2.u12),c(s2.u12,s2.u22))
    
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
    cont1.correct <- NULL
    for(j in subj){
      idx <- which(dat.err[,1]==j) ## index for subject
      idx.1 <- which(dat.err[idx,2]==1) ## index for case (for the subject)
      
      ## correct exposure
      exp.vece.corrected <- (dat.err[idx,3]-dat.err[idx[idx.1],3])-(s2.u%*%b0)[1]
      exp.vece.corrected[idx.1] <- 0 ## assign 0 to cases (since the difference between case and case is 0)
      exp.correct<-append(exp.correct,exp.vece.corrected) ## collect results
      
      ## correct covariate
      cont1.vece.corrected <- (dat.err[idx,5]-dat.err[idx[idx.1],5])-(s2.u%*%b0)[2]
      cont1.vece.corrected[idx.1] <- 0 ## assign 0 to cases (since the difference between case and case is 0)
      cont1.correct<-append(cont1.correct,cont1.vece.corrected) ## collect results
    }
    
    data.cc.correct <- cbind(dat.err,exp.correct,cont1.correct) ## corrected data
    
    ## fit model to corrected data
    cc.fit.correct <- tryCatch(clogit(case.vece ~ exp.correct+cont1.correct+strata(subj.vece), data=data.cc.correct),
                               error=function(a){list(coefficients=rep(NA,2))}) ## report NAs if convergence error
    b<-cc.fit.correct$coefficients; dif<-sum(abs(b0-b)) ## get new estimates
    if (sum(is.na(b))>0) break ## stop if error
    if ( dif<eps ) break  ## stop if reached convergence
    if ( iter>1000 ) break ## stop after 1000 iterations 
    
    b0<-b ## get b0 to compute corrected exposures
    iter<-iter+1 ## next iteration
  }
  
  return(b)
  
}


get_lambda <- function(bbeta,yy){
  
  ## get intercept
  bb0 <- glm(yy~offset(X.mat[,2:ncol(X.mat)]%*%bbeta),family=poisson(link='log'))$coefficients[1]
  ## compute lambda based on beta estimate and intercept
  new_lambda <-  exp(X.mat%*%c(bb0,bbeta))
  
  return(new_lambda)
}

### simulate data
sim.data.new<-function(i,lambda,vsize=300,delays=true_delays){  
  set.seed(1234+4343*i)
  
  y.sim <- rpois(total.hours, lambda)                 # draw outcomes for every hour, (except the first 1000 hrs for some reason)
  sumy <- sum(y.sim)                                  # total outcomes
  
  hours.5years <- hours.seq[1:total.hours]            # only taking relevant 5 years (I think)
  hours.events <- hours.5years[y.sim>0]               # hours at which an event occurs
  hours.events2 <- hours.events                       # duplicate
  y.sim.events <- y.sim[y.sim>0]                      # non-zero events only
  hours.events <- rep(hours.events, y.sim.events)     # vector of hours at which events occur, where each hour is repeated by # of events at that hour
  
  delay <- sample(delays,sumy,replace=T)              # random sample of delay time, one for each event
  hours.events.error <- hours.events+delay            # time of mismeasurement
  
  
  data.cc.true <- MLmatrix[1001+hours.events,4:18]           
  data.cc.error <- MLmatrix[1001+hours.events.error,4:18]
  
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
  
  data.cc.true <- reshape(data.cc.true, varying = 2:21, sep = "_", direction = 'long')
  data.cc.error <- reshape(data.cc.error, varying = 2:21, sep = "_", direction = 'long')
  
  data.cc.true <- data.cc.true[order(data.cc.true[,1],data.cc.true[,2],decreasing=FALSE),]
  data.cc.error <- data.cc.error[order(data.cc.error[,1],data.cc.error[,2],decreasing=FALSE),]
  
  data.cc.true <- data.cc.true[,c(1,6,3,5,4)]
  data.cc.error <- data.cc.error[,c(1,6,3,5,4)]
  
  colnames(data.cc.true) <- c("subj.vect", "case.vect", "exp.vect", "time.true", "cont.vect1")
  colnames(data.cc.error) <- c("subj.vece", "case.vece", "exp.vece", "time.error", "cont.vece1")
  
  data.cc.true <- data.cc.true[data.cc.true$exp.vect!=-99 & data.cc.true$cont.vect1!=-99,]
  data.cc.error <- data.cc.error[data.cc.error$exp.vece!=-99 & data.cc.error$cont.vece1!=-99,]
  
  ## validation sample
  valid_id <- sample(1:sumy,vsize,replace = FALSE) ## vsize is validation size
  # valid_id <- 1:sumy ## full validation sample
  
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


############
##  Loop  ##
############
results <- c() ## initialize results

for(rr in loopvec){
  
  #####################
  ## Generate Events ##
  #####################
  
  ## generate true dataset
  i<-1000+rr
  dat<-sim.data.new(i,lambda,delays=true_delays)
  
  
  ###################
  ## Main Analyses ##
  ###################
  
  ## true model (no mismeasurement)
  cc.fit.true <- clogit(case.vect ~ exp.vect+cont.vect1+strata(subj.vect), data=dat$true)
  cc.fit.true.coef <- cc.fit.true$coef
  cc.fit.true.SE <- sqrt(diag(cc.fit.true$var))
  
  ## naive method (not adjusting for errors)
  cc.fit.error <- clogit(case.vece ~ exp.vece+cont.vece1+strata(subj.vece), data=dat$error)
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
  boot.error <- boot.cs <- boot.rc <- boot.ml <- matrix(NA,nr=BB,nc=length(beta)-1, dimnames=list(NULL, c('exposure','control1')))
  
  ## Get lambdas for parametric bootstrap
  b.idx<-which(dat$error$case.vece==1) ## case index
  b.case<-dat$error$time.error[b.idx] ## case times
  b.Y<-rep(0,total.hours+1000)  ## EDITED: +1000
  for (a in 1:length(b.case)) {b.Y[b.case[a]]<-b.Y[b.case[a]]+1}
  ## compute lambdas
  b.lambda.error  <- get_lambda(bbeta=cc.fit.error.coef,yy=b.Y)  ## error (naive estimator)
  if(sum(is.na(cc.fit.cs.coef))==0){ ## if statement in case of NAs
    b.lambda.cs   <- get_lambda(bbeta=cc.fit.cs.coef,yy=b.Y)     ## conditional score 
  }  else{ b.lambda.cs <- b.lambda.error }   
  b.lambda.rc     <- get_lambda(bbeta=cc.fit.rc.coef,yy=b.Y)     ## regression calibration
  b.lambda.ml     <- get_lambda(bbeta=cc.fit.ml.coef,yy=b.Y)     ## marginalized likelihood
  
  ## bootstrap loop
  for(bb in 1:BB){
    k <- 2000+bb ## set seed
    
    ## error estimator (naive)
    b.dat<-sim.data.new(k,b.lambda.error,delays=dat$delays) ## sim data
    boot.fit.error <- clogit(case.vece ~ exp.vece+cont.vece1+strata(subj.vece), data=b.dat$error) ## fit data
    if(length(boot.fit.error$coef)==2){ boot.error[bb,] <- boot.fit.error$coef } ## collect ests
    
    ## conditional score
    b.dat<-sim.data.new(k,b.lambda.cs,delays=dat$delays) ## sim data
    sigMat <- get_sigma(b.dat) ## get new sigma
    boot.fit.cs <- fit_condscore(ddat=b.dat,beta.init=cc.fit.cs.coef,sig=sigMat) ## fit data
    if(length(boot.fit.cs)==2){ boot.cs[bb,] <- boot.fit.cs } ## collect ests
    
    ## regression calibration
    b.dat<-sim.data.new(k,b.lambda.rc,delays=dat$delays) ## sim data
    boot.fit.rc <- fit_regcal(ddat=b.dat,beta.init=cc.fit.rc.coef) ## fit data
    if(length(boot.fit.rc$est)==2){ boot.rc[bb,] <- boot.fit.rc$est } ## collect ests
    
    ## marginalized likelihood
    b.dat<-sim.data.new(k,b.lambda.ml,delays=dat$delays) ## sim data
    boot.fit.ml <- fit_MargLik(ddat=b.dat,beta.init=cc.fit.ml.coef) ## fit data
    if(length(boot.fit.ml$est)==2){ boot.ml[bb,] <- boot.fit.ml$est } ## collect ests
    
  }

  ## apply bias correction and compile results
  ests <- rbind(cc.fit.true.coef,cc.fit.error.coef,cc.fit.cs.coef,cc.fit.rc.coef,cc.fit.ml.coef, ## raw versions
               2*cc.fit.error.coef - apply(boot.error,2,mean,na.rm=T), ## bootstrap-corrected versions
               2*cc.fit.cs.coef - apply(boot.cs,2,mean,na.rm=T),  
               2*cc.fit.rc.coef - apply(boot.rc,2,mean,na.rm=T),     
               2*cc.fit.ml.coef - apply(boot.ml,2,mean,na.rm=T))
  
  sim.est <- ests[,1]       ## first column of ests is for PM2.5
  sim.temp.est <- ests[,2]  ## second column for temperature
  
  SEs <- rbind(cc.fit.true.SE,cc.fit.error.SE,apply(boot.cs,2,sd,na.rm=T),cc.fit.rc.SE,cc.fit.ml.SE, ## raw versions
               apply(boot.error,2,sd,na.rm=T), ## bootstrap-corrected versions
               apply(boot.cs,2,sd,na.rm=T),
               apply(boot.rc,2,sd,na.rm=T),
               apply(boot.ml,2,sd,na.rm=T))
  
  sim.SE <- SEs[,1]       ## first column of ests is for PM2.5
  sim.temp.SE <- SEs[,2]  ## second column for temperature
  
  
  ## collect results
  results <- rbind(results,c(sim.est,sim.SE,sim.temp.est,sim.temp.SE))
  
  print(rr)
}



## save results
colnames(results) <- paste(rep(c("est","SE","est.temp","SE.temp"),each=9),
                           rep(c("True","Error","CS","RC","ML","B.Error","B.CS","B.RC","B.ML"),times=4))
write.table(x=results, file=paste0(path,'sim_par',version_suffix,'_results',suffix,'.txt'))

  
  
  
