# Sample code: https://www4.stat.ncsu.edu/~bjreich/BSM2/Chapter8/Cloud_BART


install.packages("BART")
library(BART)
install.packages("lubridate")
library(lubridate)

filename   <- "https://www4.stat.ncsu.edu/~bjreich/BSM2/Chapter8/cloud_data_clean.csv"
dat        <- read.csv(url(filename))
complete   <- rowSums(is.na(dat))==0
dat        <- dat[complete,]
ID         <- dat$StormID
lead_time  <- dat$lead_time
basin      <- ifelse(dat$basin=="atlantic",1,0)
X          <- as.matrix(dat[,5:22])
Xs         <- scale(X)
HWRF       <- dat$HWRF
NHC        <- dat$NHC
Y          <- dat$VMAX
year       <- year(dat$Date)
lead_times <- sort(unique(lead_time))
nleads     <- length(lead_times)

test        <- year==2017
train       <- !test

pred1    <- pred2 <- rep(NA,length(Y))
varcount <- NULL 
options(warn=-1)
for(t in 1:nleads){
  tr        <- train & lead_time==lead_times[t]
  te        <- test  & lead_time==lead_times[t]
  
  temp      <- capture.output(fit1<-wbart(Xs[tr,],Y[tr],x.test=Xs[te,],
                                          printevery=Inf,ntree=20))
  pred1[te] <- fit1$yhat.test.mean
  varcount  <- cbind(varcount,colMeans(fit1$varcount))
  
  temp      <- capture.output(fit2<-wbart(Xs[tr,],Y[tr],x.test=Xs[te,],
                                          printevery=Inf,ntree=100))
  pred2[te] <- fit2$yhat.test.mean
}
options(warn=0)

mae      <- function(x){mean(abs(x),na.rm=TRUE)}
MAE_NHC  <- aggregate(Y[test]-  NHC[test],list(lead_time[test]),mae)
MAE_HWRF <- aggregate(Y[test]- HWRF[test],list(lead_time[test]),mae)
BART1    <- aggregate(Y[test]-pred1[test],list(lead_time[test]),mae)
BART2    <- aggregate(Y[test]-pred2[test],list(lead_time[test]),mae)