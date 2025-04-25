# BART on reservoir data

# libraries
# install.packages('dplyr')
library(dplyr)
# install.packages('tidyr')
library(tidyr)
# install.packages("BART")
# force everything to run on one thread: ended up session abort
# Sys.setenv(
#  OMP_NUM_THREADS       = 1,
#  OPENBLAS_NUM_THREADS  = 1,
#  MKL_NUM_THREADS       = 1
# )
# library(BART)
# install.packages("dbarts")
library(dbarts) # worked

# read data
set.seed(12)
res_chars <- read.csv('res_chars.csv')
QI_df <- read.csv('QI.csv')
QO_df <- read.csv('QO.csv')

# Transfer data to merge
QI <- as.matrix(sapply(QI_df[, -1], as.numeric))
QO <- as.matrix(sapply(QO_df[, -1], as.numeric))
rownames(QO) <- QO_df[, 1]
rownames(QI) <- QI_df[, 1]
res_chars$ReservoirID <- as.character(res_chars$X)

# Convert outflow matrix to a long tibble to merge
df_out <- as_tibble(QO,  .name_repair = "minimal") %>% 
  mutate(ReservoirID = rownames(QO)) %>% 
  pivot_longer(
    cols      = -ReservoirID,
    names_to  = "Year",
    values_to = "Outflow"
  )

# Convert inflow matrix to a long tibble to merge
df_in  <- as_tibble(QI,  .name_repair = "minimal") %>% 
  mutate(ReservoirID = rownames(QI)) %>% 
  pivot_longer(
    cols      = -ReservoirID,
    names_to  = "Year",
    values_to = "Inflow"
  )

# Merge them outflow, inflow, and res-char variables
df <- df_out %>%
  left_join(df_in, by = c("ReservoirID","Year")) %>%
  # pick only the 5 cols you need from res_char_df:
  left_join(
    res_chars %>% select(ReservoirID, Capacity, DOR, Depth, CatchArea, AvgRelease),
    by = "ReservoirID"
  ) 
df$Year <- as.numeric(sub("^X", "", df$Year))

# Create X and Y
Y <- df$Outflow
X <- model.matrix(
  ~ Inflow + Capacity + Depth + DOR + CatchArea + AvgRelease,
  data = df)

# train-test-split
test_idx  <- sample(nrow(df), size = 0.2 * nrow(df))
train_idx <- setdiff(seq_len(nrow(df)), test_idx)
X_train   <- X[train_idx, ]
Y_train   <- Y[train_idx]
X_test    <- X[test_idx, ]
Y_test    <- Y[test_idx]

# BART: aborted
fit <- wbart(
  x.train    = X_train, 
  y.train    = Y_train,
  x.test     = X_test,
  ntree      = 5,
  nskip      = 200,    # burn in
  ndpost     = 300,    # posterior draws
  printevery = 0
)

# dbarts
fit2 <- bart(
  x.train = X_train,
  y.train = Y_train,
  x.test  = X_test,
  ntree   = 5,
  verbose = FALSE
)

# Results
yhat_mean <- fit2$yhat.test.mean
yhat_samples <- fit2$yhat.test
yhat_sd <- apply(yhat_samples, 2, sd)

results <- data.frame(
  Observed   = Y_test,       
  Predicted  = yhat_mean,
  Uncertainty = yhat_sd
)


# 1. Compute R²
ss_res <- sum((Y_test - yhat_mean)^2)
ss_tot <- sum((Y_test - mean(Y_test))^2)
R2     <- 1 - ss_res/ss_tot
cat("R-squared:", round(R2, 4), "\n")

# 2. Scatter plot of true vs. predicted
plot(
  Y_test, 
  yhat_mean,
  xlab = "Observed Outflow (Y_test)",
  ylab = "Predicted Outflow (yhat_mean)",
  main = "Observed vs. Predicted Outflow",
  pch  = 16
)
abline(0, 1, lty = 3)  # 45° line for reference











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