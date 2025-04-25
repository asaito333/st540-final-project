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



# Model: ntree=5
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
  main = "Observed vs. Predicted Outflow ntree=5",
  pch  = 16
)
abline(0, 1, lty = 3)  # 45° line for reference


# feature importance
varcount_mat <- fit2$varcount

# 2. Compute the average split‐count per predictor
var_importance <- colMeans(varcount_mat)

# 3. Put into a data.frame and sort descending
var_imp_df <- data.frame(
  Variable   = colnames(X_train),
  Importance = var_importance
)
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]

# 4. Print the sorted importance
print(var_imp_df)

# 5. Barplot of variable importance
barplot(
  var_imp_df$Importance,
  names.arg  = var_imp_df$Variable,
  las        = 2,            # vertical axis labels
  cex.names  = 0.8,          # label size
  main       = "BART Variable Importance",
  ylab       = "Average Split Count"
)





# Model: ntree = 10
fit3 <- bart(
  x.train = X_train,
  y.train = Y_train,
  x.test  = X_test,
  ntree   = 10,
  verbose = FALSE
)


yhat_mean3 <- fit3$yhat.test.mean
yhat_samples3 <- fit3$yhat.test
yhat_sd3 <- apply(yhat_samples3, 2, sd)

results3 <- data.frame(
  Observed   = Y_test,       
  Predicted  = yhat_mean3,
  Uncertainty = yhat_sd3
)


# 1. Compute R²
ss_res3 <- sum((Y_test - yhat_mean3)^2)
ss_tot <- sum((Y_test - mean(Y_test))^2)
R2_3     <- 1 - ss_res3/ss_tot
cat("R-squared:", round(R2_3, 4), "\n")

# 2. Scatter plot of true vs. predicted
plot(
  Y_test, 
  yhat_mean3,
  xlab = "Observed Outflow (Y_test)",
  ylab = "Predicted Outflow (yhat_mean)",
  main = "Observed vs. Predicted Outflow ntree=10",
  pch  = 16
)
abline(0, 1, lty = 3)  # 45° line for reference

# feature importance
varcount_mat <- fit3$varcount

# 2. Compute the average split‐count per predictor
var_importance <- colMeans(varcount_mat)

# 3. Put into a data.frame and sort descending
var_imp_df <- data.frame(
  Variable   = colnames(X_train),
  Importance = var_importance
)
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]

# 4. Print the sorted importance
print(var_imp_df)

# 5. Barplot of variable importance
barplot(
  var_imp_df$Importance,
  names.arg  = var_imp_df$Variable,
  las        = 2,            # vertical axis labels
  cex.names  = 0.8,          # label size
  main       = "BART Variable Importance ntree=10",
  ylab       = "Average Split Count"
)
