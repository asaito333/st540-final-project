# BART
# temporal and spatial train/test split

library(dplyr)
library(tidyr)
library(dbarts) 

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

# 1) Temporal split: Years 0–15 for training, 17–20 for testing
temporal_train <- df$Year >= 0  & df$Year <= 15
temporal_test  <- df$Year >= 16 & df$Year <= 19

# check sizes
table(temporal_train)
table(temporal_test)

# Split X and Y
X_train_temp <- X[temporal_train, ]
Y_train_temp <- Y[temporal_train]
X_test_temp    <- X[temporal_test, ]
Y_test_temp    <- Y[temporal_test]

# 2) Spatial split: sample 80% of ReservoirID for training
all_ids    <- unique(df$ReservoirID)
n_train_id <- floor(0.8 * length(all_ids))

train_ids  <- sample(all_ids, n_train_id)
test_ids   <- setdiff(all_ids, train_ids)

# or as logical vectors:
spatial_train_ix <- df$ReservoirID %in% train_ids
spatial_test_ix  <- df$ReservoirID %in% test_ids

table(spatial_train_ix)
table(spatial_test_ix)

# Split X and Y
X_train_sp <- X[spatial_train_ix, ]
Y_train_sp <- Y[spatial_train_ix]
X_test_sp    <- X[spatial_test_ix, ]
Y_test_sp    <- Y[spatial_test_ix]


# Model: temporal, ntree=5
fit_temp <- bart(
  x.train = X_train_temp,
  y.train = Y_train_temp,
  x.test  = X_test_temp,
  ntree   = 5,
  verbose = FALSE
)

# Results
yhat_mean_temp <- fit_temp$yhat.test.mean
yhat_samples_temp <- fit_temp$yhat.test
yhat_sd_temp <- apply(yhat_samples_temp, 2, sd)

results_temp <- data.frame(
  Observed   = Y_test_temp,       
  Predicted  = yhat_mean_temp,
  Uncertainty = yhat_sd_temp
)

# R2
ss_res <- sum((Y_test_temp - yhat_mean_temp)^2)
ss_tot <- sum((Y_test_temp - mean(Y_test_temp))^2)
R2     <- 1 - ss_res/ss_tot
cat("R-squared for temporal test:", round(R2, 4), "\n")

# Plot
plot(
  Y_test_temp, 
  yhat_mean_temp,
  xlab = "Observed Outflow (Y_test)",
  ylab = "Predicted Outflow (yhat_mean)",
  main = "Observed vs. Predicted Outflow test_year=16~19",
  pch  = 16
)
abline(0, 1, lty = 3)  # 45° line for reference


# feature importance
varcount_mat <- fit_temp$varcount # change
var_importance <- colMeans(varcount_mat)
var_imp_df <- data.frame(
  Variable   = colnames(X_train_temp), # change
  Importance = var_importance
)
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]
print(var_imp_df)
barplot(
  var_imp_df$Importance,
  names.arg  = var_imp_df$Variable,
  las        = 2,            # vertical axis labels
  cex.names  = 0.8,          # label size
  main       = "BART Variable Importance for temporal validity",
  ylab       = "Average Split Count"
)



# Model: spatial, ntree=5
fit_sp <- bart(
  x.train = X_train_sp,
  y.train = Y_train_sp,
  x.test  = X_test_sp,
  ntree   = 5,
  verbose = FALSE
)

# Results
yhat_mean_sp <- fit_sp$yhat.test.mean
yhat_samples_sp <- fit_sp$yhat.test
yhat_sd_sp <- apply(yhat_samples_sp, 2, sd)

results_temp <- data.frame(
  Observed   = Y_test_sp,       
  Predicted  = yhat_mean_sp,
  Uncertainty = yhat_sd_sp
)

# R2
ss_res <- sum((Y_test_sp - yhat_mean_sp)^2)
ss_tot <- sum((Y_test_sp - mean(Y_test_sp))^2)
R2     <- 1 - ss_res/ss_tot
cat("R-squared for spatial test:", round(R2, 4), "\n")

# Plot
plot(
  Y_test_sp, 
  yhat_mean_sp,
  xlab = "Observed Outflow (Y_test)",
  ylab = "Predicted Outflow (yhat_mean)",
  main = "Observed vs. Predicted Outflow test_spatial",
  pch  = 16
)
abline(0, 1, lty = 3)  # 45° line for reference

# feature importance
varcount_mat <- fit_sp$varcount # change
var_importance <- colMeans(varcount_mat)
var_imp_df <- data.frame(
  Variable   = colnames(X_train_sp), # change
  Importance = var_importance
)
var_imp_df <- var_imp_df[order(var_imp_df$Importance, decreasing = TRUE), ]
print(var_imp_df)
barplot(
  var_imp_df$Importance,
  names.arg  = var_imp_df$Variable,
  las        = 2,            # vertical axis labels
  cex.names  = 0.8,          # label size
  main       = "BART Variable Importance for spatial validity",
  ylab       = "Average Split Count"
)

