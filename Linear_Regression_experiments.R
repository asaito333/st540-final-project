##Run with only one reservoir characteristics i.e DOR

library(rjags)

set.seed(12)
res_chars <- read.csv('res_chars.csv')

# Read QI and QO
QI_df <- read.csv('QI.csv')
QO_df <- read.csv('QO.csv')

# Transfer data to matrix
QI <- as.matrix(sapply(QI_df[, -1], as.numeric))
QO <- as.matrix(sapply(QO_df[, -1], as.numeric))

colnames(QI) <- NULL
rownames(QI) <- NULL
colnames(QO) <- NULL
rownames(QO) <- NULL

# Standardize variables
QI_scaled <- scale(QI)
QO_scaled <- scale(QO)

N <- nrow(res_chars)  # number of reservoirs
Y <- ncol(QI)   # number of years per reservoir


# Simple Linear Regression ####Normal###########################################
# Data
data_jags <- list(
  N = N,
  Y = Y,
  QI = QI_scaled,
  QO = QO_scaled
)


model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      QO[i, y] ~ dnorm(beta0 + beta1*QI[i, y], tau)  #likelihood
      QO_pred[i, y] ~ dnorm(beta0 + beta1*QI[i, y], tau)  # posterior predictive
    }
    }
  # prior
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  tau  ~ dgamma(0.1, 0.1)
}
"



# Multivariate Linear Regression ###############################################

data_jags <- list(
  N = N,
  Y = Y,
  QI = QI_scaled,
  QO = QO_scaled,
  Capacity = scale(res_chars$Capacity),
  DOR = scale(res_chars$DOR),
  Depth = scale(res_chars$Depth),
  CatchArea = scale(res_chars$CatchArea),
  AvgRelease = scale(res_chars$AvgRelease)
)

model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      QO[i, y] ~ dnorm(beta0 + beta1*QI[i, y] + beta2*Capacity[i,] + beta3*DOR[i,] + beta4*Depth[i,] + beta5*CatchArea[i,] + beta6*AvgRelease[i,], tau)  #likelihood
      QO_pred[i, y] ~ dnorm(beta0 + beta1*QI[i, y] + beta2*Capacity[i,] + beta3*DOR[i,] + beta4*Depth[i,] + beta5*CatchArea[i,] + beta6*AvgRelease[i,], tau)  # posterior predictive
    }
    }
  # prior
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  beta2 ~ dnorm(0,0.001)
  beta3 ~ dnorm(0,0.001)
  beta4 ~ dnorm(0,0.001)
  beta5 ~ dnorm(0,0.001)
  beta6 ~ dnorm(0,0.001)
  tau  ~ dgamma(0.1, 0.1)
}
"


# Simple Linear Regression ###Log normal###-> fail########################################
# Data
data_jags <- list(
  N = N,
  Y = Y,
  QI = log(QI_scaled),
  QO = log(QO_scaled)
)


model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      QO[i, y] ~ dnorm(beta0 + beta1*QI[i, y], tau)  #likelihood
      log_QO_pred[i, y] ~ dnorm(beta0 + beta1*QI[i, y], tau)  # posterior predictive
    }
    }
  # prior
  beta0 ~ dnorm(0,0.001)
  beta1 ~ dnorm(0,0.001)
  tau  ~ dgamma(0.1, 0.1)
}
"


# One-way random effect ###########################################
# Data
data_jags <- list(
  N = N,
  Y = Y,
  QO = QO_scaled
)


model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      QO[i,y] ~ dnorm(beta[i], sig2_inv) # likelihood
      QO_pred[i, y] ~ dnorm(beta[i], sig2_inv)   # posterior predictive
    }
    }
  # random effects
  for (i in 1:N){
    beta[i] ~ dnorm(mu, tau2_inv)
  }
  # priors
  mu   ~ dnorm(0, 0.0001)
  sig2_inv ~ dgamma(0.1, 0.1)
  tau2_inv ~ dgamma(0.1, 0.1)
}
"

# One-way random effect with half-Cauchy prior ###########################################
# Data
data_jags <- list(
  N = N,
  Y = Y,
  QO = QO_scaled
)


model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      QO[i,y] ~ dnorm(beta[i], sig2_inv) # likelihood
      QO_pred[i, y] ~ dnorm(beta[i], sig2_inv)   # posterior predictive
    }
    }
  # random effects
  for (i in 1:N){
    beta[i] ~ dnorm(mu, tau2_inv)
  }
  # priors
  mu   ~ dnorm(0, 0.0001)
  sig2_inv <- pow(sigma1, -2)
  tau2_inv <- pow(sigma2, -2)
  sigma1 ~ dt(0, 1, 1)T(0, )
  sigma2 ~ dt(0, 1, 1)T(0, )
}
"


# One-way random effect with hierarchical random effect ###########################################
# Data
data_jags <- list(
  N = N,
  Y = Y,
  QO = QO_scaled,
  Capacity = scale(res_chars$Capacity),
  DOR = scale(res_chars$DOR),
  Depth = scale(res_chars$Depth),
  CatchArea = scale(res_chars$CatchArea),
  AvgRelease = scale(res_chars$AvgRelease)
)


model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      QO[i,y] ~ dnorm(beta[i], sig2_inv) # likelihood
      QO_pred[i, y] ~ dnorm(beta[i], sig2_inv)   # posterior predictive
    }
    }
  # random effects
  for (i in 1:N){
    beta[i] ~ dnorm(gamma0 + gamma1*Capacity[i,] + gamma2*DOR[i,] + gamma3*Depth[i,] + gamma4*CatchArea[i,] + gamma5*AvgRelease[i,], tau2_inv)
  }
  # priors
  mu   ~ dnorm(0, 0.0001)
  sig2_inv ~ dgamma(0.1, 0.1)
  tau2_inv ~ dgamma(0.1, 0.1)
  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  gamma2 ~ dnorm(0, 0.01)
  gamma3 ~ dnorm(0, 0.01)
  gamma4 ~ dnorm(0, 0.01)
  gamma5 ~ dnorm(0, 0.01)
}
"



###############################################


jags_model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 2, #quiet=FALSE,
  n.adapt = 5000  # Adaptation (includes burn-in)
)

# Step 2: Sample from posterior
samples <- coda.samples(
  jags_model,
  variable.names = c("gamma0", "delta0", "eta0"),
  n.iter = 5000  # Increased for better convergence
)

# Step 3: Posterior predictive checks (PPC)
samples_ppc <- coda.samples(
  jags_model,
  variable.names = c("QO_pred"),
  n.iter = 2000  # Can be smaller if only for visualization
)

# library(ggplot2)

# Convert to matrix and extract predicted values
samples_mat <- as.matrix(samples_ppc)
pred_vars <- grep("QO_pred", colnames(samples_mat), value = TRUE)

# Compute means across MCMC samples
pred_means <- colMeans(samples_mat[, pred_vars])
pred_matrix <- matrix(pred_means, nrow = N, ncol = Y)

# Flatten observed and predicted
obs <- as.vector(QO_scaled)
pred <- as.vector(pred_matrix)

# Plot
df_plot <- data.frame(Observed = obs, Predicted = pred)

nse <- 1 - sum((obs - pred)^2) / sum((obs - mean(obs))^2)
r_squared <- cor(obs, pred)^2

ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Posterior Predictive Check: Observed vs Predicted Outflow",
    subtitle = paste0("NSE = ", round(nse, 3), ",  R² = ", round(r_squared, 3)),
    x = "Observed Outflow",
    y = "Posterior Predictive Mean"
  ) +
  theme_minimal(base_size = 14)


# Scatter plot of 
QO_list <- as.list(c(QO))
QI_list <- as.list(c(QI))


plot(QO_list, QI_list, xlab="Outflow", ylab="inflow", 
     xlim=c(0, 15000), ylim=c(0, 30000))