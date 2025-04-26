#Run with only one reservoir characteristics i.e DOR

library(rjags)

set.seed(12)
res_chars <- read.csv('res_chars.csv')
QI_df <- read.csv('QI.csv')
QO_df <- read.csv('QO.csv')

QI <- as.matrix(sapply(QI_df[, -1], as.numeric))
QO <- as.matrix(sapply(QO_df[, -1], as.numeric))

colnames(QI) <- NULL
rownames(QI) <- NULL
colnames(QO) <- NULL
rownames(QO) <- NULL

N <- nrow(res_chars)  # number of reservoirs
Y <- ncol(QI)   # number of years per reservoir

QI_scaled <- scale(QI)
QO_scaled <- scale(QO)

data_jags <- list(
  N = N,
  Y = Y,
  QI = QI,
  QO = QO,
  Capacity = res_chars$Capacity,
  DOR = res_chars$DOR,
  Depth = res_chars$Depth,
  CatchArea = res_chars$CatchArea,
  AvgRelease = res_chars$AvgRelease
)



model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      mu[i, y] <- beta0[i] + beta1[i] * QI[i, y]
      sigma[i, y] <- eta0 + eta1 * QI[i, y]
      
      
      alpha[i, y] <- mu[i, y]^2 / (sigma[i, y]^2 + 0.01)
      lambda[i, y] <- mu[i, y] / (sigma[i, y]^2 + 0.01)
      
      QO[i, y] ~ dgamma(alpha[i, y], lambda[i, y])  #likelihood
      QO_pred[i, y] ~ dgamma(alpha[i, y], lambda[i, y])  # posterior predictive
    
    }

    beta0[i] ~ dnorm(gamma0 + gamma1*Capacity[i] + gamma2*DOR[i] + gamma3*Depth[i] + gamma4*CatchArea[i] + gamma5*AvgRelease[i] , 0.01)
    beta1[i] ~ dnorm(delta0 + delta1*Capacity[i] + delta2*DOR[i] + delta3*Depth[i] + delta4*CatchArea[i] + delta5*AvgRelease[i] , 0.01)
  }

  eta0 ~ dgamma(0.5, 1)
  eta1 ~ dgamma(0.5, 1)

  gamma0 ~ dgamma(0.5, 1)
  gamma1 ~ dgamma(0.5, 1)
  gamma2 ~ dgamma(0.5, 1)
  gamma3 ~ dgamma(0.5, 1)
  gamma4 ~ dgamma(0.5, 1)
  gamma5 ~ dgamma(0.5, 1)

  delta0 ~ dgamma(0.5, 1)
  delta1 ~ dgamma(0.5, 1)
  delta2 ~ dgamma(0.5, 1)
  delta3 ~ dgamma(0.5, 1)
  delta4 ~ dgamma(0.5, 1)
  delta5 ~ dgamma(0.5, 1)
}
"

# Step 1: Compile model with adaptation
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

# Diagnostics
summary(samples)
gelman.diag(samples)  # Check R-hat
traceplot(samples)    # Visual inspection

library(ggplot2)

# Convert to matrix and extract predicted values
samples_mat <- as.matrix(samples_ppc)
pred_vars <- grep("QO_pred", colnames(samples_mat), value = TRUE)

# Compute means across MCMC samples
pred_means <- colMeans(samples_mat[, pred_vars])
pred_matrix <- matrix(pred_means, nrow = N, ncol = Y)

# Flatten observed and predicted
obs <- as.vector(QO)
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

