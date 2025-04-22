set.seed(12)

# Constants
N <- 200  # number of reservoirs
Y <- 20   # number of years per reservoir

# Reservoir characteristics
res_chars <- data.frame(
  Capacity = runif(N, 50, 300),
  DOR = runif(N, 0.1, 1),
  Depth = runif(N, 5, 50),
  CatchArea = runif(N, 100, 1000),
  AvgRelease = runif(N, 10, 300)
)

# Generate inflow and outflow data
QI <- matrix(runif(N * Y, 20, 200), nrow = N)
true_b0 <- rnorm(N, 50, 10)
true_b1 <- rnorm(N, 0.5, 0.1)
eta0 <- -1
eta1 <- 0.01

QO <- matrix(NA, nrow = N, ncol = Y)

for (i in 1:N) {
  for (y in 1:Y) {
    mu <- true_b0[i] + true_b1[i] * QI[i, y]
    sigma <- exp(eta0 + eta1 * QI[i, y])
    shape <- mu^2 / sigma^2
    rate <- mu / sigma^2
    QO[i, y] <- rgamma(1, shape = shape, rate = rate)
  }
}

library(rjags)

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
      log_mu[i,y] <- beta0[i] + beta1[i] * QI[i, y]
      log_sigma[i, y] <- eta0 + eta1 * QI[i, y]

      mu[i, y] <- exp(log_mu[i,y])
      sigma[i, y] <- exp(log_sigma[i, y])
      
      alpha[i, y] <- mu[i, y]^2 / sigma[i, y]^2
      lambda[i, y] <- mu[i, y] / sigma[i, y]^2
      
      QO[i, y] ~ dgamma(alpha[i, y], lambda[i, y])
      QO_pred[i, y] ~ dgamma(alpha[i, y], lambda[i, y])  # posterior predictive

    }

    beta0[i] ~ dnorm(gamma0 + gamma1*Capacity[i] + gamma2*DOR[i] + gamma3*Depth[i] + gamma4*CatchArea[i] + gamma5*AvgRelease[i], tau0)
    beta1[i] ~ dnorm(delta0 + delta1*Capacity[i] + delta2*DOR[i] + delta3*Depth[i] + delta4*CatchArea[i] + delta5*AvgRelease[i], tau1)
  }

  eta0 ~ dnorm(0, 1)
  eta1 ~ dnorm(0, 1)

  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)
  gamma2 ~ dnorm(0, 0.01)
  gamma3 ~ dnorm(0, 0.01)
  gamma4 ~ dnorm(0, 0.01)
  gamma5 ~ dnorm(0, 0.01)

  delta0 ~ dnorm(0, 0.01)
  delta1 ~ dnorm(0, 0.01)
  delta2 ~ dnorm(0, 0.01)
  delta3 ~ dnorm(0, 0.01)
  delta4 ~ dnorm(0, 0.01)
  delta5 ~ dnorm(0, 0.01)

  tau0 ~ dt(0, pow(2, -2), 1) T(0,)
  tau1 ~ dt(0, pow(2, -2), 1) T(0,)
}
"


# Step 1: Compile model with adaptation
jags_model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  n.chains = 2,
  n.adapt = 5000  # Adaptation (includes burn-in)
)

# Step 2: Optional extra burn-in (if needed)
# update(jags_model, 1000)  # Uncomment if diagnostics suggest it

# Step 3: Sample from posterior
samples <- coda.samples(
  jags_model,
  variable.names = c("gamma0", "delta0", "eta0", "tau0", "tau1"),
  n.iter = 20000  # Increased for better convergence
)

# Step 4: Posterior predictive checks (PPC)
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

ggplot(df_plot, aes(x = Observed, y = Predicted)) +
  geom_point(alpha = 0.7, color = "steelblue") +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray40") +
  labs(
    title = "Posterior Predictive Check: Observed vs Predicted Outflow",
    x = "Observed Outflow",
    y = "Posterior Predictive Mean"
  ) +
  theme_minimal(base_size = 14)

