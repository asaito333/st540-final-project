##Run with only one reservoir characteristics i.e DOR

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

log_QO <- log(QO)

N <- nrow(res_chars)  # number of reservoirs
Y <- ncol(QI)   # number of years per reservoir

data_jags <- list(
  N = N,
  Y = Y,
  QI = QI,
  log_QO = log_QO,
  DOR = res_chars$DOR
)

model_string <- "
model {
  for (i in 1:N) {
    for (y in 1:Y) {
      mu_log[i, y] <- beta0[i] + beta1[i] * QI[i, y]      # Mean of log-transformed data
      log_sigma[i, y] <- eta0 + eta1 * QI[i, y]           # Log of SD
      tau[i, y] <- pow(1 / exp(log_sigma[i, y]), 2)       # Precision

      log_QO[i, y] ~ dnorm(mu_log[i, y], tau[i, y])

      log_QO_pred[i, y] ~ dnorm(mu_log[i, y], tau[i, y])  # Prediction in log-space
      QO_pred[i, y] <- exp(log_QO_pred[i, y])             # Back-transform to original scale

    }

    beta0[i] ~ dnorm(gamma0 +  gamma1*DOR[i] , tau0)
    beta1[i] ~ dnorm(delta0 +  delta1*DOR[i] , tau1)
  }

  eta0 ~ dnorm(0, 1)
  eta1 ~ dnorm(0, 1)

  gamma0 ~ dnorm(0, 0.01)
  gamma1 ~ dnorm(0, 0.01)

  delta0 ~ dnorm(0, 0.01)
  delta1 ~ dnorm(0, 0.01)

  tau0 ~ dt(0, pow(2, -2), 1) T(0,) #Half cauchy (0,2)
  tau1 ~ dt(0, pow(2, -2), 1) T(0,)
}
"
inits <- function() {
  list(
    beta0 = rnorm(N, 0, 1),
    beta1 = rnorm(N, 0, 1),
    eta0 = 0,
    eta1 = 0,
    gamma0 = 0,
    gamma1 = 0,
    delta0 = 0,
    delta1 = 0,
    tau0 = 1,
    tau1 = 1
  )
}



# Step 1: Compile model with adaptation
jags_model <- jags.model(
  textConnection(model_string),
  data = data_jags,
  inits = inits,
  n.chains = 2,
  n.adapt = 5000  # Adaptation (includes burn-in)
)

# Step 2: Sample from posterior
samples <- coda.samples(
  jags_model,
  variable.names = c("gamma0", "delta0", "eta0", "tau0", "tau1"),
  n.iter = 20000  # Increased for better convergence
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
pred_means <- colMeans(exp(samples_mat[, pred_vars]))
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

