# Load packages
library(rstan)
library(tidyverse)

# Read in data
res_chars <- read.csv("C:/Users/Emma/Downloads/res_chars.csv")
QI_df <- read.csv("C:/Users/Emma/Downloads/QI.csv")
QO_df <- read.csv("C:/Users/Emma/Downloads/QO.csv")

# Prepare data
QI <- as.matrix(QI_df[, -1])  # remove ReservoirID column
QO <- as.matrix(QO_df[, -1])

# Dimensions
N <- nrow(QI)  # number of reservoirs
T <- ncol(QI)  # number of years

# Stan model for simplified GEV (global alpha, beta; fixed xi)
gev_stan_model <- "
functions {
  real gev_lpdf(real y, real mu, real sigma, real xi) {
    real z = (y - mu) / sigma;
    if (fabs(xi) > 1e-6) {
      real t = 1 + xi * z;
      if (t <= 0) return negative_infinity();
      return -log(sigma) - (1 + 1/xi) * log(t) - pow(t, -1/xi);
    } else {
      return -log(sigma) - z - exp(-z);
    }
  }
}
data {
  int<lower=1> N;           // number of reservoirs
  int<lower=1> T;           // number of years
  matrix[N, T] QI;          // inflow
  matrix[N, T] QO;          // outflow (response)
}
parameters {
  real alpha;
  real beta;
  real<lower=0> sigma;
}
model {
  // Priors
  alpha ~ normal(0, 10);
  beta ~ normal(0, 10);
  sigma ~ normal(10, 5);

  // Likelihood
  for (i in 1:N)
    for (t in 1:T)
      QO[i, t] ~ gev(alpha + beta * log(QI[i, t]), sigma, 0.1); // fixed xi = 0.1
}
"

# Compile Stan model
stan_model <- stan_model(model_code = gev_stan_model)

# Prepare data for Stan
data_list <- list(
  N = N,
  T = T,
  QI = QI,
  QO = QO
)

# Fit the model
fit <- sampling(
  object = stan_model,
  data = data_list,
  iter = 1500,
  chains = 4,
  seed = 42,
  control = list(adapt_delta = 0.99)
)
