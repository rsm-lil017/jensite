# Load packages
library(dplyr)
library(tidyr)
library(ggplot2)
library(MASS)
library(coda)

# Load Data
conjoint_data <- read.csv("~/Desktop/MSBA/SP/MGTA495b/hw3/conjoint_data.csv")

# Task 1: Prepare the data

# Reshape and prep the data for MNL: include respondent (i), alternative (j), and covariate (k)
conjoint_data <- conjoint_data %>%
  mutate(
    alt_id = rep(1:3, times = nrow(conjoint_data) / 3),
    brand_N = ifelse(brand == "N", 1, 0),
    brand_P = ifelse(brand == "P", 1, 0),
    ad_Yes = ifelse(ad == "Yes", 1, 0),
    resp = as.integer(resp),
    task = as.integer(task)
  ) %>%
  arrange(resp, task, alt_id)
conjoint_data <- conjoint_data %>%
  mutate(brand_N = ifelse(brand == "N", 1, 0),
         brand_P = ifelse(brand == "P", 1, 0),
         ad_Yes = ifelse(ad == "Yes", 1, 0))

X <- as.matrix(dplyr::select(conjoint_data, brand_N, brand_P, ad_Yes, price))
y <- conjoint_data$choice

# Task 2: Define Log-likelihood for MNL
softmax_utilities <- function(X, beta) {
  utilities <- X %*% matrix(beta, ncol = 1)
  N <- nrow(utilities)
  if (N %% 3 != 0) stop("Number of rows in X must be a multiple of 3")
  utilities <- matrix(utilities, ncol = 3, byrow = TRUE)
  exp_u <- exp(utilities - apply(utilities, 1, max))
  probs <- exp_u / rowSums(exp_u)
  return(as.vector(t(probs)))
}

log_likelihood <- function(beta, X, y) {
  probs <- softmax_utilities(X, beta)
  -sum(y * log(probs + 1e-12))
}

# Task 3: MLE Estimation
init_beta <- rep(0, ncol(X))
mle_result <- optim(par = init_beta, fn = log_likelihood, X = X, y = y, hessian = TRUE)
mle_beta <- mle_result$par
se <- sqrt(diag(solve(mle_result$hessian)))
z <- qnorm(0.975)
ci_lower <- mle_beta - z * se
ci_upper <- mle_beta + z * se

mle_results <- data.frame(
  Parameter = c("brand_N", "brand_P", "ad_Yes", "price"),
  Estimate = mle_beta,
  StdErr = se,
  CI_Lower = ci_lower,
  CI_Upper = ci_upper
)
print(mle_results)

# Task 4: Bayesian Estimation via MCMC
log_prior <- function(beta) {
  sum(dnorm(beta[1:3], 0, 5, log = TRUE)) + dnorm(beta[4], 0, 1, log = TRUE)
}

log_posterior <- function(beta, X, y) {
  -log_likelihood(beta, X, y) + log_prior(beta)
}

n_iter <- 11000
burn_in <- 1000
beta_chain <- matrix(0, nrow = n_iter, ncol = 4)
beta_current <- rep(0, 4)
log_post_current <- log_posterior(beta_current, X, y)
proposal_sd <- c(0.05, 0.05, 0.05, 0.005)

for (t in 1:n_iter) {
  beta_proposal <- beta_current + rnorm(4, 0, proposal_sd)
  log_post_proposal <- log_posterior(beta_proposal, X, y)
  accept_ratio <- exp(log_post_proposal - log_post_current)
  if (runif(1) < accept_ratio) {
    beta_current <- beta_proposal
    log_post_current <- log_post_proposal
  }
  beta_chain[t, ] <- beta_current
}

beta_chain_burned <- beta_chain[(burn_in+1):n_iter, ]

posterior_means <- colMeans(beta_chain_burned)
posterior_std <- apply(beta_chain_burned, 2, sd)
posterior_ci_lower <- apply(beta_chain_burned, 2, quantile, probs = 0.025)
posterior_ci_upper <- apply(beta_chain_burned, 2, quantile, probs = 0.975)

bayes_results <- data.frame(
  Parameter = c("brand_N", "brand_P", "ad_Yes", "price"),
  PosteriorMean = posterior_means,
  StdDev = posterior_std,
  CI_Lower = posterior_ci_lower,
  CI_Upper = posterior_ci_upper
)
print(bayes_results)

# Task 5: Trace plot and histogram for 'price'
par(mfrow = c(1, 2))
plot(beta_chain_burned[, 4], type = 'l', main = "Trace plot for price", xlab = "Iteration", ylab = "Beta (price)")
hist(beta_chain_burned[, 4], breaks = 30, probability = TRUE, main = "Posterior distribution for price", xlab = "Beta (price)")
lines(density(beta_chain_burned[, 4]), col = 'blue')
