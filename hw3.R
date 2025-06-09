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

# Task 5: Trace plot and histogram for 'brand_N' (Netflix)
par(mfrow = c(1, 2))
plot(beta_chain_burned[, 1], type = 'l', main = "Trace plot for brand_N", xlab = "Iteration", ylab = "Beta (brand_N)")
hist(beta_chain_burned[, 1], breaks = 30, probability = TRUE, main = "Posterior distribution for brand_N", xlab = "Beta (brand_N)")
lines(density(beta_chain_burned[, 1]), col = 'blue')")
lines(density(beta_chain_burned[, 4]), col = 'blue')"




# Install needed packages
install.packages("flexmix")
install.packages("tidyverse")
library(flexmix)
library(tidyverse)

# Read your data
yogurt <- read.csv("~/Desktop/MSBA/SP/MGTA495b/hw4/yogurt_data.csv")
keydrivers <- read.csv("~/Desktop/MSBA/SP/MGTA495b/hw4/data_for_drivers_analysis.csv")
penguins <- read.csv("~/Desktop/MSBA/SP/MGTA495b/hw4/palmer_penguins.csv")

# Step 1: Long format for choice
yogurt_long <- yogurt %>%
  pivot_longer(
    cols = starts_with("y"),
    names_to = "alt",
    names_pattern = "y(\\d+)",
    values_to = "choice"
  ) %>%
  mutate(
    alt = as.integer(alt)
  )

# Step 2: Get price and featured for each alt (using pivot_longer & join)
price_long <- yogurt %>%
  pivot_longer(
    cols = starts_with("p"),
    names_to = "alt",
    names_pattern = "p(\\d+)",
    values_to = "price"
  ) %>%
  mutate(
    alt = as.integer(alt)
  ) %>%
  select(id, alt, price)

featured_long <- yogurt %>%
  pivot_longer(
    cols = starts_with("f"),
    names_to = "alt",
    names_pattern = "f(\\d+)",
    values_to = "featured"
  ) %>%
  mutate(
    alt = as.integer(alt)
  ) %>%
  select(id, alt, featured)

# Step 3: Join them all by id and alt
yogurt_long <- yogurt_long %>%
  left_join(price_long, by = c("id", "alt")) %>%
  left_join(featured_long, by = c("id", "alt"))

library(flexmix)
# successes = choice, failures = 1-choice
yogurt_long$success <- yogurt_long$choice
yogurt_long$failure <- 1 - yogurt_long$choice

results <- list()
for (k in 2:5) {
  fm <- flexmix(cbind(success, failure) ~ price + featured | id,
                data = yogurt_long,
                k = k,
                model = FLXMRglm(family = "binomial"))
  results[[paste0("k", k)]] <- fm
  cat("\n-----\nk =", k, ", BIC =", BIC(fm), "\n")
}
best_model <- results[["k2"]]  # or whichever k has the lowest BIC
parameters(best_model)
summary(best_model)




install.packages("mlogit")
library(mlogit)

# Prepare data for mlogit
mlogit_data <- mlogit.data(
  yogurt_long,
  choice = "success",    # or whatever your choice column is (should be 0/1 or TRUE/FALSE)
  shape = "long",
  id.var = "id",
  alt.var = "alt"
)

# Fit standard MNL (conditional logit)
mnl_model <- mlogit(success ~ price + featured | 0, data = mlogit_data)
summary(mnl_model)


set.seed(42)
n <- 100
x1 <- runif(n, -3, 3)
x2 <- runif(n, -3, 3)
boundary <- sin(4 * x1) + x1
y <- ifelse(x2 > boundary, 1, 0)
dat <- data.frame(x1 = x1, x2 = x2, y = as.factor(y))

library(ggplot2)
ggplot(dat, aes(x = x1, y = x2, color = y)) +
  geom_point(size = 2) +
  stat_function(fun = function(x) sin(4 * x) + x, color = "black", linetype = "dashed") +
  labs(title = "Synthetic Data with Wiggly Boundary") +
  theme_minimal()

set.seed(7)
x1_test <- runif(n, -3, 3)
x2_test <- runif(n, -3, 3)
boundary_test <- sin(4 * x1_test) + x1_test
y_test <- ifelse(x2_test > boundary_test, 1, 0)
dat_test <- data.frame(x1 = x1_test, x2 = x2_test, y = as.factor(y_test))

euclidean_dist <- function(x, y) sqrt(sum((x - y)^2))

knn_predict <- function(train, test, k) {
  pred <- numeric(nrow(test))
  for (i in 1:nrow(test)) {
    dists <- apply(train[, c("x1", "x2")], 1, euclidean_dist, y = as.numeric(test[i, c("x1", "x2")]))
    neighbors <- order(dists)[1:k]
    pred[i] <- names(sort(table(train$y[neighbors]), decreasing = TRUE))[1]
  }
  as.factor(pred)
}

library(class)
# Use class::knn
knn_builtin <- function(k) {
  class::knn(
    train = dat[, c("x1", "x2")],
    test = dat_test[, c("x1", "x2")],
    cl = dat$y,
    k = k
  )
}

accuracy_manual <- numeric(30)
accuracy_builtin <- numeric(30)

for (k in 1:30) {
  pred_manual <- knn_predict(dat, dat_test, k)
  pred_builtin <- knn_builtin(k)
  accuracy_manual[k] <- mean(pred_manual == dat_test$y)
  accuracy_builtin[k] <- mean(pred_builtin == dat_test$y)
}

plot(1:30, accuracy_manual * 100, type = "b", pch = 16, col = "blue",
     xlab = "k (number of neighbors)", ylab = "Accuracy (%)",
     main = "KNN Test Accuracy vs. k")
lines(1:30, accuracy_builtin * 100, type = "b", pch = 1, col = "red")
legend("bottomright", legend = c("Manual KNN", "Built-in KNN"), col = c("blue", "red"), pch = c(16, 1))

optimal_k <- which.max(accuracy_builtin)
cat(sprintf("Optimal k (highest accuracy): %d (%.1f%%)\n", optimal_k, accuracy_builtin[optimal_k] * 100))


