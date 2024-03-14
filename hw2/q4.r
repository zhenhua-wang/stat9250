load("./hw2/glm dat.RData")

tau <- sqrt(10)
sample_size <- 10000
beta1 <- 0.1
beta2 <- 0.1
beta3 <- 0.1
sigma2 <- 0.1
theta_mcmc <- matrix(NA, sample_size, 4)
accept_mcmc <- matrix(0, sample_size, 3)
for (i in 1:sample_size) {
  ## old
  mu <- (1 + beta1 * X) / (1 + beta2 * exp(beta3 * X))
  pi.log <- sum(dnorm(mu, sqrt(sigma2), log = TRUE))
  ## update beta
  beta1.star <- rnorm(1, 0, tau)
  mu.star <- (1 + beta1.star * X) / (1 + beta2 * exp(beta3 * X))
  pi.log.star <- sum(dnorm(mu.star, sqrt(sigma2), log = TRUE))
  alpha.log <- min(0, pi.log.star - pi.log)
  U.log <- log(runif(1))
  if (U.log < alpha.log) {
    theta_mcmc[i, 1] <- beta1.star
    beta1 <- beta1.star
    pi.log <- pi.log.star
    accept_mcmc[i, 1] <- 1
  } else {
    theta_mcmc[i, 1] <- beta1
  }
  ## update beta2, 3
  beta2.star <- rnorm(1, 0, tau)
  beta3.star <- rnorm(1, 0, tau)
  mu.star <- (1 + beta1 * X) / (1 + beta2.star * exp(beta3.star * X))
  pi.log.star <- sum(dnorm(mu.star, sqrt(sigma2), log = TRUE))
  alpha.log <- min(0, pi.log.star - pi.log)
  U.log <- log(runif(1))
  if (U.log < alpha.log) {
    theta_mcmc[i, 2] <- beta2.star
    theta_mcmc[i, 3] <- beta3.star
    beta2 <- beta2.star
    beta3 <- beta3.star
    pi.log <- pi.log.star
    accept_mcmc[i, 2] <- 1
  } else {
    theta_mcmc[i, 2] <- beta2
    theta_mcmc[i, 3] <- beta3
  }
  ## update sigma2
  sigma2.star <- rgamma(1, shape = 1, scale = 10)
  mu.star <- (1 + beta1 * X) / (1 + beta2 * exp(beta3 * X))
  pi.log.star <- sum(dnorm(mu.star, sqrt(sigma2.star), log = TRUE))
  alpha.log <- min(0, pi.log.star - pi.log)
  U.log <- log(runif(1))
  if (U.log < alpha.log) {
    theta_mcmc[i, 4] <- sigma2.star
    sigma2 <- sigma2.star
    pi.log <- pi.log.star
    accept_mcmc[i, 3] <- 1
  } else {
    theta_mcmc[i, 4] <- sigma2
  }
}

par(mfrow = c(2, 2))
plot(1:sample_size, theta_mcmc[, 1], type = "l")
plot(1:sample_size, theta_mcmc[, 2], type = "l")
plot(1:sample_size, theta_mcmc[, 3], type = "l")
plot(1:sample_size, theta_mcmc[, 4], type = "l")

apply(accept_mcmc, 2, mean)
