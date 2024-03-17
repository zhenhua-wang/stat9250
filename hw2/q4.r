load("./hw2/glm dat.RData")

## * Block MH
block_MH <- function(X, Y,
                     sample_size, burning_size,
                     init_parameter, block_idxes,
                     logposterior, logproposal,
                     proposal_func) {
  parameters <- matrix(NA, sample_size + 1, length(init_parameter))
  parameters[1, ] <- init_parameter
  accept_rates <- matrix(0, sample_size, length(block_idxes))
  ## old target density
  parameter_current <- parameters[1, ]
  pi.log <- logposterior(Y, X, parameter_current)
  for (i in 1:sample_size) {
    for (j in seq_along(block_idxes)) {
      ## propose parameters
      parameter_star <- proposal_func(parameter_current, block_idxes[[j]])
      ## new target density
      pi.log.star <- logposterior(Y, X, parameter_star)
      ## proposal density
      q_xn_xstar <- logproposal(parameter_current,
        parameter_star, block_idxes[[j]])
      q_xstar_xn <- logproposal(parameter_star,
        parameter_current, block_idxes[[j]])
      ## update
      alpha.log <- min(0, pi.log.star + q_xn_xstar - pi.log - q_xstar_xn)
      U.log <- log(runif(1))
      if (U.log < alpha.log) {
        parameters[i + 1, ] <- parameter_star
        parameter_current <- parameter_star
        pi.log <- pi.log.star
        accept_rates[i, j] <- 1
      } else {
        parameters[i + 1, ] <- parameter_current
      }
    }
  }
  return(list(
    samples = parameters[(burning_size + 2):(sample_size + 1), ],
    accept_rates = accept_rates[(burning_size + 1):sample_size, ]))
}

## * Density
logposterior <- function(Y, X, parameter) {
  tau <- sqrt(10)
  beta1 <- parameter[1]
  beta2 <- parameter[2]
  beta3 <- parameter[3]
  sigma2 <- parameter[4]
  ## reject invalid sigma2
  if (sigma2 < 0) {
    return(log(0))
  }
  mu <- (1 + beta1 * X) / (1 + beta2 * exp(beta3 * X))
  pi.log <- sum(dnorm(Y, mu, sqrt(sigma2), log = TRUE)) +
    sum(dnorm(beta1, 0, tau, log = TRUE)) +
    sum(dnorm(beta2, 0, tau, log = TRUE)) +
    sum(dnorm(beta3, 0, tau, log = TRUE)) +
    sum(dgamma(sigma2, shape = 1, scale = 10, log = TRUE))
  return(pi.log)
}

logproposal <- function(para1, para2, idxes) {
  density <- +
    sum(dnorm(
      para1[1], para2[1],
      proposal_hyperparam$sd1, log = TRUE)) +
    sum(dnorm(
      para1[2], para2[2],
      proposal_hyperparam$sd2, log = TRUE)) +
    sum(dnorm(
      para1[3], para2[3],
      proposal_hyperparam$sd3, log = TRUE)) +
    sum(dnorm(
      para1[4], para2[4],
      proposal_hyperparam$sd4, log = TRUE))
  return(density)
}

proposal <- function(parameter, idxes) {
  parameter_star <- parameter
  for (idx in idxes) {
    if (idx == 1) {
      parameter_star[idx] <- rnorm(
        1, parameter[idx], proposal_hyperparam$sd1)
    } else if (idx == 2) {
      parameter_star[idx] <- rnorm(
        1, parameter[idx], proposal_hyperparam$sd2)
    } else if (idx == 3) {
      parameter_star[idx] <- rnorm(
        1, parameter[idx], proposal_hyperparam$sd3)
    } else {
      parameter_star[idx] <- rnorm(
        1, parameter[idx], proposal_hyperparam$sd4)
    }
  }
  return(parameter_star)
}

## * Init parameter
## pick median as init
beta1 <- - 10#1 / X[166]
beta2 <- 1#1 / Y[166] - 1
beta3 <- 1#log((1 + beta1 - Y[166]) / (Y[166] * beta2))
sig2 <- 1#sd(Y)

## * Tuning
sample_size <- 200000
burning_size <- 100000
proposal_hyperparam <- list(
  sd1 = 0.35, sd2 = 0.26, sd3 = 0.07, sd4 = 3.2)
res_mcmc <- block_MH(
  X = X, Y = Y,
  sample_size = sample_size,
  burning_size = burning_size,
  init_parameter = c(beta1, beta2, beta3, sig2),
  block_idxes = list(1, 2, 3, 4),
  logposterior = logposterior,
  logproposal = logproposal,
  proposal_func = proposal)

## * Result
theta_mcmc <- res_mcmc$samples
accept_mcmc <- res_mcmc$accept_rates
par(mfrow = c(2, 2))
plot(1:(sample_size - burning_size), theta_mcmc[, 1], type = "l")
plot(1:(sample_size - burning_size), theta_mcmc[, 2], type = "l")
plot(1:(sample_size - burning_size), theta_mcmc[, 3], type = "l")
plot(1:(sample_size - burning_size), theta_mcmc[, 4], type = "l")
mtext(paste(sprintf("accept rate %.3f", apply(accept_mcmc, 2, mean)),
  collapse = ', '),
  side = 3, line = -2, cex = 1, outer = TRUE)

apply(accept_mcmc, 2, mean)

par(mfrow = c(2, 2))
hist(theta_mcmc[, 1])
hist(theta_mcmc[, 2])
hist(theta_mcmc[, 3])
hist(theta_mcmc[, 4])

## * Prediction using posterior mean
beta1.posmean <- mean(theta_mcmc[, 1])
beta2.posmean <- mean(theta_mcmc[, 2])
beta3.posmean <- mean(theta_mcmc[, 3])
beta4.posmean <- mean(theta_mcmc[, 4])
X_grid <- seq(0.001, 10, length.out = 10000)
Y_grid <- (1 + beta1.posmean * X_grid) /
  (1 + beta2.posmean * exp(beta3.posmean * X_grid))
plot(X_grid, Y_grid, type = "l")
points(X, Y)
