load("./hw2/glm dat.RData")

## * Block MH
block_MH <- function(sample_size, burning_size,
                     init_parameter, block_idxes,
                     logposterior, logproposal,
                     proposal_func, args) {
  accept_rates <- matrix(0, sample_size, length(block_idxes))
  parameters <- matrix(NA, sample_size + 1, length(init_parameter))
  parameters[1, ] <- init_parameter
  ## old target density
  parameter_current <- parameters[1, ]
  pi.log <- logposterior(parameter_current, args)
  for (i in 1:sample_size) {
    for (j in seq_along(block_idxes)) {
      ## propose parameters
      parameter_star <- proposal_func(parameter_current, block_idxes[[j]], args)
      ## new target density
      pi.log.star <- logposterior(parameter_star, args)
      ## proposal density
      q_xn_xstar <- logproposal(parameter_current, parameter_star,
        block_idxes[[j]], args)
      q_xstar_xn <- logproposal(parameter_star, parameter_current,
        block_idxes[[j]], args)
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
logposterior <- function(parameter, args) {
  Y <- args$Y
  X <- args$X
  tau <- args$tau
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

logproposal <- function(para1, para2, idxes, args) {
  density <- +
    sum(dnorm(para1[1], para2[1], args$sd1, log = TRUE)) +
    sum(dnorm(para1[2], para2[2], args$sd2, log = TRUE)) +
    sum(dnorm(para1[3], para2[3], args$sd3, log = TRUE)) +
    sum(dnorm(para1[4], para2[4], args$sd4, log = TRUE))
  return(density)
}

proposal <- function(parameter, idxes, args) {
  parameter_star <- parameter
  for (idx in idxes) {
    if (idx == 1) {
      parameter_star[idx] <- rnorm(1, parameter[idx], args$sd1)
    } else if (idx == 2) {
      parameter_star[idx] <- rnorm(1, parameter[idx], args$sd2)
    } else if (idx == 3) {
      parameter_star[idx] <- rnorm(1, parameter[idx], args$sd3)
    } else {
      parameter_star[idx] <- rnorm(1, parameter[idx], args$sd4)
    }
  }
  return(parameter_star)
}

## * Tuning
sample_size <- 200000
burning_size <- 100000
res_mcmc <- block_MH(
  sample_size = 200000,
  burning_size = 100000,
  init_parameter = c(10, 1, 1, 1),
  block_idxes = list(1, 2, 3, 4),
  logposterior = logposterior,
  logproposal = logproposal,
  proposal_func = proposal,
  args = list(
    Y = Y, X = X, tau = sqrt(10),
    sd1 = 0.35, sd2 = 0.6, sd3 = 0.05, sd4 = 1.7))
sample_size <- dim(res_mcmc$samples)[1]

## * Evaluation
effective_size <- function(samples) {
  n <- length(samples)
  rho <- acf(samples, plot = FALSE)$acf
  rho <- rho[2:(n - 1)]
  rho[is.na(rho)] <- 0
  return(n / (1 + sum(rho)))
}

theta_mcmc <- res_mcmc$samples
accept_mcmc <- res_mcmc$accept_rates
par(mfrow = c(2, 2))
plot(1:sample_size, theta_mcmc[, 1], type = "l")
plot(1:sample_size, theta_mcmc[, 2], type = "l")
plot(1:sample_size, theta_mcmc[, 3], type = "l")
plot(1:sample_size, theta_mcmc[, 4], type = "l")
mtext(paste(
  paste(sprintf("accept rate %.3f", apply(accept_mcmc, 2, mean)),
    collapse = ', '),
  paste(sprintf("ESS %.3f", apply(theta_mcmc, 2, effective_size)),
    collapse = ', '),
  sep="\n"),
  side = 3, line = -3, cex = 1, outer = TRUE)

## * Posterior histogram
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
par(mfrow = c(1, 1))
plot(X_grid, Y_grid, type = "l")
points(X, Y)

## check normal assumption
residual.posmean <- Y - (1 + beta1.posmean * X) /
  (1 + beta2.posmean * exp(beta3.posmean * X))
plot(residual.posmean)
abline(h = 0, col = "red")

## * Importance Sampling
target_dense_IS <- function(theta, args) {
  Y <- args$Y
  X <- args$X
  tau <- args$tau
  beta1 <- theta[1]
  beta2 <- theta[2]
  beta3 <- theta[3]
  sigma2 <- theta[4]
  mu <- (1 + beta1 * X) / (1 + beta2 * exp(beta3 * X))
  pi.log <- dnorm(Y, mu, sqrt(sigma2)) *
    dnorm(beta1, 0, tau) *
    dnorm(beta2, 0, tau) *
    dnorm(beta3, 0, tau) *
    dgamma(sigma2, shape = 1, scale = 10)
  return(prod(pi.log))
}

proposal_sample_IS <- function(n, args) {
  cbind(rnorm(n, args$mu1, args$sd1),
    rnorm(n, args$mu2, args$sd2),
    rnorm(n, args$mu3, args$sd3),
    rgamma(n, shape = args$shape, scale = args$scale))
}

proposal_dense_IS <- function(theta, args) {
  dnorm(theta[, 1], args$mu1, args$sd1) *
    dnorm(theta[, 2], args$mu2, args$sd2) *
    dnorm(theta[, 3], args$mu3, args$sd3) *
    dgamma(theta[, 4], shape = args$shape, scale = args$scale)
}

importance_sampler <- function(Y, X, num_iter, args) {
  N <- round(exp(seq(2, 10, length.out = num_iter)))
  MC_est_IS <- matrix(NA, num_iter, 4)
  ## Importance sampling
  for (iter in 1:num_iter) {
    theta <- proposal_sample_IS(N[iter], args)
    g_X <- theta
    h_X <- apply(theta, 1, target_dense_IS, args = args)
    q_X <- proposal_dense_IS(theta, args)
    ratio_top <- g_X * h_X / q_X
    ratio_bot <- h_X / q_X
    ## update parameters
    MC_est_IS[iter, ] <- apply(ratio_top, 2, sum) / sum(ratio_bot)
  }
  return(MC_est_IS)
}

IS_estimate <- importance_sampler(Y, X, 100,
  list(Y = Y, X = X, tau = sqrt(5),
    mu1 = 6, mu2 = 3, mu3 = -0.4,
    sd1 = 2, sd2 = 2, sd3 = 2,
    shape = 1, scale = 10))
plot(IS_estimate[, 3])

target_dense_IS(c(6, 3, -0.4, 4), list(Y = Y, X = X, tau = sqrt(10)))
exp(logposterior(c(6, 3, -0.4, 4), list(Y = Y, X = X, tau = sqrt(10))))
