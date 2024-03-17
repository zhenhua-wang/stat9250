load("./hw2/glm dat.RData")

## * Block MH
loglikelihood <- function(Y, X, parameter) {
  beta1 <- parameter[1]
  beta2 <- parameter[2]
  beta3 <- parameter[3]
  sigma2 <- parameter[4]
  mu <- (1 + beta1 * X) / (1 + beta2 * exp(beta3 * X))
  pi.log <- sum(dnorm(Y, mu, sqrt(sigma2), log = TRUE))
  return(pi.log)
}

block_MH <- function(sample_size, burning_size,
                     X, Y,
                     init_parameter, block_idxes,
                     loglikelihood, proposal_func) {
  parameters <- matrix(NA, sample_size + 1, length(init_parameter))
  parameters[1, ] <- init_parameter
  accept_rates <- matrix(0, sample_size, length(block_idxes))
  for (i in 1:sample_size) {
    ## old density
    parameter_current <- parameters[i, ]
    pi.log <- loglikelihood(Y, X, parameter_current)
    for (j in seq_along(block_idxes)) {
      ## propose parameters
      parameter_star <- proposal_func(parameter_current, block_idxes[[j]])
      ## new density
      pi.log.star <- loglikelihood(Y, X, parameter_star)
      ## update
      alpha.log <- min(0, pi.log.star - pi.log)
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
    accept_rates = accept_rates))
}

proposal <- function(parameter, idxes) {
  parameter_star <- parameter
  for (idx in idxes) {
    if (idx == 4) {
      parameter_star[idx] <- rgamma(1, shape = 1, scale = 10)
    } else {
      parameter_star[idx] <- rnorm(1, 0, tau)
    }
  }
  return(parameter_star)
}

sample_size <- 100000
burning_size <- 50000
tau <- sqrt(10)
res_mcmc <- block_MH(
  sample_size = sample_size, burning_size = burning_size,
  X = X, Y = Y,
  init_parameter = c(0.1, 0.1, 0.1, 0.1),
  block_idxes = list(1, 2, 3, 4),
  loglikelihood = loglikelihood,
  proposal_func = proposal)

theta_mcmc <- res_mcmc$samples
accept_mcmc <- res_mcmc$accept_rates
par(mfrow = c(2, 2))
plot(1:(sample_size - burning_size), theta_mcmc[, 1], type = "l")
plot(1:(sample_size - burning_size), theta_mcmc[, 2], type = "l")
plot(1:(sample_size - burning_size), theta_mcmc[, 3], type = "l")
plot(1:(sample_size - burning_size), theta_mcmc[, 4], type = "l")

apply(accept_mcmc, 2, mean)

## * Importance sampling
N <- round(exp(seq( 2, 10, length.out=100)))
MC_est_IS <- rep( NA, length(N))
MC_se_IS <- rep( NA, length(N))
## Importance sampling
sd <- sqrt(0.5); f_X <- 1/10

for(iter in 1:length(N)){
  X <- rnorm(N[iter], mean = 5, sd = sd)
  g_X <- exp(-2*abs(X-5))*10
  q_X <- dnorm(X, mean = 5, sd = sd)
  ratio <- g_X*f_X/q_X
  MC_est_IS[iter] <- mean(ratio)
  MC_se_IS[iter] <- sqrt(var(ratio)/N[iter])
}
