## Build the sampler function
# Define the model and priors
likelihood <- function(beta1, beta2, beta3, sigma2, X, Y) {
  y_pred <- (1 + beta1 * X)/(1 + beta2 * exp(beta3 * X))
  return(dnorm(Y, mean = y_pred, sd = sqrt(sigma2), log = TRUE))
}

prior_beta <- function(beta, tau2 = 10) {
  dnorm(beta, mean = 0, sd = sqrt(tau2), log = TRUE) # Non-informative prior
}

prior_sigma2 <- function(sigma2, a, b) {
  dgamma(sigma2, shape = a, scale = b, log = TRUE)
}


## Using gibbs sampler for Beta 1 and Matropolis-Hastings for rest of the parameters
## Treat the Beta 2 and Beta 3 as the block
mix_sampler <- function(X, Y, nburn, nsim, nthin, tau2 = 10, a = 1, b = 10) {
  m <- nrow(X)
  p <- ncol(X)
  inv.sigma.beta <- 1/tau2

  # Initial values
  Beta_1 <- 10
  Beta_2 <- 1
  Beta_3 <- 1 
  sigma2 <- 1

  # Chain containers
  Beta_1.chain <- array(0, dim = c(1, nsim / nthin))
  Beta_2.chain <- array(0, dim = c(1, nsim / nthin))
  Beta_3.chain <- array(0, dim = c(1, nsim / nthin))
  sigma2.chain <- rep(0, nsim / nthin)
  acceptance_beta.chain <- rep(0, nsim / nthin)
  acceptance_sigma.chain <- rep(0, nsim / nthin)

  for (index in 1:(nsim + nburn)) {
    if (index %% 10000 == 0) cat(index, "\n")

    
    # Update Beta 2 and Beta 3

    # Propose new values
    proposed_beta_2 <- rnorm(1, mean = Beta_2, sd = 0.1)
    proposed_beta_3 <- rnorm(1, mean = Beta_3, sd = 0.1)
    
    # Compute acceptance ratio
    log_acceptance_ratio_beta <- sum(likelihood(Beta_1, proposed_beta_2, proposed_beta_3, sigma2, X, Y)) +
      prior_beta(proposed_beta_2) + prior_beta(proposed_beta_3) -
      sum(likelihood(Beta_1, Beta_2, Beta_3, sigma2, X, Y)) -
      prior_beta(Beta_2) - prior_beta(Beta_3)

    # Accept or reject for Beta 2 and Beta 3
    if(log(runif(1)) < log_acceptance_ratio_beta) {
      Beta_2 <- proposed_beta_2
      Beta_3 <- proposed_beta_3
      acceptance_beta <- 1
    } else {
      Beta_2 <- Beta_2
      Beta_3 <- Beta_3
      acceptance_beta <- 0
    }
    
    # Update Sigma2
    # Propose sigma2
    proposed_sigma2 <- rgamma(1, shape = a, scale = b)

    # Compute acceptance ratio
    log_acceptance_ratio_sigma2 <- sum(likelihood(Beta_1, Beta_2, Beta_3, proposed_sigma2, X, Y)) +
      prior_sigma2(proposed_sigma2, a, b) -
      sum(likelihood(Beta_1, Beta_2, Beta_3, sigma2, X, Y)) -
      prior_sigma2(sigma2, a, b)

    if(log(runif(1)) < log_acceptance_ratio_sigma2) {
      sigma2 <- proposed_sigma2
      acceptance_sigma <- 1
    } else {
      sigma2 <- sigma2
      acceptance_sigma <- 0
    }


    # Define new parameters
    eta <- 1/(1 + Beta_2 * exp(Beta_3 * X))
    zeta <- X * eta

    # Update Beta 1
    var.Beta_1 <- solve(t(zeta) %*% zeta / sigma2 + inv.sigma.beta)
    mean.Beta_1 <- var.Beta_1 %*% (t(zeta) %*% (Y - eta) / sigma2)
    Beta_1 <- rnorm(1, mean.Beta_1, sqrt(var.Beta_1))
    
    if (index > nburn && (index - nburn) %% nthin == 0) {
      Beta_1.chain[, (index - nburn) / nthin] <- Beta_1
      Beta_2.chain[, (index - nburn) / nthin] <- Beta_2
      Beta_3.chain[, (index - nburn) / nthin] <- Beta_3
      sigma2.chain[(index - nburn) / nthin] <- sigma2
      acceptance_beta.chain[(index - nburn) / nthin] <- acceptance_beta
      acceptance_sigma.chain[(index - nburn) / nthin] <- acceptance_sigma
    }
  }

  list(Beta_1.chain = Beta_1.chain, Beta_2.chain = Beta_2.chain, Beta_3.chain = Beta_3.chain, sigma2.chain = sigma2.chain, 
  acceptance_beta = acceptance_beta.chain, acceptance_sigma = acceptance_sigma.chain)
}


## Read the data
gls <- load("/home/vinux/Documents/Homework/glm dat.RData")

result <- mix_sampler(X, Y, 100000, 100000, 1)
plot(result$Beta_1.chain[1,], type = "l")
plot(result$Beta_2.chain[1,], type = "l")
plot(result$Beta_3.chain[1,], type = "l")
plot(result$sigma2.chain, type = "l")

mean(result$acceptance_beta)
mean(result$acceptance_sigma)
