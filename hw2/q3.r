## * Compute neighboring matrix
M1_base <- cbind(diag(3), rep(0, 3)) + cbind(rep(0, 3), -diag(3))
M1 <- rbind(
  cbind(M1_base, matrix(0, 3, 12)),
  cbind(matrix(0, 3, 4), M1_base, matrix(0, 3, 8)),
  cbind(matrix(0, 3, 8), M1_base, matrix(0, 3, 4)),
  cbind(matrix(0, 3, 12), M1_base))
M2_base <- cbind(diag(4), -diag(4))
M2 <- rbind(
  cbind(M2_base, matrix(0, 4, 8)),
  cbind(matrix(0, 4, 4), M2_base, matrix(0, 4, 4)),
  cbind(matrix(0, 4, 8), M2_base))
M <- rbind(M1, M2)

D <- function(V) {
  X <- matrix(V, 4, 4)
  sum(abs(diff(X, 1, 1))) + sum(abs(diff(t(X), 1, 1)))
}

## test
V1 <- c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
print(sprintf("For example 1, D(X) = %i", sum(abs(M %*% V1))))
V2 <- c(1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)
print(sprintf("For example 2, D(X) = %i", sum(abs(M %*% V2))))

## * Metropolis-Hastings algorithm
random_walk_metropolis <- function(true_dens, prop_func,
                                   X_init, sample_size, lambda) {
  X <- X_init
  samples <- matrix(NA, sample_size, length(X_init))
  for (i in 1:sample_size) {
    X_new <- prop_func(X)
    alpha <- min(1, true_dens(X_new, lambda) / true_dens(X, lambda))
    U <- runif(1)
    if (U < alpha) {
      samples[i, ] <- X_new
      X <- X_new
    } else {
      samples[i, ] <- X
    }
  }
  return(samples)
}

true_dens <- function(V, lambda) {
  N_X <- sum(V)
  D_X <- D(V) #sum(abs(M %*% V))
  0.8^N_X * 0.2^(16 - N_X) * exp(-lambda * D_X)
}

diag_all_one <- function(V) {
  X <- matrix(V, 4, 4)
  all(diag(X) == 1)
}

## * MH1
prop_func1 <- function(X) {
  rbinom(16, 1, 0.5)
}

V <- matrix(0, 4, 4)
sample_size <- 10000
samples <- random_walk_metropolis(true_dens, prop_func1, V, sample_size, 0.5)
plot(1:sample_size, apply(samples, 1, D),
  type = "l", col = "blue", xlab = "iteration", ylab = "N(X)"
)
print(sprintf(
  "Probability for algo 1, lambda = 0.5: %.3f",
  mean(apply(samples, 1, diag_all_one))
))

samples <- random_walk_metropolis(true_dens, prop_func1, V, sample_size, 1)
plot(1:sample_size, apply(samples, 1, D),
  type = "l", col = "blue", xlab = "iteration", ylab = "N(X)"
)
print(sprintf(
  "Probability for algo 1, lambda = 1: %.3f",
  mean(apply(samples, 1, diag_all_one))
))

## * MH2
prop_func2 <- function(X) {
  idx <- sample(16, 1, replace = TRUE, prob = rep(1 / 16, 16))
  X[idx] <- rbinom(1, 1, 0.5)
  return(X)
}

samples <- random_walk_metropolis(true_dens, prop_func2, V, sample_size, 0.5)
plot(1:sample_size, apply(samples, 1, D),
  type = "l", col = "blue", xlab = "iteration", ylab = "N(X)"
)
print(sprintf(
  "Probability for algo 2, lambda = 0.5: %.3f",
  mean(apply(samples, 1, diag_all_one))
))

samples <- random_walk_metropolis(true_dens, prop_func2, V, sample_size, 1)
plot(1:sample_size, apply(samples, 1, D),
  type = "l", col = "blue", xlab = "iteration", ylab = "N(X)"
)
print(sprintf(
  "Probability for algo 2, lambda = 1: %.3f",
  mean(apply(samples, 1, diag_all_one))
))

## * MH3
sample_size <- 100000
lambda <- 1
V <- V2
X <- matrix(V, 4, 4)
samples <- matrix(NA, sample_size, length(V))
for (i in 1:sample_size) {
  ## diagonal
  X.star <- X
  diag(X.star) <- rbinom(4, 1, 0.5)
  alpha <- min(1, true_dens(c(X.star), lambda) / true_dens(c(X), lambda))
  U <- runif(1)
  if (U < alpha) {
    samples[i, ] <- c(X.star)
    X <- X.star
  } else {
    samples[i, ] <- c(X)
  }
  ## off-diagonal
  X.star <- X
  X.star[col(X.star) != row(X.star)] <- rbinom(12, 1, 0.5)
  alpha <- min(1, true_dens(c(X.star), lambda) / true_dens(c(X), lambda))
  U <- runif(1)
  if (U < alpha) {
    samples[i, ] <- c(X.star)
    X <- X.star
  } else {
    samples[i, ] <- c(X)
  }
}

plot(1:sample_size, apply(samples, 1, D),
  type = "l", col = "blue", xlab = "iteration", ylab = "N(X)"
)

print(sprintf(
  "Probability for algo 2, lambda = 1: %.3f",
  mean(apply(samples, 1, diag_all_one))
))
