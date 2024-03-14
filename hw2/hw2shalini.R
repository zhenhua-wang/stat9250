
# Problem 3 ---------------------------------------------------------------
library(batchmeans)

X <- matrix(c(1,0,0,0,
              1,1,0,0,
              0,0,0,0,
              0,0,0,0), ncol = 4, nrow = 4)

M <- matrix(0, nrow = 24, ncol = 16)
i <- 1
for (j in 1:16) {
  M[i,j] <- 1 
  M[i,j+1] <- -1 
  M[i+1,j] <- 1 
  M[i+1,j+4] <- -1 
  i <- i + 2
}

lambda <- 0.5

prob_dist <- function(lambda, X){
  N_X <- sum(X)
  D_X <- sum(abs(M%*%as.vector(t(X))))
  prob <- 0.8^N_X * 0.2^(16-N_X) * exp(-lambda*D_X)
  return(prob)
}

Q1 <- function(X){
  X_star <- matrix(rbinom(16,1,0.5), nrow = 4, ncol = 4)
  return(X_star)
}
Q2 <- function(X_new){
  loc <- sample(16, 1, replace = TRUE, prob = rep(0.0625, 16))
  X_new[loc] <- rbinom(1,1,0.5)
  return(X_new)
}

#MCMC

N <- round(exp(seq(2, 13, length.out=100)))

MCMC <- function(Q, lambda){
  diag_prob <- array(NA, length(N))
  diag_se <- array(NA, length(N))
  
  for (j in 1:length(N)) {
    diag_sum <- array(0, N[j])
    for (i in 1:N[j]) {
      X_new <- Q(X)
      ratio <- prob_dist(lambda, X_new)/prob_dist(lambda, X)
      alpha <- min(1,ratio)
      if(runif(1) < alpha){
        X <- X_new
      }
      if(sum(diag(X)) == 4)
        diag_sum[i] <- 1
    }
    #diag_prob[j] <- sum(diag_sum)/N[j]
    diag_prob[j] <- bm(diag_sum)[1]
    diag_se[j] <- bm(diag_sum)[2]
  }
  return(list(diag_prob = diag_prob, diag_se = diag_se))
}

par(mfrow = c(2,4))
for (l in c(0.5,1)) {
  q1 <- MCMC(Q1,l)
  q2 <- MCMC(Q2,l)
  plot(N,q1$diag_prob, type = "l", main = bquote(paste("Q1 and ",lambda == .(l))), ylab = "Probability")
  plot(N,q1$diag_se, type = 'l', main = bquote(paste("Q1 and ",lambda == .(l))), ylab = "se")
  plot(N,q2$diag_prob, type = "l", main = bquote(paste("Q2 and ",lambda == .(l))), ylab = "Probability")
  plot(N,q2$diag_se, type = 'l', main = bquote(paste("Q2 and ",lambda == .(l))), ylab = "se")
}

