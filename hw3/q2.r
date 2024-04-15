library(tidyverse)

peru <- read.table("./hw3/peru.txt", header = TRUE)
Y <- peru$Systol
X <- as.matrix(peru %>% select(-Systol))
N <- nrow(X)
## initialize
loss <- c()
epoch <- 10000
eps <- 1e-4
lambda <- 0.1
## scale X and Y
Y <- scale(Y, center = TRUE, scale = TRUE)
X <- scale(X, center = TRUE, scale = TRUE)
X <- cbind(rep(1, N), X)
P <- ncol(X)
beta <- rep(0, P)
for (i in 1:epoch) {
  resid <- Y - X %*% beta
  ## pick a coordinate
  for (j in 1:P) {
    ## get residual for all except for j
    resid <- resid + X[, j] * beta[j]
    delta_linear <- mean(resid * X[, j])
    ## proceed with the negative direction
    if (delta_linear < -lambda) {
      beta[j] <- delta_linear + lambda / 2
    } else if (delta_linear > lambda) {
      beta[j] <- delta_linear - lambda / 2
    }
    ## restore full residual
    resid <- resid - X[, j] * beta[j]
  }
  ## record loss
  loss[i] <- t(Y - X%*%beta) %*% (Y - X%*%beta) + lambda * sum(abs(beta))
  ## quit when converge
  if ((length(loss) > 1) && (abs(loss[i] - loss[i - 1]) < eps)) break
}

plot(loss, type = "l")
beta <- matrix(beta, P, 1)
rownames(beta) <- colnames(X)
rownames(beta)[1] <- "(intercept)"
beta

library(glmnet)
model <- glmnet(X[, 2:10], Y, alpha = 1, lambda = lambda)
coef(model)
