library(tidyverse)

lasso_regression <- function(X, Y, lambda_init, epoch = 10000, eps = 1e-4) {
  lambda <- lambda_init
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
        beta[j] <- delta_linear + lambda
      } else if (delta_linear > lambda) {
        beta[j] <- delta_linear - lambda
      } else {
        beta[j] <- 0
      }
      ## restore full residual
      resid <- resid - X[, j] * beta[j]
    }
    ## record loss
    loss[i] <- t(Y - X%*%beta) %*% (Y - X%*%beta) + lambda * sum(abs(beta))
    ## quit when converge
    if ((length(loss) > 1) && (abs(loss[i] - loss[i - 1]) < eps)) break
  }
  return(list(beta = beta, loss = loss))
}

peru <- read.table("./hw3/peru.txt", header = TRUE)
Y <- peru$Systol
X <- as.matrix(peru %>% select(-Systol))
N <- nrow(X)
## initialize
loss <- c()
epoch <- 10000
lambda <- 0.1
result <- lasso_regression(X, Y, lambda, epoch)
loss <- result$loss
result$beta

## result
plot(loss, type = "l")
beta <- matrix(result$beta, dim(X)[2] + 1, 1)
rownames(beta) <- c("(intercept)", colnames(X))
beta

## compare with glmnet
library(glmnet)
Y <- scale(Y, center = TRUE, scale = TRUE)
X <- scale(X, center = TRUE, scale = TRUE)
model <- glmnet(X, Y, alpha = 1, lambda = lambda)
coef(model)
beta <- Matrix(coef(model), sparse = FALSE)
