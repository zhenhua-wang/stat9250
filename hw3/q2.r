library(tidyverse)

lasso_loss <- function(X, Y, beta) {
  t(Y - X%*%beta) %*% (Y - X%*%beta) + lambda * sum(abs(beta))
}

lasso_regression <- function(X, Y, lambda_init, epoch = 10000, eps = 1e-4) {
  lambda <- lambda_init
  loss <- c()
  ## scale X and Y
  Y <- scale(Y, center = TRUE, scale = TRUE)
  X <- scale(X, center = TRUE, scale = TRUE)
  X <- cbind(rep(1, length(Y)), X)
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
    loss[i] <- lasso_loss(X, Y, beta)
    ## quit when converge
    if ((length(loss) > 1) && (abs(loss[i] - loss[i - 1]) < eps)) break
  }
  return(list(beta = beta, loss = loss))
}

peru <- read.table("./hw3/peru.txt", header = TRUE)
Y <- peru$Systol
X <- as.matrix(peru %>% select(-Systol))
N <- nrow(X)

## validation
train_idx <- 1:floor(N * 0.7)
X_train <- X[train_idx, ]
Y_train <- Y[train_idx]
X_test <- X[-train_idx, ]
Y_test <- Y[-train_idx]
lambda_list <- c(0.01, 0.05, 0.1, 0.3, 0.7, 0.9)
test_loss <- c()
for (lambda in lambda_list) {
  result <- lasso_regression(X_train, Y_train, lambda = lambda, epoch = 10000)
  test_loss <- c(test_loss,
    lasso_loss(cbind(rep(1, length(Y_test)), X_test), Y_test, result$beta))
}
best_idx <- which(test_loss == min(test_loss))

## result
lambda <- lambda_list[best_idx]
result <- lasso_regression(X, Y, lambda = lambda, epoch = 10000)
loss <- result$loss
result$beta
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
lasso_loss(cbind(rep(1, N), X), Y, beta)
