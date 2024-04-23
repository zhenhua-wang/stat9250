setwd('~/Workspace/Course/stat9250/final/')
library(tidyverse)
library(caret)
library(pROC)

## * methods
# helper
standardize <- function(matrix_data, col_min, col_max) {
  # Standardize each column
  norm_matrix <- sweep(matrix_data, 2, col_min, FUN = "-")
  norm_matrix <- sweep(norm_matrix, 2, col_max - col_min, FUN = "/")
  return(norm_matrix)
}

sigmoid <- function(Z) {
  1 / (1 + exp(-Z))
}

logloss <- function(Y, X, beta) {
  - t(Y) %*% log(sigmoid(X %*% beta)) -
    t(1 - Y) %*% log(1 - sigmoid(X %*% beta))
}

jacobian <- function(Y, X, beta, lambda) {
  -t(X) %*% (Y - sigmoid(X %*% beta)) + 2 * lambda * beta
}

hessian <- function(Y, X, beta, lambda) {
  pi <- sigmoid(X %*% beta)
  W <- diag(drop(pi * (1 - pi)))
  t(X) %*% W %*% X + diag(rep(2 * lambda, ncol(X)))
}

# using Newton-Raphson
newton_raphson <- function(Y, X, beta_init,
                           batch_size = 1000,
                           lambda = 0.01,
                           epoch = 10000,
                           eps = 1e-4) {
  beta <- beta_init
  loss <- c()
  for (t in 1:epoch) {
    # Sample minibatch
    indices <- sample(1:nrow(X), nrow(X), replace = FALSE)
    num_batches <- floor(nrow(X) / batch_size)
    for (i in 1:num_batches) {
      batch_idx <- indices[((i-1)*batch_size + 1):(i * batch_size)]
      X_batch <- X[batch_idx, ]
      Y_batch <- Y[batch_idx]
      ## fitting
      beta <- beta -
        solve(hessian(Y_batch, X_batch, beta, lambda)) %*%
        jacobian(Y_batch, X_batch, beta, lambda)
    }
    loss[t] <- logloss(Y_batch, X_batch, beta) / batch_size
    if (t > 1 && abs(loss[t] - loss[t-1]) < eps) break
  }
  return(list(loss = loss, beta = beta))
}

# predict
predict_prob <- function(X, beta) {
  c(sigmoid(X %*% beta))
}

predict_label <- function(X, beta, threshold = 0.5) {
  1 * (predict_prob(X, beta) >= threshold)
}

# evaluation
get_optimal_threshold <- function(Y_test, X_test, beta) {
  Y_test_prob <- predict_prob(X_test, beta)
  roc_result <- roc(Y_test, Y_test_prob)
  optimal_idx <-
    which.max(roc_result$sensitivities + roc_result$specificities - 1)
  return(roc_result$thresholds[optimal_idx])
}

accuracy <- function(Y, Y_pred) {
  mean(Y == Y_pred)
}

## bootstrap
bootstrap <- function(optim_func, num_boot,
                      Y, X, beta_init,
                      batch_size = 1000,
                      alpha = 0.01, lambda = 0.01,
                      epoch = 10000, eps = 1e-3) {
  model_list <- list()
  n <- nrow(X)
  for (i in 1:num_boot) {
    sample_indices <- sample(1:n, n, replace = TRUE)  # Sampling with replacement
    X_boot <- X[sample_indices, ]
    Y_boot <- Y[sample_indices]
    model <- optim_func(Y_boot, X_boot, beta_init,
      batch_size = batch_size,
      lambda = lambda,
      epoch = epoch, eps = eps)
    model_list[[i]] <- model
    cat(i, "\r")
  }
  return(model_list)
}

bootstrap_CI <- function(result_boot) {
  num_param <- length(result_boot[[1]]$beta)
  num_boot <- length(result_boot)
  params <- matrix(NA, num_boot, num_param)
  for (i in 1:num_boot) {
    params[i, ] <- result_boot[[i]]$beta
  }
  return(list(beta_mean = apply(params, 2, mean),
    beta_CI = apply(params, 2, quantile, c(0.025, 0.975))))
}

bootstrap_predict <- function(Y, X, result_boot) {
  num_obs <- length(Y)
  num_boot <- length(result_boot)
  Y_pred_boot <- matrix(NA, num_boot, num_obs)
  for (i in 1:num_boot) {
    Y_pred_boot[i, ] <- predict_prob(X, result_boot[[i]]$beta)
  }
  Y_pred_prob <- apply(Y_pred_boot, 2, mean)
  roc_result <- roc(Y, Y_pred_prob)
  optimal_idx <-
    which.max(roc_result$sensitivities + roc_result$specificities - 1)
  threshold_best <- roc_result$thresholds[optimal_idx]
  Y_pred <- ifelse(Y_pred_prob >= threshold_best, 1, 0)
  return(Y_pred)
}

## * load data
load("./data/heart.RData")
X_min_cont <- apply(X_train[, 2:3], 2, min)
X_max_cont <- apply(X_train[, 2:3], 2, max)
X_train[, 2:3] <- standardize(X_train[, 2:3], X_min_cont, X_max_cont)
X_test[, 2:3] <- standardize(X_test[, 2:3], X_min_cont, X_max_cont)

## * analysis using Newton-Raphson
## training
start.time <- Sys.time()
beta_init <- rep(0, ncol(X_train))
result <- newton_raphson(Y_train, X_train,
  beta_init, epoch = 100, eps = 1e-4, batch_size = 1000, lambda = 0.05)
end.time <- Sys.time()
print(end.time - start.time)
plot(result$loss, type = "l")

## result
optimal_threshold <- get_optimal_threshold(Y_test, X_test, result$beta)
Y_pred <- predict_label(X_test, result$beta, optimal_threshold)
conf_mat <- confusionMatrix(factor(Y_pred), factor(Y_test))
conf_mat$byClass["Specificity"]
accuracy(Y_test, Y_pred)
result$beta

result_roc <- roc(Y_test, predict_prob(X_test, result$beta))
plot(result_roc)
auc(result_roc)

## bootstrap
## beta_init <- rep(0, ncol(X_train))
## result_boot <- bootstrap(newton_raphson, 100, Y_train, X_train,
##   beta_init, epoch = 100, eps = 1e-4, batch_size = 1000)
## save(result_boot, file = "./data/newton.RData")

load("./data/newton.RData")
beta_boot <- bootstrap_CI(result_boot)
beta_boot$beta_mean

Y_boot <- bootstrap_predict(Y_test, X_test, result_boot)
conf_mat <- confusionMatrix(factor(Y_boot), factor(Y_test))
conf_mat$byClass["Specificity"]
accuracy(Y_test, Y_boot)
beta_boot$beta_mean

## roc
plot(roc(Y_test, predict_prob(X_test, beta_boot$beta_mean)))

## * parameter tuning
shuffled_indices <- sample(nrow(X_train))
val_idx <- shuffled_indices[1:round(0.3 * nrow(X_train))]
X_val <- X_train[val_idx, ]
Y_val <- Y_train[val_idx]
X_valtrain <- X_train[-val_idx, ]
Y_valtrain <- Y_train[-val_idx]
alpha_list <- c(0.01, 0.03, 0.05, 0.08, 0.1, 0.3, 0.5, 0.8)
acc_list <- c()
for (alpha in alpha_list) {
  val_result <-
    newton_raphson(Y_valtrain, X_valtrain, rep(0, ncol(X_train)),
      epoch = 100, lambda = alpha, batch_size = 1000)
  optimal_threshold <- get_optimal_threshold(Y_val, X_val, val_result$beta)
  Y_pred <- predict_label(X_test, val_result$beta, optimal_threshold)
  acc_list <- c(acc_list, accuracy(Y_test, Y_pred))
  cat(alpha, "\r")
}
alpha_best <- alpha_list[which.max(acc_list)]
print(paste0("best alpha is ", alpha_best))
