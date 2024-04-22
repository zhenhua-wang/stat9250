setwd('~/Workspace/Course/stat9250/final/')
load("./data/heart.RData")
library(tidyverse)
library(caret)
library(pROC)

## * methods
# helper
sigmoid <- function(Z) {
  1 / (1 + exp(-Z))
}

logloss <- function(Y, X, beta) {
  - t(Y) %*% log(sigmoid(X %*% beta)) -
    t(1 - Y) %*% log(1 - sigmoid(X %*% beta))
}

jacobian <- function(Y, X, beta) {
  -t(X) %*% (Y - sigmoid(X %*% beta)) / length(Y)
}

hessian <- function(Y, X, beta) {
  pi <- sigmoid(X %*% beta)
  W <- diag(drop(pi * (1 - pi)))
  t(X) %*% W %*% X + diag(rep(1e-8, ncol(X)))
}

# using gradient descent
gradient_descient <- function(Y, X, beta_init,
                              batch_size = 1000,
                              alpha = 0.01, lambda = 0.01,
                              epoch = 10000, eps = 1e-3) {
  beta <- beta_init
  loss <- c()
  for (t in 1:epoch) {
    # Sample minibatch
    indices <- sample(1:nrow(X), batch_size, replace = FALSE)
    X_batch <- X[indices, ]
    Y_batch <- Y[indices]
    ## fitting
    beta <- beta -
      alpha * jacobian(Y_batch, X_batch, beta) -
      lambda * 2 * beta
    loss[t] <- logloss(Y_batch, X_batch, beta) / length(indices)
    if (t > 1 && abs(loss[t] - loss[t-1]) < eps) break
  }
  return(list(loss = loss, beta = beta))
}

# using Newton-Raphson
newton_raphson <- function(Y, X, beta_init,
                           batch_size = 1000,
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
        solve(hessian(Y_batch, X_batch, beta)) %*%
        jacobian(Y_batch, X_batch, beta)
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
      batch_size = 1000,
      epoch = 10000, eps = 1e-3)
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

## * analysis using Newton-Raphson
## load data
Y <- heart$HeartDisease
X <- heart %>% dplyr::select(-HeartDisease)
X$Smoking <- as.factor(X$Smoking)
X$Sex <- as.factor(X$Sex)
X <- one_hot(as.data.table(X))
X <- as.matrix(X)
X <- cbind(rep(1, nrow(X)), X)

## split train test
shuffled_indices <- sample(nrow(X))
train_idx <- shuffled_indices[1:round(0.7 * nrow(X))]
X_train <- X[train_idx, ]
Y_train <- Y[train_idx]
X_test <- X[-train_idx, ]
Y_test <- Y[-train_idx]

## training
start.time <- Sys.time()
## result <- gradient_descient(Y_train, X_train,
##   rep(0, 7), epoch = 1000, alpha = alpha_best, batch_size = 10000)
beta_init <- rep(0, 7)#rnorm(7, 0, 0.01)
result <- newton_raphson(Y_train, X_train,
  beta_init, epoch = 100, eps = 1e-4, batch_size = 1000)
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

## bootstrap
## beta_init <- rep(0, 7)
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

## * analysis using Gradient descent
## choose parameters
shuffled_indices <- sample(nrow(X_train))
val_idx <- shuffled_indices[1:round(0.3 * nrow(X_train))]
X_val <- X_train[val_idx, ]
Y_val <- Y_train[val_idx]
X_valtrain <- X_train[-val_idx, ]
Y_valtrain <- Y_train[-val_idx]
alpha_list <- c(0.01, 0.03, 0.05, 0.08, 0.01, 0.05, 0.1)
spec_list <- c()
for (alpha in alpha_list) {
  val_result <-
    gradient_descient(Y_valtrain, X_valtrain, rep(0, 7),
      epoch = 1000, alpha = alpha, batch_size = 1000)
  optimal_threshold <- get_optimal_threshold(Y_val, X_val, val_result$beta)
  Y_pred <- predict_label(X_test, val_result$beta, optimal_threshold)
  conf_mat <- confusionMatrix(factor(Y_pred), factor(Y_test))
  spec_list <- c(spec_list, conf_mat$byClass["Specificity"])
  cat(alpha, "\r")
}
alpha_best <- alpha_list[which.max(spec_list)]
print(paste0("best alpha is ", alpha_best))