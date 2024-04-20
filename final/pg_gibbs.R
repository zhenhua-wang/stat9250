load("/home/kinux/Documents/9250/final/heart.RData")

library(MASS)  # for mvtnorm functions
library(Matrix)  # for sparse matrix operations
library(BayesLogit)  # Polya-Gamma sampling
library(caret)  # for confusion matrix
library(pROC)  # for AUC

# Splite with training and test
set.seed(123)
train_index <- sample(1:nrow(heart), nrow(heart) * 0.7)
train <- heart[train_index,]
test <- heart[-train_index,]

# Gibbs sampling
X <- model.matrix(~., data=train[,2:3])
Y <- train$HeartDisease
Z_sex <- model.matrix(~ factor(Sex) - 1, data=train)
Z_smoke <- model.matrix(~ factor(Smoking) - 1, data=train)

GibbsSamplerLogit <- function(X, Y, Z_u, Z_v, nburn, nsim, nthin, a = 0.1, b = 0.1) {
  N <- length(Y)
  p <- dim(X)[2]
  n_u <- ncol(Z_u)
  n_v <- ncol(Z_v)
  
  
  # Initialize parameters
  Mu <- rep(0, N)  # Initialize theta with zeros
  Beta <- c(-5,0)  # Initialize beta as a column vector
  Sigma2_U <- 1
  Sigma2_V <- 1
  U <- rep(1, n_u)
  V <- rep(1, n_v)
  
  inv_sigma_u <- Diagonal(n_u)
  inv_sigma_v <- Diagonal(n_v)
  
  # Storage for parameters
  Beta.chain <- array(0, dim = c(p, nsim / nthin))
  Sigma2_u.chain <- rep(0, nsim / nthin)
  Sigma2_v.chain <- rep(0, nsim / nthin)
  U.chain <- array(0, dim = c(n_u, nsim / nthin))
  V.chain <- array(0, dim = c(n_v, nsim / nthin))
  Mu.chain <- array(0, dim = c(N, nsim / nthin))
  
  
  for (index in 1:(nsim + nburn)) {
    if (index %% 10000 == 0) cat(index, "\n")
    # Update latent variable w
    w <- rpg(N, 1, Mu)
    Omega <- Diagonal(x=w)
    UO <- t(Z_u) %*% Omega
    OU <- UO %*% Z_u
    VO <- t(Z_v) %*% Omega
    OV <- VO %*% Z_v
    
    # Update Sigma_2_u
    a_prime_u <- a + n_u / 2
    b_prime_u <- b + 0.5 * t(U) %*% U
    Sigma2_u <- 1 / rgamma(1, shape = a_prime_u, rate = b_prime_u) 
    
    # Update Sigma_2_v
    a_prime_v <- a + n_v / 2
    b_prime_v <- b + 0.5 * t(V) %*% V
    Sigma2_v <- 1 / rgamma(1, shape = a_prime_v, rate = b_prime_v)
    
    # Update Beta
    var.Beta <- solve(t(X) %*% Omega %*% X + 0.01 * Diagonal(p))
    mean.Beta <- var.Beta %*% (t(X) %*% Omega %*% ((Y - 1/2)/w - Z_u %*% U - Z_v %*% V))
    Beta <- as.vector(mvrnorm(1, mu = mean.Beta, Sigma = var.Beta))
    
    # Update U
    var.U <- solve(OU  + inv_sigma_u/Sigma2_u)
    mean.U <- var.U %*% (UO %*%((Y - 1/2)/w - X %*% Beta - Z_v %*% V))
    #U <- as.vector(rmvn(1, as.vector(mean.U), var.U))
    U <- as.vector(mvrnorm(1, mu = mean.U, Sigma = var.U))
    
    # Update V
    var.V <- solve(OV  + inv_sigma_v/Sigma2_v)
    mean.V <- var.V %*% (VO %*%((Y - 1/2)/w - X %*% Beta - Z_u %*% U))
    #V <- as.vector(rmvn(1, as.vector(mean.V), var.V))
    V <- as.vector(mvrnorm(1, mu = mean.V, Sigma = var.V))
    
    # Update mu
    Mu <- as.vector(X %*% Beta + Z_u %*% U + Z_v %*% V)
    
    if (index > nburn && (index - nburn) %% nthin == 0) {
      Beta.chain[, (index - nburn) / nthin] <- Beta
      U.chain[, (index - nburn) / nthin] <- U
      V.chain[, (index - nburn) / nthin] <- V
      Sigma2_u.chain[(index - nburn) / nthin] <- Sigma2_u
      Sigma2_v.chain[(index - nburn) / nthin] <- Sigma2_v
      Mu.chain[,(index-nburn)/nthin] <- Mu
    }
  }
  
  list(Beta.chain = Beta.chain, U.chain = U.chain, Sigma2_u.chain = Sigma2_u.chain, Mu.chain = Mu.chain, V.chain = V.chain, Sigma2_v.chain = Sigma2_v.chain)
}

nsim <- 1000
nburn <- 1000
nthin <- 1
t <- system.time(example_results <- GibbsSamplerLogit(X, Y, Z_sex, Z_smoke, nburn, nsim, nthin))
plot(example_results$Beta.chain[1,], type="l", xlab="Iteration", ylab="Beta1")
plot(example_results$Beta.chain[2,], type="l", xlab="Iteration", ylab="Beta2")
plot(example_results$Mu.chain[1,], type="l", xlab="Iteration", ylab="Mu1")
plot(example_results$Mu.chain[2,], type="l", xlab="Iteration", ylab="Mu2")


# Compute the predicted probabilities
predicted_probabilities <- plogis(rowMeans(example_results$Mu.chain))
# Compute the ROC curve
roc_obj <- roc(Y, predicted_probabilities)

# Plot ROC curve
plot(roc_obj, main="ROC Curve", col="#1c61b6")
abline(a=0, b=1, col="red")  # Adding reference line
# Youden's Index
optimal_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
optimal_threshold <- roc_obj$thresholds[optimal_idx]

cat("Optimal threshold by Youden's Index is:", optimal_threshold, "\n")

# predict the probability of heart disease in test data
X_test <- model.matrix(~., data=test[,2:3])
Beta <- rowMeans(example_results$Beta.chain)
U <- rowMeans(example_results$U.chain)
V <- rowMeans(example_results$V.chain)
Z_u_test <- model.matrix(~ factor(test$Sex) - 1, data=test)
Z_v_test <- model.matrix(~ factor(test$Smoking) - 1, data=test)
Mu_test <- X_test %*% Beta + Z_u_test %*% U + Z_v_test %*% V
predicted_prob <- plogis(Mu_test)
p <- plogis(Mu_test)
Yhat <- ifelse(p > optimal_threshold, 1, 0)
Y_test <- test$HeartDisease
table <- table(Y_test, Yhat)
rowMeans(example_results$U.chain)
rowMeans(example_results$V.chain)
# Calculating FNR
FN <- table[1,2] # Predicted 0 while actual was 1
TP <- table[2,2] # Predicted 1 while actual was 1

# Calculating FNR
FNR <- FN / (FN + TP)
print(FNR)

# Calculating FPR
FP <- table[2,1] # Predicted 1 while actual was 0 
TN <- table[1,1] # Predicted 0 while actual was 0

FPR <- FP / (FP + TN)
print(FPR)

conf_mat <- confusionMatrix(factor(Yhat), factor(Y_test))

print(conf_mat$byClass)  # Precision, Recall, F1, Sensitivity, Specificity
print(conf_mat$overall['Accuracy'])  # Accuracy
print(paste("AUC:", pROC::auc(pROC::roc(Y_test, predicted_prob))))  # AUC

## Compare the result with GLM method
logit_model <- glm(HeartDisease ~ BMI + SleepTime + (1|Sex) + (1|Smoking), data=train, family=binomial())
logit_prob <- predict(logit_model, test, type="response")
logit_Yhat <- ifelse(logit_prob > optimal_threshold, 1, 0)

conf_mat_logit <- confusionMatrix(factor(logit_Yhat), factor(Y_test))

print(conf_mat_logit$byClass)  # Precision, Recall, F1, Sensitivity, Specificity
print(conf_mat_logit$overall['Accuracy'])  # Accuracy
print(paste("AUC:", pROC::auc(pROC::roc(Y_test, logit_prob))))  # AUC
