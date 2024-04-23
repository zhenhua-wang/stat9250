library(ggplot2)
library(gridExtra)
library(MASS)
library(truncnorm)
library(coda)

load("heart.RData")

str(heart)

pdf(paste0("Descriptive.pdf"),height=4,width=6)
g1 <- ggplot(heart, aes(x = as.factor(HeartDisease), y = BMI, fill = as.factor(HeartDisease))) +
  geom_boxplot() +
  labs(x = "Heart Disease", y = "BMI") +
  theme_minimal() +
  guides(fill = FALSE)
g2 <- ggplot(heart, aes(x = as.factor(HeartDisease), y = SleepTime, fill = as.factor(HeartDisease))) +
  geom_boxplot() +
  labs(x = "Heart Disease", y = "SleepTime") +
  theme_minimal() +
  guides(fill = FALSE)
g3 <- ggplot(heart, aes(x = factor(Smoking, labels = c("No", "Yes")), fill = factor(HeartDisease))) +
  geom_bar(position = "dodge") +
  labs(x = "Smoking", y = "Count", fill = "Heart Disease") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), 
                    labels = c("0" = "No", "1" = "Yes")) +
  theme_minimal()
g4 <- ggplot(heart, aes(x = factor(Sex, labels = c("Male", "Female")), fill = factor(HeartDisease))) +
  geom_bar(position = "dodge") +
  labs(x = "Sex", y = "Count", fill = "Heart Disease") +
  scale_fill_manual(values = c("0" = "blue", "1" = "red"), 
                    labels = c("0" = "No", "1" = "Yes")) +
  theme_minimal()
grid.arrange(g1, g2, g3, g4, ncol = 2)
dev.off()



# t link model------------------------------------------------------------------

set.seed(1919)
n <- 100 #number of iterations
X <- unname(data.matrix(X_train))
X <- unname(as.matrix(heart[train_idx,-1]))
X <- cbind(1,X)
y <- Y_train
N <- nrow(X)
#A <- solve(t(X) %*% X)
beta_mc <- matrix(NA,ncol(X),n)
lambda_mc <- matrix(NA,N,n)
#Z_mc <- matrix(NA,N,n)
#lambda_mc <- matrix(NA,N,n)
Z <- matrix(NA,N,1)

df <- 8

myprobit <- glm(HeartDisease ~ BMI + SleepTime  + as.factor(Smoking) + Sex, family = binomial(link = "probit"), 
                data = heart)
#myprobit <- glm(Y_train ~ X_train, family = binomial(link = "probit"))

summary(myprobit)
beta <- as.vector(myprobit$coefficients)
lambda <- rgamma(N,df/2,2/df)

for (i in 1:n) {
  for (j in 1:N) {
    m <- t(X[j,]) %*% beta
    if(y[j] == 1){
      Z[j] <- rtruncnorm(1,a = 0, mean = m, sd = sqrt(1/lambda[j]))
    }else{
      Z[j] <- rtruncnorm(1,b = 0, mean = m, sd = sqrt(1/lambda[j]))
    }
  }
  X_new <- X * lambda
  K <- solve(t(X) %*% X_new)
  Z_new <-  Z * lambda
  beta_hat <- K %*% t(X) %*% Z_new
  beta_mc[,i] <- beta <- mvrnorm(1,mu = beta_hat, Sigma = K)
  for (k in 1:N) {
    sig <- df + (Z[k]-t(X[k,]) %*% beta)^2
    lambda_mc[k,i] <- lambda[k] <- rgamma(1,(df+1)/2,2/sig)
  }
  if(i%%100==0){
    print(i)
  }
}

mcmc <- list(beta_mc = beta_mc, lambda_mc = lambda_mc)
save(mcmc,file = paste0("mcmc.Rda"))

rownames(beta_mc) <- c("Beta 0", "Beta 1", "Beta 2", "Beta 3", "Beta 4")
beta_mc_draws <- mcmc(t(beta_mc))
pdf(paste0("beta_mc_ACF.pdf"),width=8, height=8)
autocorr.plot(beta_mc_draws)
dev.off()

load("beta_mc_2000.Rda")
n <- 2000

pdf(paste0("beta_mc_trace_new.pdf"),width=8, height=8)
par(mfrow = c(3,2))
for(i in 1:5){
  k <- i-1
  plot(500:n,beta_mc[i,500:n], xlab = "", ylab = "", main = paste("Beta", k), type = 'l')
}
dev.off()

pdf(paste0("lambda_mc_trace.pdf"),width=8, height=8)
par(mfrow = c(3,2))
for(i in 1:1){
  k <- i-1
  plot(1:n,lambda_mc[i,], xlab = "", ylab = "", main = paste("Beta", k))
}
dev.off()


# Probit model ------------------------------------------------------------


set.seed(1919)
n <- 2000 #number of iterations
#X <- unname(data.matrix(X_train))
X <- unname(as.matrix(heart[train_idx,-1]))
X <- cbind(1,X)
y <- Y_train
N <- nrow(X)
beta_mc <- matrix(NA,ncol(X),n)
lambda_mc <- matrix(NA,N,n)
Z <- matrix(NA,N,1)

myprobit <- glm(HeartDisease ~ BMI + SleepTime  + Smoking + Sex, family = binomial(link = "probit"), 
                data = heart)

summary(myprobit)
beta <- as.vector(myprobit$coefficients)

A <- solve(t(X) %*% X)
#Gibbs sampling
for (i in 1:n) {
  for (j in 1:N) {
    m <- t(X[j,]) %*% beta
    if(y[j] == 1){
      Z[j] <- rtruncnorm(1,a = 0, mean = m, sd = 1)
    }else{
      Z[j] <- rtruncnorm(1,b = 0, mean = m, sd = 1)
    }
  }
  B <- t(X) %*% Z
  mvn_mean <- A %*% B
  beta_mc[,i] <- beta <- mvrnorm(1, mvn_mean, A)
  if(i%%100==0){
    print(i)
  }
}

save(beta_mc,file = paste0("beta_mc_probit.Rda"))

rownames(beta_mc) <- c("Beta 0", "Beta 1", "Beta 2", "Beta 3", "Beta 4")
beta_mc_draws <- mcmc(t(beta_mc))
pdf(paste0("beta_mc_ACF.pdf"),width=8, height=8)
autocorr.plot(beta_mc_draws)
dev.off()

load("beta_mc_probit.Rda")

pdf(paste0("beta_mc_trace.pdf"),width=8, height=8)
par(mfrow = c(3,2))
for(i in 1:5){
  k <- i-1
  plot(1:n,beta_mc[i,], xlab = "", ylab = "", main = paste("Beta", k), type = "l")
}
dev.off()


# model evaluation --------------------------------------------------------

load("beta_mc_probit.Rda")

CI_lower <- apply(beta_mc, 1, function(x) quantile(x,0.075))
CI_upper <- apply(beta_mc, 1, function(x) quantile(x,0.975))
param_est <- apply(beta_mc, 1, mean)

library(caret)  # for confusion matrix
library(pROC)  # for AUC

predicted_probabilities <- pnorm(X %*% t(t(param_est)))
summary(predicted_probabilities)
roc_obj <- roc(y, predicted_probabilities)

pdf(paste0("roc_probit.pdf"),width=8, height=8)
plot(roc_obj, main="ROC Curve", col="#1c61b6")
dev.off()

# Calculate the optimal threshold by Youden's Index
# Youden's Index
optimal_idx <- which.max(roc_obj$sensitivities + roc_obj$specificities - 1)
optimal_threshold <- roc_obj$thresholds[optimal_idx]
cat("Optimal threshold by Youden's Index is:", optimal_threshold, "\n")

X_test <- cbind(1,heart[-c(train_idx),-1])
predicted_prob_test <- pnorm( as.matrix(X_test) %*% t(t(param_est)))
Yhat <- ifelse(predicted_prob_test > optimal_threshold, 1, 0)
Y_test <- heart[-c(train_idx),1]
table <- table(Y_test, Yhat)

conf_mat <- confusionMatrix(factor(Yhat), factor(Y_test))

print(conf_mat$byClass)  # Precision, Recall, F1, Sensitivity, Specificity
print(conf_mat$overall['Accuracy'])  # Accuracy
print(paste("AUC:", pROC::auc(pROC::roc(Y_test, predicted_prob_test))))  # AUC

