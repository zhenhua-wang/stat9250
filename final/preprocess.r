setwd("./final")
library(tidyverse)
library(mltools)
library(data.table)

heart <- read.csv("./data/heart_2020_cleaned.csv")
heart <- heart %>%
  select(HeartDisease, BMI, SleepTime, Smoking, Sex)
heart$HeartDisease <- ifelse(heart$HeartDisease == "Yes", 1, 0)
heart$Smoking <- ifelse(heart$Smoking == "Yes", 1, 0)
heart$Sex <- ifelse(heart$Sex == "Female", 1, 0)
print(sprintf("any NA? %s", any(is.na(heart))))

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

save(heart, train_idx, X_train, Y_train, X_test, Y_test,
  file = "./data/heart.RData")
