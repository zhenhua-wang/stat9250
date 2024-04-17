setwd("./final")
library(tidyverse)

heart <- read.csv("./data/heart_2020_cleaned.csv")
heart <- heart %>%
  select(HeartDisease, BMI, SleepTime, Smoking, Sex)
heart$HeartDisease <- ifelse(heart$HeartDisease == "Yes", 1, 0)
heart$Smoking <- ifelse(heart$Smoking == "Yes", 1, 0)
heart$Sex <- ifelse(heart$Sex == "Female", 1, 0)

save(heart, file = "./data/heart.RData")
print(sprintf("any NA? %s", any(is.na(heart))))
