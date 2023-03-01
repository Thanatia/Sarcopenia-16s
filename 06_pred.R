library(caret)
library(glmnet)
library(pROC)
set.seed(1024)
data_asv <- ASV_counts_pred
diagnosis_col <- "Diag"

lasso_log_reg <- function(train_data, test_data) {
  x_train <- model.matrix(Diag ~ ., data = train_data)[,-ncol(train_data)]
  x_test <- model.matrix(Diag ~ ., data = test_data)[,-ncol(test_data)]
  y_train <- as.numeric(train_data$Diag)
  y_test <- as.numeric(test_data$Diag)
  
  cvfit <- cv.glmnet(x_train, y_train, family = "binomial", alpha = 1)

  lambda_opt <- cvfit$lambda.min

  model <- glmnet(x_train, y_train, family = "binomial", alpha = 1, lambda = lambda_opt)
  

  probs <- predict(model, newx = x_test, type = "Diag")

  roc <- roc(y_test, probs)
  auc <- auc(roc)
  
  return(auc)
}

n_reps <- 100

roc_auc <- numeric(n_reps)

for (i in 1:n_reps) {
  set.seed(123) 
  train_index <- createDataPartition(data[[diagnosis_col]], p = 0.8, list = FALSE)
  train_data <- data_asv[train_index, ]
  test_data <- data_asv[-train_index, ]
  roc_auc[i] <- lasso_log_reg(train_data, test_data)
}

ln_alpha1 <- roc_auc

rf_classification <- function(train_data, test_data) {
  x_train <- model.matrix(Diag ~ ., data = train_data)[,-ncol(train_data)]
  x_test <- model.matrix(Diag ~ ., data = test_data)[,-ncol(test_data)]
  y_train <- as.factor(train_data$Diag)
  y_test <- as.factor(test_data$Diag)
  
  model <- randomForest(x = x_train, y = y_train, ntree = 300, mtry = 3)
  
  probs <- predict(model, newdata = x_test, type = "prob")[, 2]
  
  roc <- roc(y_test, probs)
  auc <- auc(roc)
  
  return(auc)
}

n_reps <- 100

roc_auc <- numeric(n_reps)


for (i in 1:n_reps) {
  set.seed(1024) 
  train_index <- createDataPartition(data[[diagnosis_col]], p = 0.8, list = FALSE)
  train_data <- data_asv[train_index, ]
  test_data <- data_asv[-train_index, ]
  roc_auc[i] <- rf_classification(train_data, test_data)
}
rf_auc <- roc_auc

library(caret)
library(xgboost)
library(pROC)

data <- read.csv("data.csv")

diagnosis_col <- "Diag"
Diag_col <- ncol(data)  # assuming Diag is last column

xgb_classification <- function(train_data, test_data, Diag_col) {
  
  x_train <- train_data[, -Diag_col]
  x_test <- test_data[, -Diag_col]
  y_train <- ifelse(train_data[, Diag_col] == 1, 1, 0)
  y_test <- ifelse(test_data[, Diag_col] == 1, 1, 0)
  
  model <- xgboost(data = as.matrix(x_train), label = y_train, max_depth = 3, eta = 0.1, nrounds = 100, objective = "binary:logistic")
  
  probs <- predict(model, newdata = as.matrix(x_test))
  
  roc <- roc(y_test, probs)
  auc <- auc(roc)
  
  return(auc)
}

n_reps <- 100

roc_auc <- numeric(n_reps)

for (i in 1:n_reps) {
  
  set.seed(123) 
  train_index <- createDataPartition(data[[diagnosis_col]], p = 0.8, list = FALSE)
  train_data <- data[train_index, ]
  test_data <- data[-train_index, ]
  
  roc_auc[i] <- xgb_classification(train_data, test_data, Diag_col)
}
xg_auc <- roc_auc



