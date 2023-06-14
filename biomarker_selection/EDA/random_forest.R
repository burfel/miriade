library(glmnet)
library(caret)
library(randomForest)

library(here)
source(here("functions", "data_processing_functions.R"))

datasets_root_directory <- define_datasets_root()

# Read data
kth_data <- read.table(file.path(datasets_root_directory, "KTH/KTH AD dataset for MIRIADE biomarker selection sample info and results.csv"),
                  sep = ";", header = T, quote = "")


# Remove missing values
kth_data_NoNA <- na.omit(kth_data)

# Split the dataset into training and testing sets
set.seed(123)
training.samples <- kth_data_NoNA$Diagnosis %>% 
  createDataPartition(p = 0.8, list = FALSE)
train.data <- kth_data_NoNA[training.samples, ]
test.data <- kth_data_NoNA[-training.samples, ]

# Scale the features
preproc <- preProcess(train, method = c("center", "scale"))
trainTransformed <- predict(preproc, train)
testTransformed <- predict(preproc, test)


# Set up the Lasso regression model using the glmnet package:
x <- model.matrix(Species ~ ., train)[,-1]
y <- as.matrix(train$Species)
lasso <- cv.glmnet(x, y, alpha = 1, nfolds = 10)

# Identify the optimal value of the penalty parameter using cross-validation:
plot(lasso)
bestLambda <- lasso$lambda.min

# Use the optimal value of the penalty parameter to fit the Lasso model on the training data:
lassoModel <- glmnet(x, y, alpha = 1, lambda = bestLambda)

# Use the Lasso model to make predictions on the test data: 
xtest <- model.matrix(Species ~ ., test)[,-1]
ytest <- as.matrix(test$Species)
predictions <- predict(lassoModel, newx = xtest)

accuracy <- mean(predictions == ytest)

# # Train the model
# rf <- randomForest(Class ~ ., data = train.data, ntree = 500, mtry = 3, importance = TRUE)

# Evaluate the model
predicted <- predict(rf, test.data)
table(predicted, test.data$Class)
accuracy <- sum(predicted == test.data$Class) / nrow(test.data)
accuracy
confusionMatrix(predictions, test$Class)