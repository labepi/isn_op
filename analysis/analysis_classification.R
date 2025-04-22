bandt_pompe_path <- "~/ufrn/bandt_pompe/"

source(paste(bandt_pompe_path, "bandt_pompe.R", sep = ""))
source(paste(bandt_pompe_path, "measures.R", sep = ""))
source(paste(bandt_pompe_path, "helpers.R", sep = ""))

library(caret)
library(randomForest)

library(MLmetrics)

# parameters
tau <- 1
D_l <- c(2, 3, 4, 5, 6, 7)

N = 10000

filename = paste("../tcp-isn/data/dataset_class_", N, ".dat", sep='')

# loading the dataset
x = read.csv(filename, header=F)

# the number of columns
x_len = ncol(x)

# removing the rows with less samples
x = x[!is.na(x[,x_len]),]

# removing names and class from the dataset
x_class = x[,x_len]
x_names = x[,x_len-1]
x = x[,-c(x_len-1, x_len)]

# adjusting columns types
x = apply(x, 2, as.numeric)

# the matrix to store the features
feats <- matrix(nrow = 0, ncol = (12 + 12 + 12 + 1))
colnames(feats) <- c(
    "H2", "C2", "F2",
    "H2diff", "C2diff", "F2diff",
    "H3", "C3", "F3",
    "H3diff", "C3diff", "F3diff",
    "H4", "C4", "F4",
    "H4diff", "C4diff", "F4diff",
    "H5", "C5", "F5",
    "H5diff", "C5diff", "F5diff",
    "H6", "C6", "F6",
    "H6diff", "C6diff", "F6diff",
    "H7", "C7", "F7",
    "H7diff", "C7diff", "F7diff",
    "label"
)

# computing the features for each time series
# looping nas amostras
for (i in 1:nrow(x)) {
    # features for each serie
    feats_i <- c()

    # computing features for all D's
    for (D in D_l)
    {
        # bandt-pomping the series
        bp <- bandt_pompe_distribution(x[i, ], D = D)
        HP <- shannon_entropy(bp$probabilities, normalized = TRUE)
        SC <- complexity(bp$probabilities, HP)
        FI <- fisher_information(bp$probabilities)
        feats_i <- c(feats_i, HP, SC, FI)

        # bandt-pomping the diff of series
        bp <- bandt_pompe_distribution(diff(x[i, ]), D = D)
        HP <- shannon_entropy(bp$probabilities, normalized = TRUE)
        SC <- complexity(bp$probabilities, HP)
        FI <- fisher_information(bp$probabilities)
        feats_i <- c(feats_i, HP, SC, FI)
    }

    # adding the label
    feats_i <- c(feats_i, x_class[i])

    # adding the features to the matrix
    feats <- rbind(feats, feats_i)
}

# adjusting the features 
feats_df <- as.data.frame(cbind(apply(feats[, 1:36], 2, as.numeric)))
feats_df$label <- as.factor(feats[, "label"])
# rownames(feats_df) = names

# list of best features
feats_list = c("F2", "F3", "F4", "F5", "F6", "F7", "F2diff", "F3diff", "F4diff", "F5diff", "F6diff", "F7diff")

# filtering the best features
feats_best = feats_df[, c(feats_list, "label")]


############################
# testing the classification
############################

set.seed(123)

# spliting train/test 80/20
samples <- createDataPartition(x_class, p = 0.8, list=FALSE)
x_train <- feats_best[samples,]
x_test <- feats_best[-samples,]
y = x_test$label

# KNN

fit <- train(label ~ ., 
                data = x_train, 
                method="knn",  
                tuneGrid=data.frame(k=1:20))

print(fit)


pred <- predict(fit, newdata = x_test)
# confusionMatrix(pred, y)
# acc = sum(pred == y)/nrow(x_test)
acc = Accuracy(pred, y)
f1 = F1_Score(y, pred)

cat("KNN - N: ", N," Acc: ", acc, " F1: ", f1, "\n")


# RF

fit <- train(label ~ ., 
                data = x_train, 
                method="rf")  

print(fit)

pred <- predict(fit, newdata = x_test)
# confusionMatrix(pred, y)
# acc = sum(pred == y)/nrow(x_test)
acc = Accuracy(pred, y)
f1 = F1_Score(y, pred)

cat("RF - N: ", N," Acc: ", acc, " F1: ", f1, "\n")

# SVM Linear

fit <- train(label ~ .,
                data = x_train,
                method="svmLinear")

print(fit)

pred <- predict(fit, newdata = x_test)
# confusionMatrix(pred, y)
# acc = sum(pred == y)/nrow(x_test)
acc = Accuracy(pred, y)
f1 = F1_Score(y, pred)

cat("SVM Linear - N: ", N," Acc: ", acc, " F1: ", f1, "\n")

# SVM Radial

fit <- train(label ~ .,
                data = x_train,
                method="svmRadial")


print(fit)

pred <- predict(fit, newdata = x_test)
# confusionMatrix(pred, y)
# acc = sum(pred == y)/nrow(x_test)
acc = Accuracy(pred, y)
f1 = F1_Score(y, pred)

cat("SVM Radial - N: ", N," Acc: ", acc, " F1: ", f1, "\n")

# SVM Poly

fit <- train(label ~ .,
                data = x_train,
                method="svmPoly")


print(fit)

pred <- predict(fit, newdata = x_test)
# confusionMatrix(pred, y)
# acc = sum(pred == y)/nrow(x_test)
acc = Accuracy(pred, y)
f1 = F1_Score(y, pred)

cat("SVM Poly - N: ", N," Acc: ", acc, " F1: ", f1, "\n")


# TODO: mais classificadores?


# Possible values are found using ‘names(getModelInfo())’. 
# See <http://topepo.github.io/caret/train-models-by-tag.html>. 
# A list of functions can also be passed for a custom model
# function. See <http://topepo.github.io/caret/using-your-own-model-in-train.html>
#
# preProcess: A string vector that defines a pre-processing of the predictor data. Current possibilities are "BoxCox", "YeoJohnson", "expoTrans", "center", "scale", "range", "knnImpute", "bagImpute", "medianImpute", "pca", "ica" and "spatialSign". The default is no pre-processing. See ‘preProcess’ and ‘trainControl’ on the procedures and how to
