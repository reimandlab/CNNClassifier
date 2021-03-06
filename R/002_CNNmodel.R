setwd("~/oicr/top2b_cnn")
library(keras)
library(onehot)
library(ggplot2)
library(tensorflow)
library(pryr)
library(pROC)


tf$compat$v1$disable_eager_execution()
set.seed(2020)


data <- read.csv("string_500_rc.csv")
# if double == 0, if triple == 1
data$type <- ifelse(nchar(data$order)<=11, 0, 1)
data_double <- data[data$type == 0, ]
data_triple <- data[data$type == 1, ]

# downsample
data <- rbind(data_triple, data_double[sample(1:sum(data$type == 0), sum(data$type == 1)), ])
# split into training set and test set
split.index <- sample(1:nrow(data), 0.8*nrow(data))
data.train <- data[split.index,]
data.test <- data[-split.index,]

c(train_data, train_labels) %<-% data.train[, c("x", "type")]
c(test_data, test_labels) %<-% data.test[, c("x", "type")]

# train: transform into array
max_len = 500
cat("n training data", nrow(data.train), "\n")
array <- array(NA, dim=c(nrow(data.train),max_len,4))
for(i in 1:nrow(data.train)){
	if (i %% 100 == 1) cat(i, " ")
	try <- train_data[i]
	split <- unlist(strsplit(try, split = NULL))
	df <- data.frame(ID = 1: length(split), base = as.factor(split))
	onehot <- predict(onehot(df), df)
	array[i, ,1:4] <- onehot[1:max_len, -1]
}

# reshape the array for CNN
train_array <- array_reshape(array,c(nrow(data.train), max_len, 4))
# ==========================================================
# test
cat("n testing data", nrow(data.test), "\n")
test_array <- array(0, dim=c(nrow(data.test),max_len,4))
  
  for(i in 1: nrow(data.test)){
	if (i %% 100 == 1) cat(i, " ")
    try <- test_data[i]
    split <- unlist(strsplit(try, split = NULL))
    df <- data.frame(ID = 1: length(split), base = as.factor(split))
    # one-hot encoding
    onehot <- predict(onehot(df), df)
    # ACGT
    test_array[i, ,1:4] <- onehot[1:max_len, -1]
  }

# reshape the array
test_array <- array_reshape(test_array, c(nrow(data.test), max_len, 4))
# ============================================================
# ======= build CNN model ====================================
# ============================================================
model <- keras_model_sequential() %>%
  # conv 1
  layer_conv_1d(
    filters = 16,
    kernel_size = 15,
    activation = "relu",
    input_shape = c(max_len, 4)
  ) %>%
  # drop out
  layer_dropout(rate = 0.5) %>%
  
  # conv 2
  layer_conv_1d(filters = 4,
                kernel_size = 15,
                activation = "relu",) %>%
  # drop out
  layer_dropout(rate = 0.5) %>%
  # max_pooling
  layer_max_pooling_1d(pool_size = 2) %>%
  # flatten
  layer_flatten() %>%
  # dense layer
  layer_dense(units = 128, activation = "relu") %>%
  layer_dense(units = 1, activation = "sigmoid")
  
# ======= compile model =================================================
set.seed(2020)
model %>% compile(
  optimizer = optimizer_rmsprop(lr = 1e-3),
  loss = "binary_crossentropy",
  metrics = c("acc")
)


history <- model %>% fit(
  train_array, train_labels,
  epochs = 20,
  batch_size = 128,
  validation_split = 0.2
)

# plot training history
fname = "data/CNN_history.pdf"
pdf(fname)
plot(history)
dev.off()
system(paste("open", fname))

# test model performance
metrics <- model %>% evaluate(test_array, test_labels)
metrics

# ROC curve
fname = "data/CNN_ROC.pdf"
pdf(fname)
gc_prob <- predict(model, x = test_array, type = "prob")
gc_pROC <- roc(response = data.test$type, predictor = gc_prob)
plot(gc_pROC, print.auc=TRUE, grid=TRUE)
gc_pROC$auc

ggroc(gc_pROC, alpha = 0.5, colour = "black", linetype = 1, size = 1)+
  theme_minimal() + ggtitle("ROC curve for CNN model") + 
  geom_segment(aes(x = 1, xend = 0, y = 0, yend = 1), color="grey", linetype="dashed")+
  annotate("text", x = 0.1, y = 0.1, vjust = 0,label = paste("AUC =",sprintf("%.3f",gc_pROC$auc)))+
  theme(plot.title = element_text(hjust = 0.5))
#=========================== CNN ENDS ==========================================================

dev.off()
system(paste("open", fname))


#save(model, file = "data/model.rsav")

