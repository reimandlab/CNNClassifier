double_index <- which(train_labels == 0)
triple_index <- which(train_labels == 1)
# prediction for double and triple sites
double_pred <- predict(model, train_array)[double_index]
triple_pred <- predict(model, train_array)[triple_index]

TP_index <- triple_index[triple_pred > 0.5]
TN_index <- double_index[double_pred < 0.5]

output <- model$output[, 1]
input_layer <- model$input
# compute gradients
grads <- k_gradients(output, input_layer)[[1]]
# mean in each filter
# pooled_grads <- k_mean(grads, axis = c(1,2))
# iterate function
iterate <-
  k_function(list(model$input),
             list(grads, input_layer[1, , ]))
# ==================================================================================
score_triple <- matrix(NA, length(TP_index), 500)
for (j in 1:length(TP_index)) {
  seq <- train_array[TP_index[j], , , drop = F]
  c(grads_value, layer_output_value) %<-% iterate(list(seq))
  # saliency score = activation * gradients
  score <- grads_value[1, , ] * layer_output_value
  for (i in 1:500) {
    score_triple[j, i] <- score[i,][which(layer_output_value[i, ] == 1)]
  }
}

# ggplot saliency score
df <-
  data.frame(
    position = 1:500,
    obs = rep(1:dim(score_triple)[1], each = 500),
    score = as.vector(t(score_triple))
  )
  

fname = "data/saliency_map.pdf"
pdf(fname)

# plot saliency map
p1 = ggplot(df, aes(position, obs, fill = score)) +
  geom_tile() +
  scale_fill_distiller(limit = c(-0.02, 0.02), palette = "RdYlBu") +
  ggtitle("Saliency Score for Triple Sites")

# average score
df <- data.frame(pos = 1:500, score = colMeans(score_triple))
# plot average score
p2 = ggplot(df, aes(x = pos, y = score)) + geom_point() +
  geom_smooth(method = 'auto') +
  ggtitle("Average Score for Triple Sites")
  
print(p1)
print(p2)


dev.off()
system(paste("open", fname))