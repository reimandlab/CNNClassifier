setwd("~/oicr/top2b_cnn")
library(BSgenome.Hsapiens.UCSC.hg19)
library(seqTools)


# data normalization
load('data/0817_sites_mapped_w100_md65_cleaned.rsav',
     verbose = T)
df <- sites


# keep only double sites (CTCF & RAD21) and triple sites (CTCF & RAD21 & TOP2B)
df2 <- df[df$general_order == "CTCF-RAD21", ]
df3 <- df[df$general_order == "CTCF-RAD21-TOP2B", ]
data <- rbind(df2, df3)
data$length <- data$end - data$start + 1
# mid point
data$mid <- round(rowMeans(data[, c("start", "end")]))
# normalize data (take the midpoint and extend both ways) 4
max_len = 500
data$start <- data$mid - max_len / 2 + 1
data$end <- data$mid + max_len / 2

ranges <- makeGRangesFromDataFrame(data[, 1:3], )
# fetch the DNA content of each of the sites
string <- getSeq(Hsapiens, ranges)

# if any sequence with N, delete
string3 <-
  data.frame(
    x =  as.data.frame(string),
    order = data$mapped_order,
    complete = NA
  )
for (i in 1:nrow(string3)) {
  string3$complete[i]  <-
    !any(strsplit(string3$x[i], split = "")[[1]] == "N")
}

string3 <- string3[string3$complete, 1:2]
string4 <-
  data.frame(x = as.data.frame(reverseComplement(string)),
             order = data$mapped_order)
for (i in 1:nrow(string4)) {
  string4$complete[i]  <-
    !any(strsplit(string4$x[i], split = "")[[1]] == "N")
}
string4 <- string4[string4$complete, 1:2]

string_500_rc <- rbind(string3, string4)
write.csv(string_500_rc, file = "string_500_rc.csv", row.names = FALSE)

