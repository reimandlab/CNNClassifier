library("ggseqlogo")
library("Biostrings")
dat <- score_triple
TP_string <- train_data[TP_index]

pos <- which(dat > quantile(dat, 0.999), arr.ind = TRUE)
# get all the subsequences
sub <- rep(NA, nrow(pos))
# collect all high score subsequences
for (i in 1:nrow(pos)){
  sub[i] <- substr(TP_string[pos[i, 1]], pos[i, 2] - 7, pos[i, 2] + 7)
}

ggseqlogo(sub)
# =============== sequence clustering ===============
library(kmer)
library(insect)
library(dendextend)
library(msa)
set.seed(2020)
subclust <- char2dna(sub, simplify = FALSE)
names(subclust) <- 1:length(sub)

fname = "data/motif_tree.pdf"
pdf(fname, width=12, height=7)
tree <- cluster(subclust, k = 5, residues = "DNA")
#hc <- as.hclust(tree)
#plot(hc, main = "cluster sequences")
plot(tree, main = "cluster sequences")
dev.off()
system(paste("open", fname))

# cut 
clusters <- 4
#ct <- cutree(hc,clusters)
ct <- cutree(tree, clusters)
for (i in 1:clusters){
  assign(paste("c",i, sep = ""), sub[which(tree == i)])
}
# ============== mutiple sequence alignment =========================
m1 <- msa(DNAStringSet(c1), method = "ClustalW", gapOpening = 20， gapExtension= 20, order="aligned")
m2 <- msa(DNAStringSet(c2), method = "ClustalW", gapOpening = 20， gapExtension= 20, order="aligned")
m3 <- msa(DNAStringSet(c3), method = "ClustalW", gapOpening = 20， gapExtension= 20, order="aligned")
m4 <- msa(DNAStringSet(c4), method = "ClustalW", gapOpening = 20， gapExtension= 20, order="aligned")

# plot alignment



con1 <- consensusMatrix(m1)
cond1 <- con1[1:4,] 
ggseqlogo(cond1)

con2 <- consensusMatrix(m2)
cond2 <- con2[1:4,] 
ggseqlogo(cond2)

con3 <- consensusMatrix(m3)
cond3 <- con3[1:4,] 
ggseqlogo(cond3)

con4 <- consensusMatrix(m4)
cond4 <- con4[1:4,] 
ggseqlogo(cond4)
# choose the middle area of the alignment
