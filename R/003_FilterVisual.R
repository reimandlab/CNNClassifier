#setwd("~/oicr/top2b_cnn")
#library(keras)
#library(onehot)
#library(ggplot2)
#library(tensorflow)
#library(pryr)
#library(pROC)
#
#print(load(file = "data/model.rsav"))

# ==============================================
input_tensor <- model$input
activations_fn <- k_function(inputs = list(input_tensor),
                             outputs = list(model$layers[[1]]$output))
activations <- activations_fn(list(train_array))[[1]]
# ================================
# Find maximally-activating sequences for motif plotting
number = 16
test_seqs = train_array
filter_len = 15
threshold=0.8

all_dna_strings <- list()
for (filt in 1:dim(activations)[3]){ 
	# activations for filter 
	act <- activations[,,filt]
	act[act < threshold*max(act)] <- 0
	pos.values <- which(act > 0,arr.ind=TRUE)
	test_seqs <- train_data

	f <- function(seqq) {
	    temp_str <- test_seqs[pos.values[seqq,1]]
	    temp_idx <- pos.values[seqq,2]
	    strseq <- substr(temp_str, temp_idx, temp_idx+(filter_len-1))
	    return(strseq)
	  }
	
	seqs <- lapply(1:nrow(pos.values), f)
	dna_strings <- do.call(c,seqs)
	all_dna_strings[[filt]] <- dna_strings
}

# number of sequences for each motif
lapply(1:number, function(x) {
  return(length(all_dna_strings[[x]]))
         }
  )
# get pwm for all motifs
library(Biostrings)
library(seqLogo)
pwms <- list()
pfms <- list()
for (i in 1:number){
  dna <- all_dna_strings[[i]]
  dna <- dna[nchar(dna)==filter_len]
  pfm <- consensusMatrix(dna)
  pwm <- PWM(pfm)
  # output
  pfms[[i]] <-consensusMatrix(dna, TRUE)[1:4,]
  pwms[[i]]<- pwm
}

names(pfms) <- names(pwms) <- paste("motif" , 1:number)
names(all_dna_strings) <- paste("motif", 1:number)
# ===================================
# plot sequenece logos using ggseqLogo
# ==================================
library(ggseqlogo)
lapply(1:number, function(i) {
  assign(paste("p", i, sep = ""),
         ggplot() + geom_logo(all_dna_strings[[i]]) + theme_logo() + ggtitle(paste("motif", i)))
})

# arrange the motifs
fname = "data/motifs.pdf"
pdf(fname)
print(ggseqlogo(all_dna_strings, ncol=4))
dev.off()
system(paste("open", fname))

