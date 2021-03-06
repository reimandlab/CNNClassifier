# hits
# ===========================
library(Biostrings)
# scan each motif through all the sequences
number = 16 # number of motifs
hits <-list()
scan <- function(i, obs){
  countPWM(pwms[[i]], train_data[obs],min.score="85%")
}

for (i in 1:number){
  print(i)
  nhit <- lapply(1:length(train_data), function(obs){scan(i,obs)})
  nhit <- do.call(c,nhit)
  hits[[i]] <- nhit
}



for (motif in 1:16) {
  print(motif)
  # ======= position for triple sites =====
  hit_triple_index <- which((train_labels == 1) & (hits[[motif]] > 0))
  start_triple <- list()
  for (i in 1:length(hit_triple_index)) {
    # print(i)
    p <-
      matchPWM(pwms[[motif]], train_data[hit_triple_index[i]], min.score = "85%")
    start_triple[[i]] <- as.data.frame(p@ranges)$start
  }
  
  pos_triple <- do.call(c, start_triple)
  
  hit_double_index <- which((train_labels == 0) & (hits[[motif]] > 0))
  start_double <- list()
  for (i in 1:length(hit_double_index)) {
    # print(i)
    p <-
      matchPWM(pwms[[motif]], train_data[hit_double_index[i]], min.score = "85%")
    start_double[[i]] <- as.data.frame(p@ranges)$start
  }
  
  pos_double <- do.call(c, start_double)
  
  df_hitpos <-
    data.frame(type = c(rep("triple", length(pos_triple)), rep("double", length(pos_double))),
               position = c(pos_triple, pos_double))
  
  assign(
    paste("p", motif, sep = ""),
    ggplot(df_hitpos, aes(x = position, fill = type)) +
      geom_histogram(
        breaks = seq(0, 500, by = 10),
        position = 'identity',
        alpha = 0.5
      ) +
      scale_fill_manual(values = c("#69b3a2", "#404080")) +
      ggtitle(paste("Motif", motif)) +
      theme(plot.title = element_text(hjust = 0.5)) +
      xlim(c(0, 500))
  )
}
# plot distribution of each motif
fname = "data/motif_positions.pdf"
pdf (fname)
library(gridExtra)
grid.arrange(grobs = list(p1, p2, p3, p4), cols = 4)
grid.arrange(grobs = list(p5, p6, p7, p8), cols = 4)
grid.arrange(grobs = list(p9, p10, p11, p12), cols = 4)
grid.arrange(grobs = list(p13, p14, p15, p16), cols = 4)
dev.off()

system(paste("open ", fname))

