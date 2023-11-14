## ----setup, include=FALSE-----------------------------------------------------
library(nodeSub)
library(ggplot2)
knitr::opts_chunk$set(echo = TRUE)

plot_phyDat <- function(phyDat_alignment) {
  vx <- as.data.frame(phyDat_alignment)

  num_entries <- nrow(vx) * ncol(vx)
  to_plot <- matrix(nrow = num_entries, ncol = 3)

  cnt <- 1
  for (i in 1:nrow(vx)) {
    for (j in 1:ncol(vx)) {
      to_plot[cnt, ] <- c(i, j, vx[i, j])
      cnt <- cnt + 1
    }
  }
  colnames(to_plot) <- c("x", "y", "base")
  to_plot <- tibble::as_tibble(to_plot)
  to_plot$x <- as.numeric(to_plot$x)
  to_plot$y <- as.numeric(to_plot$y)
  ggplot(to_plot, aes(x = x, y = y, fill = base)) +
    geom_tile()
}


## ----start--------------------------------------------------------------------
seq_length <- 30
sub_rate <- 1 / seq_length

input_tree <- TreeSim::sim.bd.taxa(n = 10,
                                   numbsim = 1,
                                   lambda = 1,
                                   mu = 0.1,
                                   complete = TRUE)[[1]]

normal_alignment <- sim_normal(input_tree,
                               l = seq_length,
                               rate = sub_rate)

plot_phyDat(normal_alignment$alignment)

## ----sim linkedunlinked-------------------------------------------------------
unlinked_alignment <- sim_unlinked(input_tree,
                                   rate1 = sub_rate,
                                   rate2 = sub_rate,
                                   l = seq_length,
                                   node_time = 0.5)

plot_phyDat(unlinked_alignment$alignment)

linked_alignment <- sim_linked(input_tree,
                               rate = sub_rate,
                               node_mut_rate_double = sub_rate * sub_rate,
                               node_time = 0.5,
                               l = seq_length)

plot_phyDat(linked_alignment$alignment)

## ----sim explicit-------------------------------------------------------------
unlinked_explicit <- sim_unlinked_explicit(input_tree,
                                   rate1 = sub_rate,
                                   rate2 = sub_rate,
                                   l = seq_length,
                                   node_time = 0.5)

plot_phyDat(unlinked_explicit$alignment)

normal_explicit <- sim_normal_explicit(input_tree,
                                        l = seq_length,
                                        rate = sub_rate)

plot_phyDat(normal_explicit$alignment)

