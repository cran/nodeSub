## ----setup, include=FALSE-----------------------------------------------------
library(ggplot2)
library(magrittr)
library(tidyr)
library(nodeSub)

knitr::opts_chunk$set(echo = TRUE)

## ----create alignment---------------------------------------------------------
seq_length <- 100
sub_rate <- 1 / seq_length

input_tree <- TreeSim::sim.bd.taxa(n = 100,
                                   numbsim = 1,
                                   lambda = 1,
                                   mu = 0.1,
                                   complete = TRUE)[[1]]

target_alignment <- sim_unlinked(phy = input_tree,
                                 rate1 = sub_rate,
                                 rate2 = sub_rate,
                                 l = seq_length,
                                 node_time = 0.3)

## ----twin alignment-----------------------------------------------------------
comp_alignment <- create_equal_alignment(
                        input_tree = geiger::drop.extinct(input_tree),
                        # can only work on trees without extinct branches
                        sub_rate = sub_rate,
                        alignment_result = target_alignment)

## ----infer_phylogeny----------------------------------------------------------

if (beastier::is_beast2_installed()) {

  node_posterior <- infer_phylogeny(target_alignment$alignment,
                                    "node_posterior",
                                    clock_prior = beautier::create_strict_clock_model(clock_rate_param =    beautier::create_clock_rate_param(value = sub_rate)),
                                    chain_length = 1e5,
                                    burnin = 0.1,
                                    working_dir = getwd(),
                                    sub_rate = sub_rate)

  reference_posterior <- infer_phylogeny(comp_alignment$alignment,
                                         "reference_posterior",
                                         burnin = 0.1,
                                         clock_prior = beautier::create_strict_clock_model(clock_rate_param = beautier::create_clock_rate_param(value = comp_alignment$adjusted_rate)),
                                         chain_length = 1e5,
                                         working_dir = getwd,
                                         sub_rate = sub_rate) 
}

## ----stats--------------------------------------------------------------------
if (beastier::is_beast2_installed()) {
node_stats <- calc_sum_stats(node_posterior$all_trees,
                             input_tree)
ref_stats  <- calc_sum_stats(reference_posterior$all_trees,
                            input_tree)

node_stats$differences$method <- "node_sub"
ref_stats$differences$method  <- "reference"

results <- rbind(node_stats$differences,
                 ref_stats$differences)
}

## ----plot results-------------------------------------------------------------
if (beastier::is_beast2_installed()) {
  results %>%
  tidyr::gather(key = "statistic", value = "val", -method) %>%
  dplyr::filter(statistic != "num_tips") %>%  
  ggplot(aes(x = method , y = val)) +
    geom_boxplot() +
    facet_wrap(~statistic, scales = "free")
}

