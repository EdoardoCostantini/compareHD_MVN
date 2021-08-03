# Project:   compareHD_MVN
# Objective: Pool the results of a simulation study given its id folder
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-08-03

## Make sure we have a clean environment:
rm(list = ls())

## Initialize the environment:
source("./init.R")

# Load Results ----------------------------------------------------------

  inDir <- "../output/"
  files <- list.files(inDir)
  target_tar <- files[length(files)]
  output <- readTarGz(target_tar)

# Restructure Results -----------------------------------------------------
# list of conditions containing results for every repetition

  out <- do.call(rbind, output$out)
  gg_shape <- reshape2::melt(out, id.var = colnames(out)[1:6])

# Add reference value and coverage

  true_vec <- out %>%
    # Subset
    filter(grepl("OG", method)) %>%
    group_by(p, pm, par) %>%
    summarize(ref = mean(est))

  out_augmented <- merge(out, true_vec, by = c("pm", "p", "par"))
  out_augmented_sc <- out_augmented[, c(colnames(out), "ref")]
  out_augmented_sc$par <- factor(out_augmented_sc$par,
                                 levels = unique(out$par))
  output_sr <- out_augmented_sc %>% arrange(rp, p, pm, method, par)

# Add coverage
  output_sr <- output_sr %>%
    mutate(CIC = ci.lower < ref & ref < ci.upper)

# Shape for plots
  gg_shape <- reshape2::melt(output_sr, id.var = colnames(out)[1:6])