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
  ref <- out %>%
    # Subset
    filter(grepl("OG", method)) %>%
    group_by(p, pm, par) %>%
    summarize(ref = mean(est))

  out_plus <- merge(out, ref, by = c("pm", "p", "par"))
  out_plus_sc <- out_plus[, c(colnames(out), "ref")]
  out_plus_sc$par <- factor(out_plus_sc$par,
                                 levels = unique(out$par))
  out_sr <- out_plus_sc %>% arrange(rp, p, pm, method, par)

# Add coverage
  out_ggready <- out_sr %>%
    mutate(CIC = ci.lower < ref & ref < ci.upper)

# Shape for plots
  gg_shape <- reshape2::melt(out_ggready, id.var = colnames(out)[1:6])
  head(gg_shape)
  head(out_ggready)