# Project:   compareHD_MVN
# Objective: Script to plot the results of the simulation
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-08-03

# Load Results ----------------------------------------------------------

  # inDir <- "../output/"
  # runName <- "20210727_170009"
  #
  # # Read output
  # gg_shape <- readRDS(paste0(inDir, runName, "_res.rds"))
  #
  # # Support Functions
  # source("init.R")

# Plots -------------------------------------------------------------------

  ## Obtain plots
  parm <- unique(gg_shape$par)[1] # what parameter to plot
  result <- unique(gg_shape$variable)[2] # what result to plot
  met_cond <- paste0(unique(gg_shape$method),
                     collapse = "|")
  y_cond <- unique(gg_shape$pm)
  x_cond <- unique(gg_shape$p)

  plot1 <- gg_shape %>%
    # Subset
    filter(grepl(parm, par)) %>%
    filter(grepl(met_cond, method)) %>%
    filter(grepl(result, variable)) %>%

    # Main Plot
    ggplot(aes(x = variable,
               y = value,
               group = method,
               fill = method)) +
    geom_boxplot() +

    # Grid
    facet_grid(rows = vars(factor(p,
                                  labels = paste0("pm = ", y_cond))),
               cols = vars(factor(pm,
                                  labels = paste0("p = ", x_cond)))) +
    # Format
    theme(text = element_text(size = 15),
          plot.title = element_text(hjust = 0.5),
          axis.text = element_text(size = 10),
          axis.text.x = element_blank(),
          axis.title = element_text(size = 15)) +
    labs(title = paste0(result, " for ", parm),
         x     = NULL,
         y     = NULL)

  plot1