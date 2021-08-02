# Project:   compareHD_MVN
# Objective: initialization script
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-07-30

# Packages ----------------------------------------------------------------

  pack_list <- c("parallel",
                 "MASS",
                 "lavaan",
                 "rlecuyer",
                 "forcats",
                 "stringr",
                 "ggplot2",
                 "dplyr",
                 "mice")

  lapply(pack_list, library, character.only = TRUE, verbose = FALSE)

# Load Functions ----------------------------------------------------------

  # Subroutines
  all_subs <- paste0("./subroutines/",
                     list.files("./subroutines/"))
  lapply(all_subs, source)

  # Functions
  all_funs <- paste0("./functions/",
                     list.files("./functions/"))
  lapply(all_funs, source)

  # Helper
  all_help <- paste0("./helper/",
                     list.files("./helper/"))
  lapply(all_help, source)

# Fixed Parameters --------------------------------------------------------

  # Empty List
  parms    <- list()

  # Simulation
  parms$rps <- 1e3
  parms$rep_counter <- 0
  parms$seed <- 2021
  parms$nStreams <- 1000

  # Data generation
  parms$n <- 2e2 # sample size
  parms$item_mean <- 5 # 5 # true item mean
  parms$item_var  <- (2.5)^2 # true item variance

  # Variable Maps
  parms$vmap_miss <- list(ta = 1:4, # TArget ampute
                          mp = 5:6) # mar predictors
  parms$rho_high <- .6 # correlation for highly correlated variables
  parms$rho_low  <- .1 # correlation for correlated variables

  # Imputation Routine
  parms$mice_ndt <- 2
  parms$mice_iters <- 3
  parms$blasso_ndt <- parms$mice_ndt
  parms$blasso_iters <- 3

  # Storing location
  parms$outDir <- "../output/"

# Conditions --------------------------------------------------------------

  p   <- c(40, 200)
  pm <- c(.1, .3)
  method <- c("durr.gaus", "iurr.gaus",
              "blasso",
              "norm", # bridge
              "pcr.boot",
              "cart", "rf",
              "OG", "cca")

  conds <- expand.grid(p, pm, method)
    colnames(conds) <- c("p", "pm", "method")