# Project:   compareHD_MVN
# Objective: subroutine runCell a single condition for a single repetition
# Author:    Edoardo Costantini
# Created:   2021-07-30
# Modified:  2021-07-30

runCell <- function(rp, cond, parms, fs) {

# Example Internals -------------------------------------------------------

  # rp   = 1
  # cond = conds[1, ]

# Data Generation ---------------------------------------------------------

  ## Data according to condition
  data_OG <- dataGen(parms = parms, cond = cond)

  ## Impose Missingness according to condition
  data_miss <- amputePerVar(data = data_OG,
                            targets = parms$vmap_miss$ta,
                            preds = parms$vmap_miss$mp,
                            pm = cond$pm,
                            type = "high")

# Imputation --------------------------------------------------------------

  if(!cond$method %in% c("OG", "cca")){
    if(cond$method %in% "blasso"){
      mids_out <- imputeBlasso(boys,
                               m = parms$blasso_ndt,
                               maxit = parms$blasso_iters)
    } else {
      mids_out <- mice(data_miss,
                       method = as.vector(cond$method),
                       m = parms$mice_ndt,
                       maxit = parms$mice_iters,
                       printFlag = TRUE)
    }
  }

# Analyze and pool --------------------------------------------------------

  # Are we doing OG/cca or MI?
  if(cond$method %in% c("OG", "cca")){
    # Are we doing OG or cca?
    if(cond$method %in% "OG"){
      data_noMI <- data_original
    } else {
      data_noMI <- na.omit(data_miss)
    }
      # Fit on complete data
      og_sat_fit <- miFitSat(mi_data = list(data_noMI),
                             model = satModWrite(names(data_noMI[,
                                                         parms$vmap_miss$ta]))
      )[[1]]
      og_est_all <- parameterEstimates(og_sat_fit, standardized = TRUE)
      og_est <- og_est_all[, c("std.all", "est", "ci.lower", "ci.upper")]
      result <- cbind(par = apply(og_est_all[, 1:3], 1, paste0, collapse = ""),
                      og_est,
                      fmi = NA,
                      riv = NA)
  } else {
    mi_sat_fits <- miFitSat(mi_data = complete(mids_out, "all"),
                            model = satModWrite(names(data_miss[,
                                                        parms$vmap_miss$ta]))
    )
    mi_sat_pool <- miPool(mi_fits = mi_sat_fits,
                          m = parms$mice_ndt,
                          N = parms$n)
    result <- mi_sat_pool
  }

  # Attach descriptor
  row.names(cond) <- NULL # to avoid a warning in the cbind
  result <- cbind(rp = rp, cond, result)

# Store Output ------------------------------------------------------------

  ## Store Results
  saveRDS(result,
          file = paste0(fs$outDir,
                        "rep_", rp, "_", cond$tag,
                        ".rds")
  )

}