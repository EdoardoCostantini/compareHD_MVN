# Project:   compareHD_MVN
# Objective: MI with blasso method non-mice function w/ mids output
# Author:    Edoardo Costantini
# Created:   2021-08-02
# Modified:  2021-08-02

imputeBlasso <- function(Z, m, maxit){

  # Prep data ---------------------------------------------------------------
  # Z = data_miss
  # Z = nhanes
  # m = parms$blasso_ndt
  # maxit = parms$blasso_iters

    tryCatch({

      # Dry run of mice for sum internal elements
      mice_dry <- mice(Z, m = 1, maxit = 0)
      method <- replace(mice_dry$method, mice_dry$method != "", "blasso")

      # Save NA action state
      default.na.action <- options('na.action')

      # Change to pass for model.matrix generation
      options(na.action = 'na.pass')

      # Generate model matrix
      Z_mm <- model.matrix(~., Z)[, -1]

      # Revert NA action
      options(na.action = default.na.action$na.action)

      # Get rid of constants (might happen for 1 unpopular dummy code)
      const <- names(which(apply(Z_mm, 2, var) == 0))
      Z_mm  <- Z_mm[, !colnames(Z_mm) %in% const]

      # Find Dummies that have 95% of observations in 1 category
      tabular <- apply(Z_mm, 2, table)

      # Select only dummies
      tabular <- tabular[sapply(tabular, length) == 2]

      # Vector of dummy names to discard
      dum.disc <- lapply(tabular, function(x) {
        x[1] / sum(x) > .95 | x[1] / sum(x) < 1-.95
      })
      dum.disc <- names(which(dum.disc == TRUE))
      Z_mm    <- Z_mm[, !colnames(Z_mm) %in% dum.disc]

      # Find collinear variables
      coll_vars <- findCollinear(Z_mm)
      Z_mm  <- data.frame(Z_mm[, !colnames(Z_mm) %in% coll_vars])

      # Generate observed/missing data locator object O
      O <- as.data.frame(!is.na(Z_mm))

      # Obtain data characteristics of model.matrix data version
      target_n    <- sum(colMeans(O) < 1)
      target_id <- names(which(colMeans(O) < 1))
      nmis        <- colSums(!O)
      nr       <- colSums(!O[,colMeans(O) < 1])

      # To store imputed values and check convergence
      imp_blasso_val <- vector("list", m)
        names(imp_blasso_val) <- seq(1:m)

      # mids element: imp
      imp <- lapply(1:ncol(Z), function(j){
        data.frame(
          matrix(data = NA,
                 nrow = nmis[j],
                 ncol = m),
          check.names = FALSE
        )
      }
      )

      # mids element: chainMean and chainVar
      chainMean <- chainVar <- array(NaN,
                                     dim = c(ncol(Z), maxit, m),
                                     dimnames = list(colnames(Z),
                                                     NULL,
                                                     paste(rep("Chain", m), 1:m)))

      # For one chain of imputatuions
      for (cc in 1:m){
        # cc <- 1
        # Initialize data per chain
        Zm <- initializeImp(Z_mm)
        for(j in 1:target_n){
          J <- which(colnames(Zm) %in% target_id[j])
          wy <- !O[, J]
          chainMean[J, 1, cc] <- mean(Zm[wy, J])
          chainVar[J, 1, cc] <- var(Zm[wy, J])
        }

        # Storing multiply imputed datasets (in from the last chain)
        imp_blasso_dat <- vector("list", maxit)
          names(imp_blasso_dat) <- seq(1:maxit)
        imp_blasso_dat$`1` <- Zm # imputations from initialization

        # Empty storing objects for MCMC samples
        # For each, I need one slot per imputation model / target variable
        # regression coefficients
        beta_out <- lapply(1:target_n, matrix, data = NA, nrow = maxit,
                           ncol = ncol(Zm)-1)
        for (i in 1:target_n) beta_out[[i]][1, ] <- rep(0, ncol(Zm)-1)

        # Error variance
        sigma_out <- lapply(1:target_n, matrix, data = NA, nrow = maxit,
                          ncol = 1)
        for (i in 1:target_n) sigma_out[[i]][1, ] <- 1

        # tau parameter in Hans 2009 and 10 blasso model
        tau_out <-  lapply(1:target_n, matrix, data = NA, nrow = maxit,
                           ncol = 1)
        for (i in 1:target_n) tau_out[[i]][1, ] <- 1

        # phi parameter in Hans 2010 blasso model
        phi_out <-  lapply(1:target_n, matrix, data = NA, nrow = maxit,
                           ncol = 1)
        for (i in 1:target_n) phi_out[[i]][1, ] <- .5

        # Imputed scores
        imp_out <- lapply(target_id, function(x) {
          matrix(data = NA, nrow = maxit, ncol = nr[x],
                 dimnames = list(NULL, rownames(Zm[!O[, x],]) ))
        })
        for (i in 1:target_n) imp_out[[i]][1, ] <- Zm[!O[, target_id[i]],
                                                      target_id[i]]

        # Loop across iterations
        for (i in 2:maxit) {
          Z_out <- Z # each iteration will output 1 dataset
          print(paste0("BLASSO - Chain: ",
                       cc, "/", i, "; Iter: ",
                       i, "/", maxit,
                       " at ", Sys.time()))
          # Loop across variables (cycle)
          for (j in 1:target_n) {
            J <- which(colnames(Zm) %in% target_id[j])
            zj_obs <- Zm[O[, J], J]
            Z_obs  <- Zm[O[, J], -J]
            Z_mis  <- Zm[!O[, J], -J]

            # Scaled versions
            s_zj_obs <- scale(zj_obs)
            s_Z_obs  <- scale(Z_obs)
            s_Z_mis  <- scale(Z_mis,
                              center = colMeans(Z_obs),
                              scale = apply(Z_obs, 2, sd))

            # Storing draws
            beta_m  <- beta_out[[j]][i-1, ]
            sigma_m <- sqrt(sigma_out[[j]][i-1])
            tam_m   <- tau_out[[j]][i-1]
            phi_m   <- phi_out[[j]][i-1]

            # step 1: update j-th imputation model parameters
            pdraw <- blasso::blasso.vs(Y = s_zj_obs, X = s_Z_obs,
                                       iters = 1,
                                       burn  = 0,
                                       beta  = beta_m,
                                       beta.prior = "scaled",
                                       sig2 = sigma_m, sig2prior = c(a = .1, b = .1),
                                       tau  = tam_m,   tauprior  = c(r = .01, s = .01),
                                       phi  = phi_m,   phiprior  = c(h = 1, g = 1),
                                       fixsig = FALSE, fixtau = FALSE, fixphi = FALSE,
                                       noisy = FALSE)

            # step 2: sample from predictive distribution (predict on beta draws scale)
            s_pdraw_zj_imp <- rnorm(nrow(s_Z_mis),
                                    mean = (s_Z_mis %*% as.vector(pdraw$beta)),
                                    sd = sqrt(pdraw$sig2) )
            # and rescale
            pdraw_zj_imp <- s_pdraw_zj_imp * sd(zj_obs) + mean(zj_obs)

            # Store draws
            beta_out[[j]][i, ] <- as.vector(pdraw$beta)
            sigma_out[[j]][i]    <- pdraw$sig2
            tau_out[[j]][i]    <- ifelse(pdraw$tau > 0, pdraw$tau, tau_out[[j]][i-1])
            # In unlikely scenario where Tau is smaller than 0, keep previous draw
            # [not sound]
            phi_out[[j]][i]    <- pdraw$phi
            imp_out[[j]][i, ]  <- pdraw_zj_imp

            # Append imputation (for next iteration)
            Zm[!O[, J], J] <- pdraw_zj_imp
            imp[[J]][, cc] <- pdraw_zj_imp
            chainMean[J, i, cc] <- mean(pdraw_zj_imp)
            chainVar[J, i, cc] <- var(pdraw_zj_imp)
          }
          # Replace imputed columns with imputated ones
          Z_out[, target_id] <- Zm[, target_id]
          imp_blasso_dat[[i]] <- Z_out
        }
        imp_blasso_val[[cc]] <- imp_out
      }

      # Cast as mids object
      midsobj <- list(
        data = Z, # Original (incomplete) data set.
        imp = imp,
        m = m,
        where = mice_dry$where,
        blocks = mice_dry$blocks,
        call = "mice(data = Z)",
        nmis = mice_dry$nmis,
        method = method,
        predictorMatrix = mice_dry$predictorMatrix,
        visitSequence = mice_dry$visitSequence,
        formulas = mice_dry$formulas,
        post = mice_dry$post,
        blots = mice_dry$blots,
        ignore = mice_dry$ignore,
        seed = mice_dry$seed,
        iteration = maxit,
        lastSeedValue = .Random.seed,
        chainMean = chainMean,
        chainVar = chainVar,
        loggedEvents = mice_dry$loggedEvents,
        version = packageVersion("mice"),
        date = Sys.Date()
      )
      oldClass(midsobj) <- "mids"
      return(midsobj)

    }, error = function(e){
      err <- paste0("Original Error: ", e)
      print(err)
      return(mice_dry)
    }
    )
}