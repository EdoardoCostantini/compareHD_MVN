# Project:   compareHD_MVN
# Objective: helper functions for the various algorithms used in this project
# Author:    Edoardo Costantini
# Created:   2021-08-02
# Modified:  2021-08-02

## make list of collinear variables to remove
## Credits: mice R package
findCollinear <- function(x, threshold = 0.999, ...) {
  nvar <- ncol(x)
  x <- data.matrix(x)
  r <- !is.na(x)
  nr <- apply(r, 2, sum, na.rm = TRUE)
  ord <- order(nr, decreasing = TRUE)
  xo <- x[, ord, drop = FALSE]
  varnames <- dimnames(xo)[[2]]
  z <- suppressWarnings(cor(xo, use = "pairwise.complete.obs"))
  hit <- outer(seq_len(nvar), seq_len(nvar), "<") & (abs(z) >= threshold)
  out <- apply(hit, 2, any, na.rm = TRUE)
  return(varnames[out])
}

## Generates a starting filled in dataset for the gibbs sampler
initializeImp <- function(data){
  ## Internals
  # data <- boys

  ## Body
  # Extract variable types
  target_names <- names(which(colSums(is.na(data)) != 0)) # names target variables
  target_n <- length(target_names) # number of variables needing imputation
  vartypes <- sapply(1:target_n, function(j){
    class(data[, target_names[j]])[1]
  })
  cont_vars <- target_names[vartypes == "numeric" | vartypes == "integer"]
  fact_vars <- target_names[vartypes == "factor" | vartypes == "ordered"]

  # Impute random draw from observed data for for numeric variables
  if(length(cont_vars) != 0){
    for(j in cont_vars){
      var_j <- data[, j]
      wy <- is.na(var_j) # missing
      ry <- !wy # observed
      data[wy, j] <- sample(var_j[ry], sum(wy), replace = TRUE)
    }
  }

  # Impute random draw from observed data for (any) factor
  if( length(fact_vars) != 0 ){
    for (j in fact_vars) {
      j <- fact_vars[1]
      var_j <- data[, j]
      wy <- is.na(var_j) # missing
      data[wy, j] <- sample(levels(var_j), sum(wy), replace = TRUE)
    }
  }

  return(data)  # an initialized dataset
}