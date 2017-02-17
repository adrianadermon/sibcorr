#' Estimate sibling and cousin correlations.
#'
#' @param formula Three-part formula describing the outcome variable, control variables to be regressed out, and identifiers for the individual, immediate family, and extended family
#' @param data Estimation data set.
#' @param weight Select one of four weighting schemes.
#' @details The formula must be specified as \code{outcome ~ controls | individual + family + ext_family},
#' where \code{individual} is an individual identifier, \code{family} is a family (sibling group) identifier,
#' and \code{ext_family} is an extended family (cousin group) identifier.
#' If \{ext_family} is omitted, the sibling correlation is estimated -
#' otherwise, the cousin correlation is estimated
#'
#' The formula does not handle functions on the left hand side.
#' This means that any transformations of the outcome variable must be performed before estimation.
#' @return The estimated sibling or cousin correlation coefficient
#' and number of individuals, sibling or cousin pairs, and families or extended families
#' @import data.table

# Define function to calculate correlations
sibcorr <- function(formula, data, weight = 4) {

  # Put formula in Formula format
  formula <- Formula::Formula(formula)

  # Parse formula
  y <- all.vars(terms(formula, rhs = 0))
  controls <- all.vars(terms(formula, lhs = 0, rhs = 1))
  controls_transformed <- labels(terms(formula, lhs = 0, rhs = 1))
  identifiers <- all.vars(terms(formula, lhs = 0, rhs = 2))

  # Use number of id variables to select sibling or cousin correlation
  len_id <- length(identifiers)

  # Ensure that the correct number of id variables have been given
  if (len_id != 2 & len_id != 3) stop("formula must include two or three identifier variables")


  # Make copy of data table so that changes aren't brought out of function scope
  dt <- copy(data)

  # Convert to data table
  setDT(dt)

  # Rename id variables
  setnames(dt, y, "y")
  setnames(dt, identifiers[1], "id1")
  setnames(dt, identifiers[2], "id2")
  if (len_id == 3) {
    setnames(dt, identifiers[3], "id3")
  }

  # Set key
  setkey(dt, key = "id2")

  # Prepare data for analysis
  #--------------------------

  if (length(controls) != 0) {
    # Residualize outcome by regressing on control variables,
    # and take residuals as new outcome variable
    formula <- as.formula(paste("y ~", paste(controls_transformed, collapse = " + ")))

    X <- model.matrix(formula, data = dt)
    y <- dt$y

    res <- RcppEigen::fastLmPure(X, y)$residuals

    dt[, e := res]

    # Drop control variables
    dt[, append("y", controls) := NULL]
    setnames(dt, "e", "y")
  }

  # Estimate variance
  #------------------

  # Get number of observations used
  n_ind <- nrow(dt)
  # Calculate variance
  variance <- var(dt$y)
  # The variance function uses n-1 in the denominator, but Solon
  # uses n - correct for this
  variance <- variance * (n_ind - 1) / n_ind

  # Reshape data into sibling or cousin pairs
  #------------------------------------------

  # Get rid of unneeded variables
  keeplist <- c("id2", "id1", "y")
  if (len_id == 3) {
    keeplist <- append("id3", keeplist)
  }
  dt <- dt[, keeplist, with = FALSE]

  # Create all sibling/cousin combinations
  if (len_id == 2) {
    byvar <- "id2"
  } else {
    byvar <- "id3"
  }
  dt <- merge(dt, dt, by = byvar, suffix = c(".1", ".2"), allow.cartesian = TRUE)

  # Drop siblings for cousin correlation
  if (len_id == 3) dt <- subset(dt, id2.1 != id2.2)

  # Calculate family size
  dt[ , n := .N, by = byvar]
  dt[ , n := sqrt(n)]

  # Drop self matches
  dt <- subset(dt, id1.1 != id1.2)

  # Drop duplicate observations
  dt <- subset(dt, id1.1 < id1.2)

  # Count number of unique families
  n_fams <- nrow(unique(dt, by = byvar))

  # Now we have a dataset with one copy of each unique sibling/cousin pair

  # Calculate weights
  #------------------

  # Weighting schemes (Hällsten, footnote 8)
  # w1 (all families weighted equally, strongest down-weighting of large families) = [1/2n(n - 1)]^{-1}
  # w2(strong down-weighting) = [1/2(n - 1)]^{-1}
  # w3 (weak down-weighting) = [1/2[n(n - 1)]^{1/2}]^{-1}
  # w4 (equal cousin pair weight, no down-weighting of large families) = 1

  # Calculate the selected weights
  if (weight == 1) {
    # Calculate w1 weights
    w <- (1/2 * dt$n * (dt$n - 1))
  } else if (weight == 2) {
    # Calculate w2 weights
    w <- (1/2 * (dt$n - 1))
  } else if (weight == 3) {
    # Calculate w3 weights
    w <- (1/2 * (dt$n * (dt$n - 1))^(1/2))
  } else if (weight == 4) {
    w <- rep(1, nrow(dt))
  } else
    stop("Weight must be an integer between 1 and 4")


  # Estimate covariances
  #---------------------

  # Get rid of unneeded variables
  dt <- dt[, .(y.1, y.2)]

  # Calculate weighted covariance
  covariance <- cov.wt(dt, wt = 1/w, method = "ML")[["cov"]]["y.1", "y.2"]

  # Calculate correlation
  rho <- covariance / variance

  # Get number of sibling pairs
  n_pairs <- nrow(dt)

  results <- c(correlation = rho, n_individuals = n_ind, n_pairs = n_pairs, n_families = n_fams)
  return(results)
}
