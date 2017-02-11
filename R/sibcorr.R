#' Estimate sibling and cousin correlations
#'
#' @param data Estimation data set.
#' @param y Outcome variable
#' @param id1 Individual identifier
#' @param id2 Family identifier
#' @param id3 Extended family (cousin group) identifier
#' @param weight Select one of four weighting schemes.
#' @param controls Control variables to regress out before estimating correlation.
#' @param cousins Estimate cousin correlation if TRUE, otherwise sibling correlation.
#' @import data.table

# Define function to calculate correlations
sibcorr <- function(data, y, id1, id2, id3, weight = 4, controls = NULL, cousins = FALSE) {

  # Make copy of data table so that changes aren't brought out of function scope
  dt <- copy(data)

  # Convert to data table
  setDT(dt)

  # Rename id variables
  setnames(dt, y, "y")
  setnames(dt, id1, "id1")
  setnames(dt, id2, "id2")
  if (cousins == TRUE) {
    setnames(dt, id3, "id3")
  }

  # Set key
  setkey(dt, key = "id2")

  # Prepare data for analysis
  #--------------------------

  if (is.null(controls) == FALSE) {
    # Residualize outcome by regressing on control variables,
    # and take residuals as new outcome variable
    formula <- as.formula(paste("y ~", paste(controls, collapse = " + ")))

    X <- model.matrix(formula, data = dt)
    y <- dt$y

    res <- RcppEigen::fastLmPure(X, y)$residuals

    dt[, e := res]

    # Drop control variables
    dt[, append("y", controls) := NULL]
    setnames(dt,"e", "y")
  }

  # Reshape data into sibling or cousin pairs
  #------------------------------------------

  # Get rid of unneeded variables
  keeplist <- c("id2", "id1", "y")
  if (cousins == TRUE) {
    keeplist <- append("id3", keeplist)
  }
  dt <- dt[, keeplist, with = FALSE]

  # Create all sibling/cousin combinations
  if (cousins == FALSE) {
    byvar <- "id2"
  } else {
    byvar <- "id3"
  }
  dt <- merge(dt, dt, by = byvar, suffix = c(".1", ".2"), allow.cartesian = TRUE)

  # id1 now contains at least one copy of each individual - we tag one of each
  setDT(dt, key = "id1.1")

  # Drop siblings for cousin correlation
  if (cousins == TRUE) dt <- subset(dt, id2.1 != id2.2)

  # Get number of observations used
  n_v <- nrow(unique(dt, by = "id1.1"))
  # Calculate variance
  variance <- var(
    unique(dt, by = "id1.1")$y.1
  )
  # The variance function uses n-1 in the denominator, but Solon
  # uses n - correct for this
  variance <- variance * (n_v - 1) / n_v


  # Calculate family size
  dt[ , n := .N, by = byvar]
  dt[ , n := sqrt(n)]

  # Drop self matches
  dt <- subset(dt, id1.1 != id1.2)

  # Drop duplicate observations
  dt <- subset(dt, id1.1 < id1.2)

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

  return(rho)
}
