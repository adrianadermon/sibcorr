#' Estimate sibling and cousin correlations.
#'
#' \code{sibcorr} estimates sibling or cousin correlations with optional block bootstrap standard errors.
#'
#' @param formula Three-part formula describing the outcome variable, control variables to be regressed out, and identifiers for the individual, immediate family, and extended family
#' @param data Estimation data set.
#' @param weight Select one of four weighting schemes.
#' @param restriction Put a restriction on the pairs to be used for estimating the covariance
#' @param variance Estimate variance on all individuals ("pre", the default) or only restricted sample ("post")
#' @param reps Number of bootstrap replications - if set to zero, confidence intervals are not estimated
#' @param ci_level Set level for bootstrap confidence interval
#' @param cores Set number of processor cores to use for bootstrap estimation
#' @return The estimated sibling or cousin correlation coefficient,
#' number of individuals, sibling or cousin pairs, and families or extended families,
#' and optionally a bootstrap confidence interval
#' @details The formula must be specified as \code{outcome ~ controls | individual + family + ext_family},
#' where \code{individual} is an individual identifier, \code{family} is a family (sibling group) identifier,
#' and \code{ext_family} is an extended family (cousin group) identifier.
#' If \code{ext_family} is omitted, the sibling correlation is estimated -
#' otherwise, the cousin correlation is estimated
#'
#' The formula does not handle functions on the left hand side.
#' This means that any transformations of the outcome variable must be performed before estimation.
#'
#' The restriction is specified as a vector with the first element giving the name of the variable to restrict on.
#' If the second element is 0, only pairs with the same value are used.
#' If the second element is a positive integer, only pairs with that specific difference are used.
#' If the second element is "unequal", only pairs with different values for the variable are used.
#' If two integers are given (as second and third elements), and the first is smaller than the second,
#' only pairs with a difference within that range (inclusive) are used;
#' if the second is smaller than the first,
#' only pairs with a difference outside that range (exclusive) are used.
#'
#' @import data.table
#' @import foreach


# Define function to bootstrap the correlation
sibcorr <- function(formula, data, weight = 4, restriction = NULL, variance = "pre", reps = 0, ci_level = 95, cores = 1) {

  # Put formula in Formula format
  forml <- Formula::Formula(formula)

  # Parse formula
  y <- all.vars(terms(forml, rhs = 0))
  controls <- all.vars(terms(forml, lhs = 0, rhs = 1))
  identifiers <- all.vars(terms(forml, lhs = 0, rhs = 2))

  # Use number of id variables to select sibling or cousin correlation
  len_id <- length(identifiers)

  # Ensure that the correct number of id variables have been given
  if (len_id != 2 & len_id != 3) stop("formula must include two or three identifier variables")

  # Make copy of data so that changes aren't brought out of function scope
  # Keep only relevant variables
  keep_vars <- c(y, controls, identifiers)
  if (is.null(restriction) == FALSE) {
    keep_vars <- append(restriction[1], keep_vars)
  }
  if ("data.table" %in% class(data)) {
    dt <- copy(data[, keep_vars, with = FALSE])
  } else {
    dt <- copy(data[, keep_vars])
  }

  # Convert to data table
  setDT(dt)

  # Drop incomplete rows
  dt <- na.omit(dt)

  # Get point estimate
  rho <- sc(formula, data = dt, weight = weight, restriction = restriction, variance = variance)

  if (reps == 0) {
    result <- rho
  } else if (reps > 0) {
    # Sample rows with replacement
    if (len_id == 2) {
      byvar <- identifiers[2]
    } else {
      byvar <- identifiers[3]
    }
    groups <- unique(dt, by = byvar)[, get(byvar)]

    # Draw a bootstrap sample and estimate the sibling correlation
    bs_draw <- function(...) {
      # Sample clusters with replacement
      clusters <- sample(groups, length(groups), replace = TRUE)
      # Get bootstrap sample
      resample <- dt[J(clusters), allow.cartesian = TRUE, on = byvar]
      # Calculate correlation
      result <- sc(..., data = resample)[1]
      return(result)
    }

    # Perform bootstrap using one or many processor cores
    if (cores == 1) {
      print(paste("Performing", reps, "bootstrap replications...", sep = " "))

      # Create progress bar
      pb <- txtProgressBar(min = 0, max = reps, style = 3)

      bs <- foreach(i = 1:reps,
                    .combine = "c") %do% {
        # Update progress bar
        setTxtProgressBar(pb, i)
        # Perform a bootstrap draw
        bs_draw(formula, weight = weight, restriction = restriction, variance = variance)
      }

      # Close progress bar
      close(pb)
    } else {
      cl <- parallel::makeCluster(cores)
      doParallel::registerDoParallel(cl)

      print(paste("Performing", reps, "bootstrap replications on", cores, "processor cores", sep = " "))
      bs <- foreach(
        i = 1:reps,
        .combine = "c",
        .export = c("sc", "data"),
        .packages = c("data.table")
      ) %dopar% {
        # Perform a bootstrap draw
        bs_draw(formula, weight = weight, restriction = restriction, variance = variance)
      }

      parallel::stopCluster(cl)
    }

    # Calculate confidence interval quantiles
    ci_limits <- c((1 - ci_level/100)/2, 1 - (1 - ci_level/100)/2)

    # Calculate bootstrap confidence interval
    ci <- quantile(bs, ci_limits)

    result <- c(rho, ci)
  }
  return(result)
}
