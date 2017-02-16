#' Estimate sibling and cousin correlations with boostrap confidence intervals.
#'
#' \code{sibcorr_bs} estimates sibling or cousin correlations with block bootstrap standard errors.
#'
#' @inheritParams sibcorr
#' @param ... Other options for \code{\link{sibcorr}}
#' @param reps Number of bootstrap replications
#' @param ci_level Set level for bootstrap confidence interval
#' @param cores Set number of processor cores to use for bootstrap estimation
#' @return The estimated sibling or cousin correlation coefficient and bootstrap confidence interval
#' @details The formula must be specified as \code{outcome ~ controls | individual + family + ext_family},
#' where \code{individual} is an individual identifier, \code{family} is a family (sibling group) identifier,
#' and \code{ext_family} is an extended family (cousin group) identifier.
#' The extended family identifier is only required if \code{cousins = TRUE}.
#'
#' The formula can not include functions. This means that all variable transformations,
#' including setting some variables as factors, must be performed before estimation.
#' @import data.table
#' @import foreach

# Define function to bootstrap the correlation
sibcorr_bs <- function(formula, data, ..., cousins = FALSE, reps = 50, ci_level = 95, cores = 1) {

  # Parse formula
  identifiers <- all.vars(formula(Formula::Formula(formula), lhs = 0, rhs = 2))

  # Make copy of data table so that changes aren't brought out of function scope
  dt <- copy(data)

  # Convert to data table
  setDT(dt)

  # Get point estimate
  rho <- sibcorr(formula, data = dt, cousins = cousins, ...)

  # Sample rows with replacement
  if (cousins == FALSE) {
    byvar <- identifiers[2]
  } else {
    byvar <- identifiers[3]
  }
  groups <- unique(dt, by = byvar)[, get(byvar)]


  # Draw a bootstrap sample and estimate the sibling correlation
  bs_draw <- function(...) {
    # Sample clusters with replacement
    clusters <- sample(groups, length(groups), replace = TRUE)
    # Get resampled dataset
    resample <- dt[get(byvar) %in% clusters]
    # Calculate correlation
    result <- sibcorr(..., data = resample)
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
      bs_draw(formula, cousins = cousins, ...)
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
      .export = c("sibcorr"),
      .packages = c("data.table")
    ) %dopar% {
      # Perform a bootstrap draw
      bs_draw(formula, cousins = cousins, ...)
    }

    parallel::stopCluster(cl)
  }

  # Calculate confidence interval quantiles
  ci_limits <- c((1 - ci_level/100)/2, 1 - (1 - ci_level/100)/2)

  # Calculate bootstrap confidence interval
  ci <- quantile(bs, ci_limits)

  result <- c(rho, ci)

  #names(result) <- c("estimate", "ci_lower", "ci_upper")

  return(result)
}
