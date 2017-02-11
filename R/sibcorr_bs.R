#' Estimate sibling and cousin correlations
#'
#' @inheritParams sibcorr
#' @param ... Other options for \code{\link{sibcorr}}
#' @param reps Number of bootstrap replications.
#' @param ci_level Set level for bootstrap confidence interval.
#' @param cores Set number of processor cores.
#' @import data.table
#' @import foreach

# Define function to bootstrap the correlation
sibcorr_bs <- function(data, id1, id2, id3 = "", ..., cousins = FALSE, reps = 50, ci_level = 95, cores = 1) {

  # Make copy of data table so that changes aren't brought out of function scope
  dt <- copy(data)

  # Convert to data table
  setDT(dt)

  # Get point estimate
  rho <- sibcorr(dt, id1 = id1, id2 = id2, id3 = id3, cousins = cousins, ...)

  # Sample rows with replacement
  if (cousins == FALSE) {
    byvar <- id2
  } else {
    byvar <- id3
  }
  groups <- unique(dt, by = byvar)[, get(byvar)]


  # Draw a bootstrap sample and estimate the sibling correlation
  bs_draw <- function(id1 = id1, id2 = id2, id3 = id3, cousins = cousins, ...) {
    # Sample clusters with replacement
    clusters <- sample(groups, length(groups), replace = TRUE)
    # Get resampled dataset
    resample <- dt[get(byvar) %in% clusters]
    # Calculate correlation
    result <- sibcorr(resample, id1 = id1, id2 = id2, id3 = id3, cousins = cousins, ...)
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
      bs_draw(id1 = id1, id2 = id2, id3 = id3, cousins = cousins, ...)
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
      bs_draw(id1 = id1, id2 = id2, id3 = id3, cousins = cousins, ...)
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
