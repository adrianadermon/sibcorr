#' Estimate sibling and cousin correlations
#'
#' @param data Estimation data set.
#' @param weight Select one of four weighting schemes.
#' @param controls Control variables to regress out before estimating correlation.
#' @param cousins Estimate cousin correlation if TRUE, otherwise sibling correlation.
#' @param reps Number of bootstrap replications.
#' @param ci_level Set level for bootstrap confidence interval.

# Define function to bootstrap the correlation
sibcorr_bs <- function(data, id1, id2, id3, ..., cousins = FALSE, reps = 50, ci_level = 95) {

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

  bs <- list()

  print(paste("Performing", reps, "bootstrap replications...", sep = " "))

  # Create progress bar
  pb <- txtProgressBar(min = 0, max = reps, style = 3)

  for (i in 1:reps) {
    # Update progress bar
    setTxtProgressBar(pb, i)
    # Sample clusters with replacement
    clusters <- sample(groups, length(groups), replace = TRUE)
    # Get resampled dataset
    resample <- dt[get(byvar) %in% clusters]
    # Calculate correlation
    result <- sibcorr(resample, id1 = id1, id2 = id2, id3 = id3, cousins = cousins, ...)
    # Store result
    bs <- append(bs, result)
  }

  # Close progress bar
  close(pb)

  # Flatten list
  bs <- unlist(bs)

  # Calculate confidence interval quantiles
  ci_limits <- c((1 - ci_level/100)/2, 1 - (1 - ci_level/100)/2)

  # Calculate bootstrap confidence interval
  ci <- quantile(bs, ci_limits)

  result <- c(rho, ci)

  #names(result) <- c("estimate", "ci_lower", "ci_upper")

  return(result)
}
