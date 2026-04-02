#' Generation-time noise utilities
#'
#' @description
#' Functions that add stochastic jitter to parameters during generation.
#' This is NOT differential privacy -- it creates realistic variability
#' so outputs are not exact replicas of the parameter tables.
#'
#' All noise is seed-controlled and optional.

#' Jitter a probability vector on the logit scale, then re-normalize
#'
#' @param probs Numeric vector of probabilities
#' @param noise A noiseConfig object
#' @return Jittered probability vector summing to 1
#' @keywords internal
jitterProbabilities <- function(probs, noise) {
  if (noise$probabilityJitterSd <= 0 || length(probs) == 0) {
    return(probs)
  }
  eps <- 1e-8
  p <- pmin(pmax(probs, eps), 1 - eps)
  logits <- log(p / (1 - p))
  logits <- logits + rnorm(length(logits), 0, noise$probabilityJitterSd)
  expLogits <- exp(logits - max(logits))
  expLogits / sum(expLogits)
}


#' Add uniform jitter to a duration in days
#'
#' @param baseDays Integer duration
#' @param noise A noiseConfig object
#' @return Jittered duration (minimum 1)
#' @keywords internal
jitterDuration <- function(baseDays, noise) {
  if (noise$durationJitterDays <= 0L) {
    return(baseDays)
  }
  delta <- sample(
    seq(-noise$durationJitterDays, noise$durationJitterDays),
    size = 1
  )
  max(1L, baseDays + delta)
}


#' Shift an administration day within a cycle
#'
#' @param day Integer day offset
#' @param noise A noiseConfig object
#' @return Jittered day (minimum 0)
#' @keywords internal
jitterCycleDay <- function(day, noise) {
  if (noise$cycleJitterDays <= 0L) {
    return(day)
  }
  delta <- sample(
    seq(-noise$cycleJitterDays, noise$cycleJitterDays),
    size = 1
  )
  max(0L, day + delta)
}


#' Decide whether to drop an exposure record (missingness)
#'
#' @param noise A noiseConfig object
#' @return Logical
#' @keywords internal
shouldDropExposure <- function(noise) {
  if (noise$missingnessRate <= 0) {
    return(FALSE)
  }
  runif(1) < noise$missingnessRate
}
