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
#' @param noise A noise_config object
#' @return Jittered probability vector summing to 1
#' @keywords internal
jitter_probabilities <- function(probs, noise) {
  if (noise$probability_jitter_sd <= 0 || length(probs) == 0) {
    return(probs)
  }
  eps <- 1e-8
  p <- pmin(pmax(probs, eps), 1 - eps)
  logits <- log(p / (1 - p))
  logits <- logits + rnorm(length(logits), 0, noise$probability_jitter_sd)
  # Softmax-style normalization
  exp_logits <- exp(logits - max(logits))
  exp_logits / sum(exp_logits)
}


#' Add uniform jitter to a duration in days
#'
#' @param base_days Integer duration
#' @param noise A noise_config object
#' @return Jittered duration (minimum 1)
#' @keywords internal
jitter_duration <- function(base_days, noise) {
  if (noise$duration_jitter_days <= 0L) {
    return(base_days)
  }
  delta <- sample(
    seq(-noise$duration_jitter_days, noise$duration_jitter_days),
    size = 1
  )
  max(1L, base_days + delta)
}


#' Shift an administration day within a cycle
#'
#' @param day Integer day offset
#' @param noise A noise_config object
#' @return Jittered day (minimum 0)
#' @keywords internal
jitter_cycle_day <- function(day, noise) {
  if (noise$cycle_jitter_days <= 0L) {
    return(day)
  }
  delta <- sample(
    seq(-noise$cycle_jitter_days, noise$cycle_jitter_days),
    size = 1
  )
  max(0L, day + delta)
}


#' Decide whether to drop an exposure record (missingness)
#'
#' @param noise A noise_config object
#' @return Logical
#' @keywords internal
should_drop_exposure <- function(noise) {
  if (noise$missingness_rate <= 0) {
    return(FALSE)
  }
  runif(1) < noise$missingness_rate
}
