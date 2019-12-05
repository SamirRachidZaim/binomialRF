#' a parallelized version of correlbinom
#'
#' \code{.fast_correlbinom} experimental section on plotting feature selection
#' @param rho what percentage of data are you subsampling
#' @param successprob calculateBinomialP output
#' @param trials ntrees
#' @param precision how many bits in the pmf
#' @param model which model to use
#' @return a ggplot2 object comparison of different models


.fast_correlbinom <- function (rho, successprob, trials, precision = 1024, model = "kuk") 
{
  
  
  if (is.na(rho) | is.na(successprob) | is.na(trials) | is.na(precision) | 
      is.na(model)) {
    stop("NA input detected.")
  }
  if (!is.numeric(rho) | !is.numeric(successprob) | !is.numeric(trials) | 
      !is.numeric(precision)) {
    stop("Non-numeric input detected.")
  }
  if (length(rho) != 1 | length(successprob) != 1 | length(trials) != 
      1 | length(precision) != 1) {
    stop("Non-scalar input detected.")
  }
  if (min(rho, successprob) < 0 | max(rho, successprob) > 1) {
    stop("Rho and success probability must both fall within the unit interval.")
  }
  if (trials <= 0) {
    stop("Number of trials must be positive.")
  }
  if (precision < 2) {
    stop("Precision must be at least two.")
  }
  if (model != "witt" & model != "kuk") {
    stop("Model must equal 'witt' or 'kuk'.")
  }
  if (as.integer(trials) != trials | as.integer(precision) != 
      precision) {
    stop("Trials and precision must both be integer-valued.")
  }
  rho <- Rmpfr::mpfr(rho, precision)
  successprob <- Rmpfr::mpfr(successprob, precision)
  if (model == "witt") {
    pows <- methods::new("mpfr", mapply(`^`, 1 - rho, 0:(trials - 
                                                           1)))
    pow.sum <- cumsum(pows[1:(trials - 1)]) * rho
    cond.prob <- c(Rmpfr::mpfr(1, precision), successprob, 
                   pow.sum + pows[2:trials] * successprob)
    prob.diag <- cumprod(cond.prob)
  }
  else if (model == "kuk") {
    prob.diag <- successprob^((1 - rho) * (0:trials))
  }
  gen.elem <- function(i) {
    coeffs <- Rmpfr::chooseMpfr(trials + 1 - i, 0:(trials + 
                                                     1 - i)) * (-1)^(0:(trials + 1 - i))
    return(sum(prob.diag[i:(trials + 1)] * coeffs) * Rmpfr::chooseMpfr(trials, 
                                                                       i - 1))
  }
  return(Rmpfr::asNumeric(new("mpfr", unlist(parallel::mclapply(1:(trials + 
                                                           1), gen.elem)))))
}

.create_ref_distributions <- function(calculateDistributions){

  ## PROB =1/10 
  system.time(pmf_N500_Rho63 <- .fast_correlbinom(rho = .63,trials = 500, successprob = .1, precision = 1024))
  system.time(pmf_N1000_Rho63 <- .fast_correlbinom(rho = .63,trials = 1000, successprob = .1, precision = 1024))
  system.time(pmf_N2000_Rho63 <- .fast_correlbinom(rho = .63,trials = 2000, successprob = .1, precision = 1024))
  
  pmf_list_prob0.1 <- list(pmf_N500_Rho63=pmf_N500_Rho63,
                           pmf_N1000_Rho63=pmf_N1000_Rho63,
                           pmf_N2000_Rho63=pmf_N2000_Rho63
  )
  
  ## PROB=1/100
  system.time(pmf_N500_Rho63 <- .fast_correlbinom(rho = .63,trials = 500, successprob = .01, precision = 1024))
  system.time(pmf_N1000_Rho63 <- .fast_correlbinom(rho = .63,trials = 1000, successprob = .01, precision = 1024))
  system.time(pmf_N2000_Rho63 <- .fast_correlbinom(rho = .63,trials = 2000, successprob = .01, precision = 1024))
  
  pmf_list_prob0.01 <- list(pmf_N500_Rho63=pmf_N500_Rho63,
                            pmf_N1000_Rho63=pmf_N1000_Rho63,
                            pmf_N2000_Rho63=pmf_N2000_Rho63
  )
  
  ## PROB=1/1000
  system.time(pmf_N500_Rho63  <- .fast_correlbinom(rho = .63,trials = 500, successprob =  .001, precision = 1024))
  system.time(pmf_N1000_Rho63 <- .fast_correlbinom(rho = .63,trials = 1000, successprob = .001, precision = 1024))
  system.time(pmf_N2000_Rho63 <- .fast_correlbinom(rho = .63,trials = 2000, successprob = .001, precision = 1024))
  
  pmf_list_prob0.001 <- list(pmf_N500_Rho63=pmf_N500_Rho63,
                             pmf_N1000_Rho63=pmf_N1000_Rho63,
                             pmf_N2000_Rho63=pmf_N2000_Rho63
  )
  
  
  pmf_list = list(prob0.1=pmf_list_prob0.1,
                  prob0.01=pmf_list_prob0.01,
                  prob0.001=pmf_list_prob0.001)
  
  save(pmf_list, file = '~/Dropbox/Samir/binomialRF/Data/corr_binom_distributions.RData')
}

