#' Determine the sample size for Bayesian two-stage trial design
#' of ordinal endpoints without proportional odds assumption
#'
#' @description
#'
#' Obtain estimated sample size based on user-specified type I
#' error, power and effect size defined by the odds ratio between
#' the treatment and control groups, without the proportional
#' odds (PO) assumption.
#'
#'
#' @param nmax the maximum sample size for searching to get the desirable power
#' @param or_alt effect size to be detected (under H_1) in
#' terms of odds ratio
#' @param pro_ctr distribution of clinical categories for the
#' control group
#' @param U the desirability of each outcome level
#' @param alpha the desirable type I error rate to be controlled
#' @param power the desirable power to be achieved
#' @param ntrial the number of simulated trials
#' @param method whether the statistical test for interim/final analysis is Bayesian or
#' Frequentist. method = "Frequentist" for Frequentist approach; method = "Bayesian"
#' for Bayesian approach
#'
#'
#' @details
#'   Grid search of sample size is used for guarantee a desirable type I error rate.
#'   The upper limitation is 200, and lower limitation default is sample size 50
#'   for the control and treatment groups at each stage. Default increment of the
#'   sequence is 50.
#'
#'   For the parameter estimation section, we have two options, and can be selected using
#'   the method argument.Two following options are available: (i) method = "Frequentist",
#'   (ii) method = "Bayesian". If method = "Frequentist", parameters are estimated via package
#'   ordinal, which is based on frequentist method, while method = "Bayesian", parameters are
#'   estimated through Bayesian model.
#'
#'   Specifically, the numerical utilities U reflect the desirability of each outcome
#'   level. To do this, in our example, we first set U[1] = 100 and U[5] = 0, and then
#'   asked physicians to specify numerical values for the intermediate levels, that
#'   reflect their desirability relative to the best and worst levels.
#'
#'   Please note, in our example, argument ntrial =  5 is for the time saving purpose.
#'
#'
#'
#' @return ss_npo() returns recommended sample size for each
#' of two groups for the interim and final stages, by assuming 1:1
#' equal randomization for the two groups at each stage; and corresponding power.
#'
#' @export
#'
#' @examples
#' set.seed(123)
#' ss_npo(nmax = 100, or_alt = c(1.8,1.5,1,1,1),
#'        pro_ctr = c(0.5,0.1,0.05,0.2,0.05,0.1), U = c(100,80,65,25,10,0),
#'        alpha = 0.05, power = 0.8, ntrial = 5, method = "Frequentist")
#'


#----------------------------------------------------------------------------
## sample size calculator
ss_npo = function(nmax, or_alt, pro_ctr, U,
                      alpha, power,ntrial, method){

  N = 200
  # under null, calculate thresholds
  cf_grid        = seq(0.6, 0.7, by=0.05)
  threshold_grid = seq(0.75, 0.85, by=0.05)

  or_null  = rep(1, length(pro_ctr)-1)
  or.mat = matrix(rep(or_null, ntrial), nrow = ntrial,
                  ncol = length(or_null), byrow = TRUE)
  output = c()

  for (cf in cf_grid){
    for (threshold in threshold_grid){

      out = multiple_trial_npo(or.mat, sd = 0.2, pro_ctr, U, n = N, cf=cf,
                               threshold=threshold, method = method)
      rr = c(cf, threshold, out)
      output = rbind(output, rr)
      colnames(output) = c("cf", "threshold", "PET(%)", "alpha", "avgss")
      results = as.data.frame(output)
    }
  }
  index = min(which(abs(results$alpha-alpha)==min(abs(results$alpha-alpha))))
  vec = c(results[index,c(1,2)])
  thrsh = c(vec$cf, vec$threshold)
  names(thrsh) = c("futility", "superority")

  # calculate power
  or.mat = matrix(rep(or_alt, ntrial), nrow = ntrial,
                  ncol = length(or_alt), byrow = TRUE)
  n_grid = seq(50, nmax, by = 50)
  output = c()

  for (n in n_grid){
    out = multiple_trial_npo(or.mat, sd = 0.2, pro_ctr, U, n=n, cf=vec$cf,
                             threshold=vec$threshold, method = method)
    rr = c(n, out)
    output = rbind(output, rr)
    colnames(output) = c("samplesize", "PET(%)", "Power", "avgss")
  }
  results = list()
  index = min(which(abs(output[,3]-power)==min(abs(output[,3]-power))))
  results$sample_size = output[index, 1]
  results$power = output[index, 3]
  results$threshold = thrsh
  return(results)
}
