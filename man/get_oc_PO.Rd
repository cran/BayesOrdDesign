% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/OC_func_PO.R
\name{get_oc_PO}
\alias{get_oc_PO}
\title{Generate operating characteristics for Bayesian two-stage trial design
of ordinal endpoints with proportional odds assumption}
\usage{
get_oc_PO(alpha, pro_ctr, nmax, fixed_es, ormax, fixed_ss, ntrial, method)
}
\arguments{
\item{alpha}{the desirable type I error rate to be controlled}

\item{pro_ctr}{distribution of clinical categories for the
control group}

\item{nmax}{the maximum sample size for operating characteristics}

\item{fixed_es}{fixed effect size when simulate the OC for various sample
size}

\item{ormax}{the maximum effect size for OC}

\item{fixed_ss}{fixed sample size when simulate the OC for various effect
size}

\item{ntrial}{the number of simulated trials}

\item{method}{whether the statistical test for interim/final analysis is Bayesian or
Frequentist. method = "Frequentist" for Frequentist approach; method = "Bayesian"
for Bayesian approach}
}
\value{
get_oc_PO() returns the operating characteristics of design as a
table, including: (1) user-defined value, either sample size or effect size
(2) corresponding power (3) average sample size
}
\description{
Obtain operating characteristics (OC) of the Bayesian two-stage trial
design of ordinal endpoints with proportional odds assumption.
}
\details{
Grid search of sample size is used for guarantee a desirable type I error rate.
  The upper limitation is 200, and lower limitation default is sample size 50
  for the control and treatment groups at each stage. Default increment of the
  sequence is 10.

  For the parameter estimation section, we have two options, and can be selected using
  the method argument.Two following options are available: (i) method = "Frequentist",
  (ii) method = "Bayesian". If method = "Frequentist", parameters are estimated via package
  ordinal, which is based on frequentist method, while method = "Bayesian", parameters are
  estimated through Bayesian model.

  Two types of operating characteristics can be implemented through this function.

  Please note, in our example, argument ntrial = 5 is for the time saving purpose.
}
\examples{

get_oc_PO(alpha = 0.05, pro_ctr = c(0.58,0.05,0.17,0.03,0.04,0.13),
          ormax = 1.5, fixed_ss = 150,
          ntrial = 5, method = "Frequentist")


get_oc_PO(alpha = 0.05, pro_ctr = c(0.58,0.05,0.17,0.03,0.04,0.13),
          nmax = 200, fixed_es = 1.5,
          ntrial = 5, method = "Frequentist")

}
