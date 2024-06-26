% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/4_maximum_likelihood_estimation.R
\name{ug_control}
\alias{ug_control}
\title{Optimization Control Parameters Passed to optim}
\usage{
ug_control(
  method = "BFGS",
  maxit = 10000,
  start = NULL,
  trace = FALSE,
  hessian = FALSE,
  ...
)
}
\arguments{
\item{method}{the method to be used. See "Details" in \code{\link[stats]{optim}}. The default
method (\code{"BFGS"}) is a quasi-Newton method (also known as a variable metric algorithm),
specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.}

\item{maxit}{the maximum number of iterations of the algorithm. Defaults to \code{10000}.}

\item{start}{An optional vector of initial values to start the algorithm
iterations. The input order should be \code{(beta, gamma, lambda1, lambda2)}, where
\code{beta} is a \emph{k}-dimensional vector with the coefficients associated with the mean,
\code{gamma} is an \emph{l}-dimensional vector with the coefficients associated with the
dispersion parameter, and \code{lambda1} and \code{lambda2} are the parameters associated with
the parametric link functions, if any. Initial guesses for \code{lambda1} or \code{lambda2}
should only be passed if the mean-related or dispersion-related link function is a parametric
link function.}

\item{trace}{non-negative integer. If positive, tracing information on the progress of the
optimization is produced. Higher values may produce more tracing information: for method
\code{"L-BFGS-B"} there are six levels of tracing. (To understand exactly what these do see
 the source code: higher levels give more detail.)}

\item{hessian}{logical; it specifies whether the asymptotic covariance matrix of the estimates
should be obtained based on the numerically differentiated Hessian provided by \code{optim}.
If \code{FALSE}, Fisher's expected information matrix is used.}

\item{...}{further arguments passed to \code{\link[stats]{optim}}.}
}
\value{
A list with the arguments specified.
}
\description{
Optimization parameters passed to \code{\link[stats]{optim}} for the maximum likelihood estimation
    in the unit gamma regression with parametric link functions. This function acts in the same
    spirit as \code{\link[betareg]{betareg.control}} from the \code{betareg} package. Its primary
    purpose is to gather all the optimization control arguments in a single function.
}
\author{
Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
}
