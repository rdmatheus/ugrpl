# This function is a adaptation of the betareg.control() function in betareg package

#' Optimization Control Parameters Passed to optim
#'
#' Optimization parameters passed to \code{\link[stats]{optim}} for the maximum likelihood estimation
#'     in the unit gamma regression with parametric link functions. This function acts in the same
#'     spirit as \code{\link[betareg]{betareg.control}} from the \code{betareg} package. Its primary
#'     purpose is to gather all the optimization control arguments in a single function.
#'
#' @param method the method to be used. See "Details" in \code{\link[stats]{optim}}. The default
#'     method (\code{"BFGS"}) is a quasi-Newton method (also known as a variable metric algorithm),
#'     specifically that published simultaneously in 1970 by Broyden, Fletcher, Goldfarb and Shanno.
#' @param maxit the maximum number of iterations of the algorithm. Defaults to \code{5,000}.
#' @param start An optional vector of initial values to start the algorithm
#'     iterations. The input order should be \code{(beta, gamma, lambda1, lambda2)}, where
#'     \code{beta} is a \emph{k}-dimensional vector with the coefficients associated with the mean,
#'     \code{gamma} is an \emph{l}-dimensional vector with the coefficients associated with the
#'     dispersion parameter, and \code{lambda1} and \code{lambda2} are the parameters associated with
#'     the parametric link functions, if any. Initial guesses for \code{lambda1} or \code{lambda2}
#'     should only be passed if the mean-related or dispersion-related link function is a parametric
#'     link function.
#' @param trace non-negative integer. If positive, tracing information on the progress of the
#'     optimization is produced. Higher values may produce more tracing information: for method
#'     \code{"L-BFGS-B"} there are six levels of tracing. (To understand exactly what these do see
#'      the source code: higher levels give more detail.)
#' @param hessian logical; it specifies whether the asymptotic covariance matrix of the estimates
#'      should be obtained based on the numerically differentiated Hessian provided by \code{optim}.
#'      If \code{FALSE}, Fisher's expected information matrix is used.
#' @param ... further arguments passed to \code{\link[stats]{optim}}.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @return A list with the arguments specified.
#' @export
#'
ug_control <- function(method = "BFGS", maxit = 5000, start = NULL,
                       trace = FALSE, hessian = FALSE, ...)
{
  rval <- list(method = method, maxit = maxit, hessian = hessian,
               trace = trace, start = start)

  rval <- c(rval, list(...))

  if (!is.null(rval$fnscale))
    warning("fnscale must not be modified")

  rval$fnscale <- -1

  if (is.null(rval$reltol))
    rval$reltol <- .Machine$double.eps^(1/1.2)

  rval
}


ug_mle <- function(y, X, Z = NULL, link = "aordaz", sigma.link = NULL,
                   control = ug_control(...), ...)
{

  ## Z matrix definition
  if (is.null(Z)) Z <- matrix(1, nrow = n, ncol = 1L)

  ## Convenience variables
  k <- ncol(X)
  l <- ncol(Z)
  n <- length(y)
  g <- make.plink

  ## Link functions
  if (is.null(sigma.link))  sigma.link <- if (l == 1) "identity" else "aordaz"

  g1.inv <- g(link)$linkinv
  g1.pl  <- g(link)$plink

  g2.inv <- g(sigma.link)$linkinv
  g2.pl  <- g(sigma.link)$plink

  ## Control list
  method <- control$method
  maxit <- control$maxit
  hessian <- control$hessian
  trace <- control$trace
  start  <- control$start
  control$method <- control$hessian <- control$start <- NULL

  ## Initial guess
  lambda10 <- lambda20 <- NULL
  if (g1.pl) lambda10 <- 1
  if (g2.pl) lambda20 <- 1

  if(is.null(start)){
    ystar <- g(link)$linkfun(y, lambda10)
    ystar <- ystar[which(is.finite(ystar))]
    Xaux <- X[which(is.finite(ystar)), ]
    start <- c(solve(t(Xaux)%*%Xaux)%*%t(Xaux)%*%ystar,
               g(sigma.link)$linkfun(0.3, lambda20),  rep(0, l - 1),
               lambda10, lambda20)
  }

  ## Estimates
  opt <- suppressWarnings(stats::optim(par = start,
                                       fn = ll_ug,
                                       gr = U_ug,
                                       y = y, X = X, Z = Z,
                                       link = link,
                                       sigma.link = sigma.link,
                                       method = method,
                                       control = control,
                                       hessian = hessian))

  if (opt$convergence > 0)
    warning(cat("optimization failed to converge\n"))

  opt$start <- start
  opt
}
