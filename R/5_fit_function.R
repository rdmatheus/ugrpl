#' @name ugrpl
#'
#' @title Unit Gamma Regression with Parametric Link Functions
#'
#' @description Fit the unit gamma regression with parametric link functions via maximum likelihood
#'  for a parameterization of this distribution that is indexed by mean and dispersion parameters.
#'
#' @param formula simbolic description of the model, of type
#'  \code{y ~ x} for covariates in the mean model only or \code{y ~ x | z}
#'   to enter covariates in the dispersion model. See details below.
#' @param data,subset,na.action  arguments controlling formula processing via
#' \link[stats]{model.frame}.
#' @param link,sigma.link character specification of the link function for the
#'  mean and the dispersion submodels, respectively. For the mean submodel,
#'  Aranda-Ordaz (\code{"aordaz"}) is the default link function. If it is not assumed that
#'  the dispersion submodel depends on covariates (i.e., is \code{sigma.link} is missing),
#'  the default link is the identity (\code{"identity"}), otherwise, the Aranda-Ordaz link function.
#'  The list of currently available links can be consulted in the details below.
#' @param y,x logicals. If \code{TRUE} the corresponding components of the fit,
#'      response and model matrices, are returned.
#' @param control a list of control parameters passed as arguments for
#'     the \code{\link[stats]{optim}} function specified via \code{\link{ug_control}}.
#' @param ... arguments passed to \code{\link{ug_control}}.
#'
#' @details This implementation uses a parameterization of the unit gamma distribution indexed by
#'     the mean (\code{mu}) and a dispersion parameter (\code{sigma}) in which both \code{mu} and
#'     \code{sigma} take values on (0, 1), see \code{\link{dugamma}}. It is assumed that the mean
#'     depends on covariates through a link function, which can belong to a family of parametric link
#'     functions or be one of the customarily used link functions (e.g., logit or probit). Moreover,
#'     since the dispersion parameter sigma takes values on (0, 1), it is equally possible to consider
#'     a regression structure using a parametric link function.
#'
#'     The basic formula is of type \code{y ~ x1 + x2 + ... + xk} which specifies the model for the
#'     mean response only. Following the syntax of the \code{betareg} package (Cribari-Neto and Zeileis, 2010),
#'     the model for the dispersion index, say in terms of \code{z1, z2, ..., zl}, is specified as
#'     \code{y ~ x1 + x2 + ... + xk | z1 + z2 + ... + zl} using functionalities inherited from package
#'     \code{Formula} (Zeileis and Croissant, 2010).
#'
#'     We assume that the link functions belonging to a parametric family are indexed by a positive
#'     parameter \code{lambda}. When the link function does not belong to this family (e.g., the logit
#'     function), then, by default, \code{lambda = NULL}. The available link functions are
#'
#'   \tabular{llc}{
#'  \bold{Link function}  \tab \bold{Abbreviation} \tab \bold{Is it a parametric link function?}\cr
#'  Logit \tab \code{"logit"} \tab \code{FALSE} \cr
#'  Probit \tab \code{"probit"} \tab \code{FALSE} \cr
#'  Cauchit \tab \code{"cauchit"} \tab \code{FALSE} \cr
#'  Log-Log \tab \code{"loglog"} \tab \code{FALSE} \cr
#'  Complement log-log \tab \code{"cloglog"} \tab \code{FALSE} \cr
#'  Identity \tab \code{"identity"} \tab \code{FALSE} \cr
#'  Aranda-Ordaz \tab  \code{"aordaz"} \tab \code{TRUE} \cr
#'  Power logit \tab  \code{"plogit"} \tab \code{TRUE} \cr
#'  Power pobit \tab \code{"pprobit"} \tab \code{TRUE}  \cr
#'  Power cauchit \tab \code{"pcloglog"} \tab \code{TRUE}  \cr
#'  Power log-log \tab \code{"ploglog"} \tab \code{TRUE}  \cr
#'  Power complement log-log \tab \code{"pcloglog"} \tab \code{TRUE}\cr
#'  Reversal power logit \tab \code{"rplogit"} \tab \code{TRUE}  \cr
#'  Reversal power pobit \tab \code{"rpprobit"} \tab \code{TRUE}  \cr
#'  Reversal power cauchit \tab \code{"rpcauchit"} \tab \code{TRUE}  \cr
#'  Reversal power log-log \tab \code{"rploglog"} \tab \code{TRUE} \cr
#'  Reversal power complement log-log \tab \code{"rpcloglog"} \tab \code{TRUE} \cr
#'  Reversal Aranda-Ordaz \tab \code{"raordaz"} \tab \code{TRUE}
#'  }
#'
#'
#' @return The \code{ugrpl} function returns an object of class \code{"ugrpl"},
#'      which consists of a list with the following components:
#'  \describe{
#'     \item{coefficients}{ a list containing the elements "mean" and
#'         "dispersion" that consist of the estimates of the coefficients
#'          associated with the mean and the dispersion, respectively.}
#'     \item{lambda}{ a list with the estimates of parameters associated with
#'         parametric link functions. If a non-parametric link is used,
#'         \code{NULL} is returned.}
#'     \item{sigma}{ a vector with the fitted dispersion parameters.}
#'     \item{link, sigma.link}{ link function specified for the mean and the dispersion, respectively.}
#'     \item{fitted.values}{ a vector with the fitted means.}
#'     \item{logLik}{ log-likelihood of the fitted model.}
#'     \item{vcov}{ asymptotic covariance matrix of the maximum likelihood
#'         estimators of all parameters in the model. By default, the asymptotic
#'         covariance matrix is based on Fisher's information matrix, but can
#'         be obtained from the Hessian matrix (obtained numerically via \code{\link[stats]{optim}})
#'         if \code{hessian = TRUE}.}
#'     \item{residuals}{ a vector of quantile residuals.}
#'     \item{convergence}{ logical; if \code{TRUE}, it indicate successful convergence.}
#'     \item{start}{ the initial values of the optimization algorithm.}
#'     \item{optim}{ output from the \link[stats]{optim} call for maximizing the log-likelihood.}
#'     \item{control}{ the control arguments passed to the \link[stats]{optim} call.}
#'     \item{nobs}{ number of observations.}
#'     \item{df.null}{ residual degrees of freedom in the null model (constant mean and dispersion), that is, n - 2.}
#'     \item{df.residual}{ residual degrees of freedom in the fitted model.}
#'     \item{call}{ the function call.}
#'     \item{formula}{ the formula used to specify the model in \code{ugrpl}.}
#'     \item{terms}{ a list with elements "mean", "dispersion" and "full"
#'          containing the terms objects for the respective models,}
#'     \item{y}{ the response vector (if \code{y = TRUE}).}
#'     \item{x}{ a list with elements "mean" and "dispersion" containing the
#'          model matrices from the respective models (if \code{x = TRUE}).}
#'    }
#'
#'
#' @references Cribari-Neto F, Zeileis A (2010). Beta Regression in R. \emph{Journal of Statistical
#'     Software}, \bold{34}, 1-24
#' @references Zeileis A, Croissant Y (2010). Extended Model Formulas in R: Multiple Parts and Multiple
#'     Responses. \emph{Journal of Statistical Software}, \bold{34}, 1-13.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
NULL

# Main function ------------------------------------------------------------------------------------
#' @rdname ugrpl
#' @export
ugrpl <- function(formula, data, subset, na.action, link = "aordaz", sigma.link,
                  control = ug_control(...), y = TRUE, x = TRUE,...)
{

  ret.y <- y
  ret.x <- x

  ## Call
  cl <- match.call()
  if (missing(data))  data <- environment(formula)

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  mf <- mf[c(1L, m)]
  mf$drop.unused.levels <- TRUE

  ## Formula
  formula <- Formula::as.Formula(formula)
  if (length(formula)[2L] < 2L) {
    formula <- Formula::as.Formula(formula(formula), ~ 1)
  } else {
    if (length(formula)[2L] > 2L) {
      formula <- Formula::Formula(formula(formula, rhs = 1:2))
      warning("formula must not have more than two RHS parts")
    }
  }

  mf$formula <- formula

  ## Evaluate model.frame
  mf[[1L]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())

  ## Extract terms, model matrix, response
  mt <- stats::terms(formula, data = data)
  mtX <- stats::terms(formula, data = data, rhs = 1L)
  mtZ <- stats::delete.response(stats::terms(formula, data = data, rhs = 2L))

  y <- stats::model.response(mf, "numeric")
  X <- stats::model.matrix(mtX, mf)
  Z <- stats::model.matrix(mtZ, mf)

  ## Some conditions
  if (length(y) < 1)
    stop("empty model")
  if (any(y < 0 & y > 1))
    stop("invalid dependent variable, all observations must be in (0, 1)")

  ## Useful variables
  g <- make.plink
  k <- NCOL(X)
  l <- NCOL(Z)
  n <- length(y)

  ## Link functions
  if (is.null(sigma.link))  sigma.link <- if (l == 1) "identity" else "aordaz"

  ## Fit
  opt <- ug_mle(y, X, Z, link = link, sigma.link = sigma.link, control = control)

  ## Convergence status
  convergence <- opt$convergence == 0

  ## Starting values
  start <- opt$start
  opt$start <- NULL

  ## Coefficients
  beta <- opt$par[1:k]
  names(beta) <- colnames(X)
  gamma <- opt$par[1:l + k]
  names(gamma) <- colnames(Z)

  ## Lambdas
  lambda1 <- NULL
  lambda2 <- NULL
  nlambdas <- length(opt$par) - k - l
  if (g(link)$plink == TRUE) lambda1 <-  as.numeric(opt$par[k + l + 1])
  if (g(sigma.link)$plink == TRUE) lambda2 <- as.numeric(opt$par[k + l + nlambdas])

  ## Covariance matrix
  if(control$hessian){
    vcov <- chol2inv(Rfast::cholesky(-opt$hessian))
  }else{
    vcov <- K_ug(par = c(beta, gamma, lambda1, lambda2), X = X, Z = Z,
                 link = link, sigma.link = sigma.link, inverse = TRUE)
  }

  ## Fitted values
  mu  <- c(g(link)$linkinv(X%*%beta, lambda1))
  sigma <- c(g(sigma.link)$linkinv(Z%*%gamma, lambda2))

  ## Log-likelihood
  logLik <- opt$value

  ## set up return value
  out <- list(coefficients = list(mean = beta, dispersion = gamma),
              lambda = list(lambda1 = lambda1, lambda2 = lambda2),
              sigma = sigma,
              link = link,
              sigma.link = sigma.link,
              fitted.values = structure(mu, .Names = names(y)),
              logLik = logLik,
              vcov = vcov,
              residuals = as.numeric(stats::qnorm(pugamma(y, mu, sigma))),
              convergence = convergence,
              start = start,
              optim = opt,
              control = control,
              nobs = n,
              df.null = n - 2,
              df.residual = n - k - l - nlambdas)


  ## Further model information
  out$call <- cl
  out$formula <- formula
  out$terms <- list(mean = mtX, dispersion = mtZ, full = mt)
  if(ret.y) out$y <- y
  if(ret.x) out$x <- list(mean = X, dispersion = Z)

  class(out) <- "ugrpl"
  out

}
