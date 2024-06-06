#' @name ugrpl-methods
#' @title Methods for 'ugrpl' objects
#'
#' @param x,object an object of class \code{ugrpl}.
#' @param k numeric, the penalty per parameter to be used; the default
#'     \code{k = 2} is the classical AIC.
#' @param formula a model formula or terms object or an R object.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL


## Model frame
#' @export
#' @rdname ugrpl-methods
model.frame.ugrpl <- function(formula, ...) {
  formula$terms <- formula$terms$full
  formula$call$formula <- formula$formula <- formula(formula$terms)
  NextMethod()
}

## Model matrix
#' @export
#' @rdname ugrpl-methods
model.matrix.ugrpl <- function(object,
                               model = c("mean", "dispersion"), ...) {
  model <- match.arg(model, c("mean", "dispersion"))
  rval <- if(!is.null(object$x[[model]])) object$x[[model]]
  else stats::model.matrix(object$terms[[model]], stats::model.frame(object))
  return(rval)
}

# Parameter estimates
#' @rdname ugrpl-methods
#' @export
#' @param model a character indicating which parameter coefficients are
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"dispersion"} model. If \code{"full"} (default), a list with
#'   coefficients for the \code{mean} and for the \code{dispersion}
#'   model is returned.
coef.ugrpl <- function(object,
                       model = c("full", "mean", "dispersion"), ...) {

  model <- match.arg(model)

  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  Z <- stats::model.matrix(object$terms$dispersion, stats::model.frame(object))
  k <- ncol(X)
  l <- ncol(Z)

  sigma.link <- object$sigma.link

  ## Names
  mean_names <- colnames(X)
  if (l > 1) {
    disp_names <- colnames(Z)
  }else{

    if (sigma.link == "identity"){
      disp_names <- "sigma"
    }else{
      if (make.plink(sigma.link)$plink == TRUE){
        disp_names <- "g(sigma, lambda)"
      }else{
        disp_names <- "g(sigma)"
      }

    }
  }

  mu <- object$coef$mean
  names(mu) <- mean_names

  sigma <- object$coef$disp
  names(sigma) <- disp_names

  switch(model,
         "full"       = list(mean = mu, dispersion = sigma),
         "mean"       = mu,
         "dispersion" = sigma)
}

#  Variance-covariance matrix
#' @rdname ugrpl-methods
#' @export
vcov.ugrpl <- function(object,
                       model = c("full", "mean", "dispersion", "lambda"), ...) {

  model <- match.arg(model)
  covm <- object$vcov

  X <- stats::model.matrix(object$terms$mean, stats::model.frame(object))
  Z <- stats::model.matrix(object$terms$dispersion, stats::model.frame(object))
  k <- ncol(X)
  l <- ncol(Z)

  link <- object$link
  sigma.link <- object$sigma.link

  ## Names
  mean_names <- colnames(X)
  if (l > 1) {
    disp_names <- colnames(Z)
  }else{

    if (sigma.link == "identity"){
      disp_names <- "sigma"
    }else{
      if (make.plink(sigma.link)$plink == TRUE){
        disp_names <- "g(sigma, lambda)"
      }else{
        disp_names <- "g(sigma)"
      }

    }
  }

  switch(model,
         "mean" = {
           out <- covm[seq.int(length.out = k), seq.int(length.out = k), drop = FALSE]
           rownames(out) <- colnames(out) <- mean_names
           out
         },
         "dispersion" = {
           out <- covm[seq.int(length.out = l) + k,
                       seq.int(length.out = l) + k, drop = FALSE]
           rownames(out) <- colnames(out) <- disp_names
           out
         },
         "lambda" = {
           if (is.null(object$lambda$lambda1) &&
               is.null(object$lambda$lambda2)){
             print(paste("Parametric link functions were not used"))
           }else{
             nlambdas <- length(c(object$lambda$lambda1,
                                  object$lambda$lambda2))

             out <- object$vcov[(k + l + 1):(k + l + nlambdas),
                                (k + l + 1):(k + l + nlambdas)]
             rownames(out) <- colnames(out) <- c(if(make.plink(link)$plink) "lambda1" else NULL,
                                                 if(make.plink(sigma.link)$plink) "lambda2" else NULL)
             out
           }
         },
         "full" = {
           out <- covm
           rownames(out) <- colnames(out) <- c(mean_names, disp_names,
                                               if(make.plink(link)$plink) "lambda1" else NULL,
                                               if(make.plink(sigma.link)$plink) "lambda2" else NULL)
           out
         })

}


# Log-likelihood
#' @rdname ugrpl-methods
#' @export
logLik.ugrpl <- function(object, ...) {
  structure(object$logLik,
            df = object$nobs - object$df.residual,
            class = "logLik")
}


# AIC
#' @export
#' @rdname ugrpl-methods
AIC.ugrpl <- function(object, ..., k = 2) {

  AIC <- - 2 * object$logLik + k * length(object$optim$par)

  class(AIC) <- "AIC"
  return(AIC)
}

# Residuals
#' @name residuals.bergreg
#' @title Extract Model Residuals for a Unit Gamma Regression with Parametric Link Functions
#'
#' @param object a \code{'ugrpl'} object.
#' @param type character; it specifies which residual should be extracted.
#'     The available arguments are \code{"quantile"} (default), \code{"pearson"},
#'     \code{"response"} (raw residuals, y - mu), and \code{"deviance"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @export
#'
residuals.ugrpl <- function(object,
                            type = c("quantile", "pearson",
                                     "response", "deviance"), ...)
{
  ## raw response residuals and desired type
  res <- object$residuals
  type <- match.arg(type)

  ## extract fitted information
  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
  mu <- stats::fitted(object)
  sigma <- object$sigma

  if(type == "response") return(as.numeric(y - mu))

  res <- switch(type,
                "pearson" = {
                  phi <- 1 / sigma^2 - 1
                  as.numeric(res / sqrt(mu * ((1 / (2 - mu^(1 / phi))^phi) - mu)))
                },

                "deviance" = {
                  ll <- function(mu, sigma) dugamma(y, mu, sigma, log.p = TRUE)
                  as.numeric(sign(res) * sqrt(2 * abs(ll(y, sigma) - ll(mu, sigma))))
                },

                "quantile" = {
                  as.numeric(stats::qnorm(pugamma(y, mu, sigma)))
                })

  res
}

# Print
#' @rdname ugrpl-methods
#' @export
print.ugrpl <- function(x, ...)
{

  ## Covenience variables
  n <- x$nobs
  k <- length(x$coefficients$mean)
  l <- length(x$coefficients$dispersion)
  nlambdas <- length(c(x$lambda$lambda1, x$lambda$lambda2))

  ## Links
  link <- x$link
  sigma.link <- x$sigma.link

  ## Names
  mean_names <- colnames(stats::model.matrix(x$terms$mean, stats::model.frame(x)))
  if (l > 1) {
    disp_names <- colnames(stats::model.matrix(x$terms$dispersion,
                                               stats::model.frame(x)))
  }else{

    if (sigma.link == "identity"){
      disp_names <- "sigma"
    }else{
      if (make.plink(sigma.link)$plink == TRUE){
        disp_names <- "g(sigma, lambda)"
      }else{
        disp_names <- "g(sigma)"
      }

    }
  }

  cat("\nCall:\n")
  print(x$call)
  if(!x$convergence) {
    cat("\nmodel did not converge\n")
  }else{
    cat("\nMean Coefficients:\n")
    tab_mean <- (x$coefficients)$mean
    names(tab_mean) <- mean_names
    print(tab_mean)
    cat("\nDispersion Coefficients:\n")
    tab_disp <- (x$coefficients)$dispersion
    names(tab_disp) <- disp_names
    print(tab_disp)
  }

  invisible(x)
}

# Summary
#' @rdname ugrpl-methods
#' @export
summary.ugrpl <- function(object, ...)
{
  ## Covenience variables
  n <- object$nobs
  k <- length(object$coefficients$mean)
  l <- length(object$coefficients$dispersion)
  nlambdas <- length(c(object$lambda$lambda1,
                       object$lambda$lambda2))

  ## Links
  link <- object$link
  sigma.link <- object$sigma.link

  ## Names
  mean_names <- colnames(stats::model.matrix(object$terms$mean, stats::model.frame(object)))
  if (l > 1) {
    disp_names <- colnames(stats::model.matrix(object$terms$dispersion,
                                               stats::model.frame(object)))
  }else{

    if (sigma.link == "identity"){
      disp_names <- "sigma"
    }else{
      if (make.plink(sigma.link)$plink == TRUE){
        disp_names <- "g(sigma, lambda)"
      }else{
        disp_names <- "g(sigma)"
      }

    }
  }

  ## Summary for residuals
  res <- stats::residuals(object, type = "quantile")
  skewness <- mean( ((res - mean(res))/stats::sd(res))^3 )
  kurtosis <- mean( ((res - mean(res))/stats::sd(res))^4 )
  TAB.residuals <- round(cbind(mean(res), stats::sd(res),
                               skewness, kurtosis), 6)
  colnames(TAB.residuals) <- c("Mean", "Sd", "Skewness", "Kurtosis")
  rownames(TAB.residuals) <- " "

  # Summary for mu
  est.beta <- object$coef$mean
  se.beta <- sqrt(diag(object$vcov))[1:k]
  zval.beta <- est.beta/se.beta
  pval.beta <- 2 * stats::pnorm(abs(zval.beta), lower.tail = FALSE)

  TAB.mu <- cbind(Estimate = est.beta,
                  `Std. Error` = se.beta,
                  `z value` = zval.beta,
                  `Pr(>|z|)` = pval.beta)
  rownames(TAB.mu) <- mean_names

  # Summary for sigma
  est.gamma <- object$coef$disp
  se.gamma <- sqrt(diag(object$vcov)[1:l + k])
  zval.gamma <- est.gamma/se.gamma
  pval.gamma <- 2*stats::pnorm(abs(zval.gamma), lower.tail = FALSE)

  TAB.sigma <- cbind(Estimate = est.gamma,
                     `Std. Error` = se.gamma,
                     `z value` = zval.gamma,
                     `Pr(>|z|)` = pval.gamma)
  rownames(TAB.sigma) <- disp_names

  # Lambdas
  link <- object$link
  lambda1 <- object$lambda$lambda1
  TAB.lambda1 <- cbind(Estimate = lambda1,
                       `Std. Error` = sqrt(diag(object$vcov)[k + l + 1]))
  rownames(TAB.lambda1) <- "lambda1"

  sigma.link <- object$sigma.link
  lambda2 <- object$lambda$lambda2
  TAB.lambda2 <- cbind(Estimate = lambda2,
                       `Std. Error` = sqrt(diag(object$vcov)[k + l + nlambdas]))
  rownames(TAB.lambda2) <- "lambda2"

  logLik <- object$logLik

  # Pseudo-R2
  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
  fit0 <- try(ug_mle(y, X = matrix(1, n), link = "identity", sigma.link = "identity"), silent = TRUE)

  if (unique(grepl("Error", fit0))){
    r2G <- NULL
  } else{
    r2G <- 1 - exp(- 2 * (object$logLik - fit0$value) / n)
  }

  r2 <- stats::cor(make.plink(object$link)$linkfun(y, object$lambda$lambda1),
                   make.plink(object$link)$linkfun(object$fitted.values, object$lambda$lambda1))^2

  # Upsilon statistics
  Upsilon <- mean(abs(sort(res) - EnvStats::evNormOrdStats(n = n)))

  out <- list(call = object$call,
              residuals = TAB.residuals,
              link = object$link,
              sigma.link = object$sigma.link,
              mean = TAB.mu,
              dispersion = TAB.sigma,
              lambda1 = lambda1,
              lambda2 = lambda2,
              summ_l1 = TAB.lambda1,
              summ_l2 = TAB.lambda2,
              logLik = logLik,
              GpseudoR2 = r2G,
              pseudoR2 = r2,
              Upsilon = Upsilon,
              AIC = stats::AIC(object),
              BIC = stats::AIC(object, k = log(n)))

  class(out) <- "summary.ugrpl"
  out
}

# Print summary
#' @rdname ugrpl-methods
#' @export
print.summary.ugrpl <- function(x, ...)
{

  Link_mean <- make.plink(x$link)$name
  Link_disp <- make.plink(x$sigma.link)$name

  cat("Call:\n")
  print(x$call)

  cat("\nSummary for quantile residuals:\n")
  print(x$residuals)

  cat("\nMean model with", Link_mean, "link function:\n")

  if (!is.null(x$lambda1)){
    cat("\nLink parameter:\n")
    print(x$summ_l1)
  }

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$mean)

  cat("\nDispersion model with", Link_disp, "link function:\n")

  if (!is.null(x$lambda2)){
    cat("\nLink parameter:\n")
    print(x$summ_l2)
  }

  cat("\nCoefficients:\n")
  stats::printCoefmat(x$dispersion)

  cat(
    "\nLog-likelihood:", x$logLik,
    "\nUpsilon statistic:", round(x$Upsilon, 4), "(reference value: ~ 0)")

  if (!is.null(x$GpseudoR2)) {
    cat("\nGeneralized pseudo-R2:", round(x$GpseudoR2, 4), "(reference value: ~ 1)")
  } else {
    cat("\nPseudo-R2:", round(x$pseudoR2, 4), "(reference value: ~ 1)")
  }

  cat("\nAIC:", x$AIC, "and BIC:", x$BIC)

  invisible(x)
}

#  Predict
#' Predict Method for Unit Gamma Regression with Parametric Link Functions
#'
#' Obtains predictions from a fitted unit gamma regression object.
#'
#' @param object a \code{'ugrpl'} object.
#' @param newdata optionally, a data frame in which to look for variables
#'     with which to predict. If omitted, the fitted linear predictors are
#'     used.
#' @param type the type of prediction required. The default is on the scale of
#'     the response variable \code{("response")}, that is, the fitted values
#'     (fitted means). The alternative \code{"link"} is on the scale of
#'     the linear predictors. The
#'     specification \code{"dispersion"} provides the fitted dispersion parameters,
#'     while \code{"variance"} provides the fitted variances. Finally,
#'     the option \code{"quantile"} gives the fitted quantiles in the order
#'     specified via \code{'at'}.
#' @param at the order of the quantile to be predicted if
#'     \code{type = "quantile"}. The default is to predict the median,
#'     that is, \code{at = 0.5}.
#' @param na.action function determining what should be done with missing
#'     values in \code{newdata}. The default is to predict \code{NA}.
#' @param ...  arguments passed to or from other methods.
#'
#' @return A vector of predictions.
#' @export
predict.ugrpl <- function(object, newdata = NULL,
                          type = c("response", "link", "dispersion",
                                   "variance", "quantile"),
                          na.action = stats::na.pass, at = 0.5, ...)
{
  type <- match.arg(type)

  qfun <- function(at, mu, sigma) {
    rval <- sapply(at, function(p) qugamma(rep(p, length(mu)), mu, sigma))
    if(length(at) > 1L) {
      if(NCOL(rval) == 1L)
        rval <- matrix(rval, ncol = length(at),
                       dimnames = list(unique(names(rval)), NULL))

      colnames(rval) <- paste("q_", at, sep = "")
    } else {
      rval <- drop(rval)
    }
    rval
  }

  link <- object$link
  sigma.link <- object$sigma.link
  lambda1 <- object$lambda$lambda1
  lambda2 <- object$lambda$lambda2

  if(missing(newdata)) {

    mu <- as.numeric(object$fitted.values)
    sigma <- as.numeric(object$sigma)

    rval <- switch(type,
                   "response" = {
                     mu
                   },
                   "link" = {
                     make.plink(link)$linkfun(mu, lambda1)
                   },
                   "dispersion" = {
                     sigma
                   },

                   "variance" = {
                     phi <- 1 / sigma^2 - 1
                     mu * ((1 / (2 - mu^(1 / phi))^phi) - mu)
                   },
                   "quantile" = {
                     qfun(at, mu, sigma)
                   })

    return(rval)

  } else {

    tnam <- switch(type,
                   "response" = "mean",
                   "link" = "mean",
                   "dispersion" = "dispersion",
                   "variance" = "full",
                   "quantile" = "full")

    mf <- stats::model.frame(stats::delete.response(object$terms[[tnam]]),
                             newdata, na.action = na.action)
    newdata <- newdata[rownames(mf), , drop = FALSE]

    if(type %in% c("response", "link", "variance", "quantile")){
      X <- stats::model.matrix(stats::delete.response(object$terms$mean), mf)
      mu <- make.plink(link)$linkinv(drop(X %*% object$coefficients$mean), lambda1)
    }


    if(type %in% c("dispersion", "variance", "quantile")){
      Z <- stats::model.matrix(object$terms$dispersion, mf)
      sigma <- make.plink(sigma.link)$linkinv(drop(Z %*% object$coefficients$dispersion), lambda2)
    }

    rval <- switch(type,
                   "response" = {
                     mu
                   },
                   "link" = {
                     make.plink(link)$linkfun(mu, lambda1)
                   },
                   "dispersion" = {
                     sigma
                   },
                   "variance" = {
                     phi <- 1 / sigma^2 - 1
                     mu * ((1 / (2 - mu^(1 / phi))^phi) - mu)

                   },
                   "quantile" = {
                     qfun(at, mu, sigma)
                   }
    )
    return(rval)

  }
}
