# Influence ----------------------------------------------------------------------------------------

#' Leverage and Local Influence for the Unit Gamma Regression
#'
#'
#'
#' @param object an \code{"ugrpl"} object.
#' @param scheme character; it specifies the perturbation scheme. Currently, the following
#'     perturbation schemes are available: case-weight perturbation (\code{"case-weight"});
#'     perturbation of the mean covariates (\code{"mean"}); perturbation of the dispersion
#'     covariates (\code{"dispersion"}); and simultaneous perturbation of the mean and dispersion
#'     covariates (\code{"both"}). If the mean and the dispersion submodels share the same
#'     covariate, then it is not possible to use \code{scheme = "mean"} or \code{scheme = "dispersion"}.
#' @param covariate character; it specifies the name of the continuous variable used in the covariate
#'     perturbation scheme. If \code{scheme = "both"}, then it must be a vector of size 2, where
#'     the first element is the name of the perturbed covariate associated with the mean and the
#'     second element is the name of the perturbed covariate associated with the dispersion.
#'     If \code{scheme = "mean"} or \code{scheme = "dispersion"}, then the mean and dispersion
#'     submodels cannot share the same perturbed covariate.
#' @param plot 	logical. If \code{TRUE}, the graphs of local influence and generalized leverage for each
#'      sampled observation are presented.
#' @param ask logical. If \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed \link{plot}.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{dmax}{local influence under the selected perturbation scheme.}
#'   \item{total.LI}{total local influence.}
#'   \item{scheme}{name of the considered perturbation scheme}
#'   \item{leverage}{the diagonal elements of the generalized leverage matrix.}
#'  }
#'
#' @details
#' Leverage points are observations that have a disproportionate weight in their fitted value.
#'     These points, in general, present a different behavior from the other points in relation
#'     to the values of the explanatory variables and can strongly influence the estimates of the
#'     regression coefficients. On the other hand, local influence analysis should be considered
#'     when there is an interest in investigating the sensitivity of postulated assumptions under
#'     small perturbations in the model or data.
#'
#' The leverage measure considered is based on the work of Wei et al. (1998) and local influence
#'     perturbation schemes are commonly considered schemes, see, for instance Cook (1986).
#'
#' @references
#'
#' Cook, R. D. (1986). Assessment of local influence. \emph{Journal of the Royal Statistical
#'     Society  B}, \bold{48}, 133--155.
#'
#' Wei, B. C., Hu, Y. Q., & Fung, W. K. (1998). Generalized leverage and its applications.
#'     \emph{Scandinavian Journal of Statistics}, \bold{25}, 25--37.
#'
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
influence <- function(object,
                      scheme = c("case-weight", "mean", "dispersion", "both"),
                      covariate,
                      plot = TRUE,
                      ask = prod(graphics::par("mfcol")) < length(which) &&
                        grDevices::dev.interactive(),...){

  scheme <- match.arg(scheme, c("case-weight", "mean", "dispersion", "both"))

  ## Model specifications
  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y
  X <- stats::model.matrix(object, model = "mean")
  Z <- stats::model.matrix(object, model = "dispersion")

  k <- NCOL(X)
  l <- NCOL(Z)
  n <- length(y)

  link <- object$link
  sigma.link <- object$sigma.link

  ## Fitted means and dispersions
  mu <- object$fitted.values
  sigma <- object$sigma

  ## Parameters
  beta <- object$coefficients$mean
  gamma <- object$coefficients$dispersion
  lambda1 <- object$lambda$lambda1
  lambda2 <- object$lambda$lambda2

  ## Covenience transformations
  mu_star <- mu^(sigma^2 / (1 - sigma^2))
  y_star <- 1 + sigma^2 * mu_star * log(y) / (1 - sigma^2)

  d <- mu_star / (1 - mu_star)
  mu_dag <- digamma((1 - sigma^2) / sigma^2) - log(d)
  y_dag <- log(-log(y))

  a <- 1 / (mu * (1 - mu_star)^2)
  b <- 2 * mu * log(mu) / (sigma * (1 - sigma^2))

  ## Matrices
  T1 <- diag(as.vector(make.plink(link)$mu.eta(X%*%beta, lambda1)))
  T1_prime <- diag(as.vector(make.plink(link)$mu2.eta2(X%*%beta, lambda1)))
  T2 <- diag(as.vector(make.plink(sigma.link)$mu.eta(Z%*%gamma, lambda2)))
  T2_prime <- diag(as.vector(make.plink(sigma.link)$mu2.eta2(Z%*%gamma, lambda2)))

  A <- diag(as.vector(a))
  B <- diag(as.vector(b))

  D_star <- diag(as.vector(y_star - mu_star))
  D_dag  <- diag(as.vector(y_dag - mu_dag))

  Sm3 <- diag(as.vector(1 / sigma^3))

  if(scheme == "case-weight"){

    ## Delta
    Delta <- rbind(t(X)%*%T1%*%A%*%D_star,
                   t(Z)%*%T2%*%(A%*%B%*%D_star - 2 * Sm3%*%D_dag))

    if (make.plink(link)$plink){
      rho1 <- make.plink(link)$rho(X%*%beta, lambda1)
      Delta_l1 <- t(rho1)%*%A%*%D_star
    }else{
      Delta_l1 <- NULL
    }

    if (make.plink(sigma.link)$plink){
      rho2 <- make.plink(sigma.link)$rho(Z%*%gamma, lambda2)
      Delta_l2 <- t(rho2)%*%(A%*%B%*%D_star - 2 * Sm3%*%D_dag)
    }else{
      Delta_l2 <- NULL
    }

    Delta <- rbind(Delta, Delta_l1, Delta_l2)

  } else if(scheme == "mean"){

    if(any(colnames(Z) == covariate)) warning(paste0("the mean and dispersion submodels shares the same covariate: ",
                                                     covariate, "."))

    r <- which(colnames(X) == covariate)
    sd_x <- stats::sd(X[, r])
    cr <- matrix(as.numeric(1:k == r))


    V_star <- diag(c( (a^2 * (1 - mu_star) / (1 - sigma^2)) * (mu_star * (1 - mu_star) -
                                                                 (1 - y_star)*(2*sigma^2 + mu_star - 1)) - a^2 * (1 - mu_star)^2))

    C_star <- diag(c(  -a * (2*(1 - y_star)/sigma + b * (sigma^2) *
                               (1 - y_star - y_star * mu_star + mu_star^2)/
                               (mu * (1-mu_star)))/(1 - sigma^2)   ))

    Dlmu <- diag(c(a * (y_star - mu_star)))

    Delta <- rbind(sd_x * (beta[r] * t(X)%*%(V_star%*%(T1^2) + Dlmu%*%T1_prime) +
                             cr%*%t(y_star - mu_star)%*%A%*%T1),
                   sd_x * beta[r] * t(Z)%*%C_star%*%T1%*%T2)

    if (make.plink(link)$plink){
      rho1 <- make.plink(link)$rho(X%*%beta, lambda1)
      rho1_eta <- make.plink(link)$rho.eta(X%*%beta, lambda1)

      Delta_l1 <- beta[r] * sd_x * (t(rho1)%*%V_star%*%T1 + t(rho1_eta)%*%A%*%D_star)
    }else{
      Delta_l1 <- NULL
    }

    if (make.plink(sigma.link)$plink){
      rho2 <- make.plink(sigma.link)$rho(Z%*%gamma, lambda2)
      rho2_eta <- make.plink(sigma.link)$rho.eta(Z%*%gamma, lambda2)

      Delta_l2 <- beta[r] * sd_x * (t(rho2)%*%C_star%*%T1)
    }else{
      Delta_l2 <- NULL
    }

    Delta <- rbind(Delta, Delta_l1, Delta_l2)

  } else if(scheme == "dispersion"){

    if(any(colnames(X) == covariate)) warning(paste0("the mean and dispersion submodels shares the same covariate: ",
                                                     covariate, "."))

    r <- which(colnames(Z) == covariate)
    sd_z <- stats::sd(Z[, r])
    cr <- matrix(as.numeric(1:l == r))

    phi <- as.vector(((1-(sigma^2))/(sigma^2)))
    Q_star <- diag(c(   6 * (y_dag - mu_dag) / (sigma^4) -
                          4 * trigamma(phi) / (sigma^6) -
                          (a^2) * (b^2) * (sigma^2) * (1 - y_star - mu_star + mu_star^2 -
                                                         mu_star^3 + y_star * (mu_star^2)) / (1 - sigma^2) -
                          a * b * (4 + 3 * mu_star * (sigma^2) - 3 * y_star * (sigma^2) -
                                     3 * mu_star - y_star)/(sigma * (1 - sigma^2)) ))
    Dlsigma <- diag(c(a*b*(y_star - mu_star) - 2*(y_dag - mu_dag)/(sigma^3)))

    Delta <- rbind(sd_z * gamma[r] * t(X)%*%C_star%*%T1%*%T2,
                   sd_z * (gamma[r] * t(Z)%*%(Q_star%*%(T2^2) + Dlsigma%*%T2_prime) +
                             cr%*%t(A%*%B%*%(y_star - mu_star) - 2 * diag(Sm3) * (y_dag - mu_dag))%*%T2))

    if (make.plink(link)$plink){
      rho1 <- make.plink(link)$rho(X%*%beta, lambda1)
      rho1_eta <- make.plink(link)$rho.eta(X%*%beta, lambda1)

      Delta_l1 <- gamma[r] * sd_z * (t(rho1)%*%C_star%*%T2)
    }else{
      Delta_l1 <- NULL
    }

    if (make.plink(sigma.link)$plink){
      rho2 <- make.plink(sigma.link)$rho(Z%*%gamma, lambda2)
      rho2_eta <- make.plink(sigma.link)$rho.eta(Z%*%gamma, lambda2)

      Delta_l2 <- gamma[r] * sd_z * (t(rho2)%*%Q_star%*%T2 +
                                       t(rho2_eta)%*%(A%*%B%*%D_star - 2 * Sm3%*%D_dag))
    }else{
      Delta_l2 <- NULL
    }

    Delta <- rbind(Delta, Delta_l1, Delta_l2)

  } else if(scheme == "both"){

    r <- which(colnames(X) == covariate[1])
    sd_x <- stats::sd(X[, r])
    cr <- matrix(as.numeric(1:k == r))

    s <- which(colnames(Z) == covariate[2])
    sd_z <- stats::sd(Z[, s])
    cs <- matrix(as.numeric(1:l == s))

    phi <- as.vector(((1-(sigma^2))/(sigma^2)))

    V_star <- diag(c( (a^2 * (1 - mu_star) / (1 - sigma^2)) * (mu_star * (1 - mu_star) -
                                                                 (1 - y_star)*(2*sigma^2 + mu_star - 1)) - a^2 * (1 - mu_star)^2))

    C_star <- diag(c(  -a * (2*(1 - y_star)/sigma + b * (sigma^2) *
                               (1 - y_star - y_star * mu_star + mu_star^2)/
                               (mu * (1-mu_star)))/(1 - sigma^2)   ))

    Q_star <- diag(c(   6 * (y_dag - mu_dag) / (sigma^4) -
                          4 * trigamma(phi) / (sigma^6) -
                          (a^2) * (b^2) * (sigma^2) * (1 - y_star - mu_star + mu_star^2 -
                                                         mu_star^3 + y_star * (mu_star^2)) / (1 - sigma^2) -
                          a * b * (4 + 3 * mu_star * (sigma^2) - 3 * y_star * (sigma^2) -
                                     3 * mu_star - y_star)/(sigma * (1 - sigma^2)) ))

    Dlmu <- diag(c(a * (y_star - mu_star)))
    Dlsigma <- diag(c(a*b*(y_star - mu_star) - 2*(y_dag - mu_dag)/(sigma^3)))


    Delta <- rbind(sd_x * (beta[r] * t(X)%*%(V_star%*%(T1^2) + Dlmu%*%T1_prime) +
                             cr%*%t(y_star - mu_star)%*%A%*%T1) +
                     gamma[s] * sd_z * t(X)%*%C_star%*%T1%*%T2,
                   sd_z * (gamma[s] * t(Z)%*%(Q_star%*%(T2^2) + Dlsigma%*%T2_prime) +
                             cs%*%t(A%*%B%*%(y_star - mu_star) - 2 * diag(Sm3) * (y_dag - mu_dag))%*%T2) +
                     beta[r] * sd_x * t(Z)%*%C_star%*%T1%*%T2)

    if (make.plink(link)$plink){
      rho1 <- make.plink(link)$rho(X%*%beta, lambda1)
      rho1_eta <- make.plink(link)$rho.eta(X%*%beta, lambda1)

      Delta_l1 <- beta[r] * sd_x * (t(rho1)%*%V_star%*%T1 + t(rho1_eta)%*%A%*%D_star) +
        gamma[s] * sd_z * t(rho1)%*%C_star%*%T2

    }else{
      Delta_l1 <- NULL
    }

    if (make.plink(sigma.link)$plink){
      rho2 <- make.plink(sigma.link)$rho(Z%*%gamma, lambda2)
      rho2_eta <- make.plink(sigma.link)$rho.eta(Z%*%gamma, lambda2)

      Delta_l2 <- gamma[r] * sd_z * (t(rho2)%*%Q_star%*%T2 +
                                       t(rho2_eta)%*%(A%*%B%*%D_star - 2 * Sm3%*%D_dag)) +
        beta[r] * sd_x * t(rho2)%*%C_star%*%T1

    }else{
      Delta_l2 <- NULL
    }

    Delta <- rbind(Delta, Delta_l1, Delta_l2)

  }

  ## Hessian matrix
  J <- J_ug(c(beta, gamma, lambda1, lambda2), y, X, Z, link, sigma.link)

  d_max <- abs(eigen(-t(Delta)%*%solve(-J)%*%Delta)$vec[,1])

  totalLI <- NULL
  for(i in 1:n){
    totalLI[i] <- 2 * abs(t(Delta[,i])%*%solve(-J)%*%Delta[,i])
  }

  totalLI <- abs((totalLI - mean(totalLI)) / stats::sd(totalLI))

  # generalized leverage
  D1 <- diag(as.vector(a * sigma^2 * mu_star / (y * (1 - sigma^2))))
  D2 <- diag(as.vector(a * b * sigma^2 * mu_star / (y * (1 - sigma^2)) - 2 / (sigma^3 * y * log(y))))

  Ltheta <- cbind(T1%*%X, matrix(0L, n, l))
  Ltheta_y <- rbind(t(X)%*%T1%*%D1, t(Z)%*%T2%*%D2)

  if (make.plink(link)$plink){
    rho1 <- make.plink(link)$rho(X%*%beta, lambda1)
    Ll1 <- matrix(rho1)
    Ll1_y <- t(rho1)%*%D1
  }else{
    Ll1 <- NULL
    Ll1_y <- NULL
  }

  if (make.plink(sigma.link)$plink){
    Ll2 <- matrix(0L, n, 1)
    Ll2_y <- t(rho2)%*%D2
  }else{
    Ll2 <- NULL
    Ll2_y <- NULL
  }

  Ltheta <- cbind(Ltheta, Ll1, Ll2)
  Ltheta_y <- rbind(Ltheta_y, Ll1_y, Ll2_y)

  GL <- diag(Ltheta%*%solve(-J)%*%Ltheta_y)

  # Plot

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  Scheme <- c("Case-weight perturbation", "Mean covariate perturbation",
              "Dispersion covariate perturbation",
              "Simultaneous mean and dispersion covariate perturbation")[scheme == c("case-weight", "mean", "dispersion", "both")]


  if(plot == TRUE){
    plot(d_max, type = "h", main = Scheme, las = 1, xlab = "Index", ylab = "Local influence", ...)
    plot(GL, type = "h", las = 1, cex = 0.8, main = "Generalized leverage", xlab = "Index", ylab = expression(GL[ii]), ...)
  }

  list(dmax = d_max, total.LI = totalLI, scheme = Scheme, leverage = GL)
}



# Plot ---------------------------------------------------------------------------------------------

# Plot
#' @name plot.ugrpl
#' @title Diagnostic plots for the Unit Gamma regression with parametric link
#'  function
#' @param x an object of class \code{"ugrpl"}.
#' @param which numeric. If a subset of the plots is required, specify a subset
#'     of the numbers \\code{1:8}. See details below.
#' @param type character; it specifies which residual should be produced
#'     in the summary plot. The available arguments are \code{"quantile"} (default),
#'      \code{"pearson"}, \code{"response"} (raw residuals, y - mu), and \code{"deviance"}.
#' @param ask logical. If \code{TRUE}, the user is asked before each plot.
#' @param ... further arguments passed to or from other methods.
#'
#' @details The \code{which} argument can be used to select a subset of currently eight supported
#'     types of displays. The displays are:
#'
#' \describe{
#'   \item{\code{which = 1}:}{residuals vs fitted values}
#'   \item{\code{which = 2}:}{residuals vs indices of observations}
#'   \item{\code{which = 3}:}{normal probability plot of the residuals. Recall that only the
#'       quantile residuals have a standard normal distribution.}
#'   \item{\code{which = 4}:}{local influence versus indices of observations under the case-weights
#'       perturbation scheme.}
#'   \item{\code{which = 5}:}{generalized leverage versus indices of observations.}
#'   \item{\code{which = 6}:}{response vs fitted values.}
#'   \item{\code{which = 7}:}{autocorrelation plot of the residuals.}
#'   \item{\code{which = 8}:}{partial autocorrelation plot of the residuals.}
#'  }
#'
#' @export
#'
plot.ugrpl <- function(x,
                       which = 1:4,
                       type = c("quantile", "pearson",
                                "response", "deviance"),
                       ask = prod(graphics::par("mfcol")) < length(which) &&
                         grDevices::dev.interactive(), ...)
{

  if(!is.numeric(which) || any(which < 1) || any(which > 8))
    stop("`which' must be in 1:8")

  ## Reading
  type <- match.arg(type, c("quantile", "pearson", "response", "deviance"))
  res <- stats::residuals(x, type = type)

  ## Legends
  types <- c("quantile", "pearson", "response", "deviance")
  Types <- c("Quantile residuals", "Pearson residuals",
             "Raw response residuals", "Deviance residuals")
  Type <- Types[type == types]

  ## Graphical parameters setting
  if (ask) {
    op <- graphics::par(ask = TRUE)
    on.exit(graphics::par(op))
  }

  ## Plots to shown
  show <- rep(FALSE, 8)
  show[which] <- TRUE

  ## Residuals versus Fitted values
  if (show[1]){
    graphics::plot(stats::fitted(x), res, xlab = "Fitted values", ylab = Type, pch = "+", ...)
    graphics::abline(h = 0, col = "gray50", lty = 2)
  }

  ## Residuals versus index observation
  if (show[2]){
    n <- x$nobs
    graphics::plot(1:n, res, xlab = "Index", ylab = Type, pch = "+", ...)
    graphics::abline(h = 0, col= "gray50", lty = 2)
  }

  ## Normal probability plot
  if(show[3]) {
    if (type == "quantile") {
      car::qqPlot(res, id = FALSE, pch = 3, grid = FALSE,
                  xlab = "Normal quantiles", ylab = Type)
    } else {
      stats::qqnorm(res, pch = "+", xlab = "Normal quantiles", ylab = Type, main = "", ...)
      graphics::abline(0, 1, col = "gray50", lty = 2)
    }

  }

  ## Local influence
  if(show[4]){
    inf <- influence(x, plot = FALSE)
    plot(inf$dmax, type = "h", main = paste(inf$scheme, "scheme"),
         las = 1, xlab = "Index", ylab = "Local influence", ...)
  }

  ## Generalized leverage
  if(show[5]){
    inf <- influence(x, plot = FALSE)
    plot(inf$leverage, type = "h", las = 1, cex = 0.8,
         main = "Generalized leverage", xlab = "Index", ylab = expression(GL[ii]), ...)
  }

  ## Fitted versus response
  if(show[6]) {
    y <- if(is.null(x$y)) stats::model.response(stats::model.frame(x)) else x$y
    graphics::plot(y, stats::fitted(x), pch = "+", xlab = "Observed values", ylab = "Predicted values", ...)
    graphics::abline(0, 1, lty = 2, col = "gray")
  }

  ## ACF of residuals
  if(show[7]) {
    stats::acf(res, main = " ", xlab = "Lags", ylab = paste("Sample ACF of", type, "residuals"))
  }

  ## PACF of residuals
  if(show[8]) {
    stats::pacf(res, main = " ", xlab = "Lags", ylab = paste("Sample PACF of", type, "residuals"))
  }

}

# RESET test ---------------------------------------------------------------------------------------
#' RESET-Type Goodness-of-Fit Test
#'
#' \code{reset} performs the RESET-type goodness-of-fit test to verify the correct specification
#'     of the unit gamma regression with parametric link function.
#'
#' @param object an object of class \code{"ugrpl"}.
#'
#' @return A list with the following components:
#' \describe{
#'   \item{statistic:}{the value the likelohood ratio test statistic.}
#'   \item{p.value:}{the p-value for the test.}
#'  }
#'
#' @details ...
#'
#' @export
#'
reset <- function(object){

  y <- if(is.null(object$y)) stats::model.response(stats::model.frame(object)) else object$y

  eta <- stats::predict(object, type = "link")


  X <- cbind(stats::model.matrix(object$terms$mean, stats::model.frame(object)), eta^2)
  Z <- cbind(stats::model.matrix(object$terms$dispersion, stats::model.frame(object)), eta^2)

  link <- object$link
  sigma.link <- object$sigma.link

  # Initial values for the new fit
  inits <- c(object$coefficients$mean, 0L, object$coefficients$dispersion, 0L)

  # Log-likelihood with fixed link parameters
  ll_aux <- function(par, y, X, Z, link = link, sigma.link = sigma.link)
  {

    k <- ncol(X); l <- ncol(Z); n <- length(y)

    beta <- as.vector(par[1:k])
    gamma <- as.vector(par[1:l + k])
    lambda1 <- NULL
    lambda2 <- NULL

    nlambdas <- make.plink(link)$plink + make.plink(sigma.link)$plink
    if (make.plink(link)$plink)
      lambda1 <- object$lambda$lambda1
    if (make.plink(sigma.link)$plink)
      lambda2 <- object$lambda$lambda2

    mu <- make.plink(link)$linkinv(X%*%beta, lambda1)
    sigma <- make.plink(sigma.link)$linkinv(Z%*%gamma, lambda2)

    if (any(mu <= 0) | any(mu >= 1) | any(sigma <= 0) | any(sigma >= 1)){
      NaN
    }else{
      ll <- suppressWarnings(dugamma(y, mu, sigma, log.p = TRUE))

      if (any(!is.finite(ll)))
        NaN
      else
        -sum(ll)
    }
  }

  opt <- stats::optim(inits, ll_aux, y = y, X = X, Z = Z, link = link, sigma.link = sigma.link)

  # Likelihood ratio statistics
  LR <- 2 * (-opt$value - object$logLik)

  # Critical value
  critical <- stats::qchisq(1 - 0.05, 2L)

  # p-Value
  p.value <- stats::pchisq(LR, 2L, lower.tail = FALSE)

  cat("\nREST-type likelihood ratio test\n",
      "\nDegrees of freedom:", 2L,
      "\nStatistic:", LR,
      "\np-value:", p.value)

}


