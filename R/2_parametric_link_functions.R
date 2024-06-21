#' Create Non-parametric and Parametric Link Functions
#'
#' @description This function has the same spirit as the \code{\link[stats]{make.link}}
#'      \{stats\} function which creates links for GLM families. However, it
#'      additionally allows the use of parametric link functions as, for
#'      example, the Aranda-Ordaz link function.
#'
#' @param link character; see details to view the current available link
#'  functions.
#'
#' @return Given a link abbreviation, it returns a link function
#'      (\code{linkfun}); an inverse link function (\code{linkinv}); the first and the second
#'       derivatives of \code{mu} with respect to \code{eta}  (\code{mu.eta} and \code{mu2.eta2},
#'       respectively); the first and the second derivatives of \code{mu} with respect to
#'       \code{lambda} (\code{rho} and \code{rho2}, respectively); a logical value which is
#'       \code{TRUE} if the link function belongs to a parametric family (\code{plink}); and
#'       the lowercase name of link function (\code{name}). More specifically, it returns a list
#'       with the following components:
#'
#' \describe{
#'   \item{linkfun}{Link function \code{function(mu, lambda)}.}
#'   \item{linkinv}{Inverse link function \code{function(eta, lambda)}.}
#'   \item{mu.eta}{Derivative \code{function(eta, lambda)} \emph{dmu/deta}.}
#'   \item{mu2.eta2}{Second order derivative \code{function(eta, lambda)} \emph{d2 mu/d eta2}.}
#'   \item{rho}{Derivative \code{function(eta, lambda)} \emph{dmu/dlambda}. If
#'   the link function does not belongs to a parametric family, then it returns
#'   \code{NULL}.}
#'   \item{rho2}{Second order derivative \code{function(eta, lambda)} \emph{d2 mu/d lambda2}. If
#'   the link function does not belongs to a parametric family, then it returns
#'   \code{NULL}.}
#'   \item{rho.eta}{Second order derivative \code{function(eta, lambda)} \emph{d2 mu/dlambda deta}. If
#'   the link function does not belongs to a parametric family, then it returns
#'   \code{NULL}.}
#'   \item{plink}{logical; if \code{TRUE}, the link function belongs to a parametric
#'    family.}
#'   \item{name}{A name to be used for the link.}
#'  }
#'
#' @details We assume that the link functions belonging to a parametric
#'     family are indexed by a positive parameter \code{lambda}.
#'     When the link function does not belong to this family
#'     (e.g., the logit function), then, by default, \code{lambda = NULL}.
#'     Otherwise, \code{lambda} must be specified.
#'
#'     The available link functions are
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
#'  Power cauchit \tab \code{"pcauchit"} \tab \code{TRUE}  \cr
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
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
#' @examples
#' ### Non-parametric link functions
#' curve(make.plink("logit")$linkinv(x), xlim = c(-10, 10), ylim = c(0, 1),
#'       main = "Non-parametric link functions",
#'       xlab = expression(eta), ylab = "Inverse")
#' curve(make.plink("probit")$linkinv(x), add = TRUE, col = 2)
#' curve(make.plink("cauchit")$linkinv(x), add = TRUE, col = 3)
#' curve(make.plink("loglog")$linkinv(x), add = TRUE, col = 4)
#' curve(make.plink("cloglog")$linkinv(x), add = TRUE, col = 6)
#' legend("bottomright", legend = c("logit", "probit", "cauchit", "loglog", "cloglog"),
#'        lty = 1, col = c(1,2,3,4,6), bty = "n")
#'
#' ### Aranda-Ordaz
#' curve(make.plink("aordaz")$linkinv(x, 0.5), xlim = c(-10, 10), ylim = c(0, 1),
#'       main = "Aranda-Ordaz link function", xlab = expression(eta), ylab = "Inverse")
#' curve(make.plink("aordaz")$linkinv(x, 1), add = TRUE, col = 2)
#' curve(make.plink("aordaz")$linkinv(x, 3), add = TRUE, col = 3)
#' curve(make.plink("aordaz")$linkinv(x, 6), add = TRUE, col = 4)
#' legend("bottomright", legend = c(expression(lambda == 0.5),
#'                                  expression(lambda == 1),
#'                                  expression(lambda == 3),
#'                                  expression(lambda == 6)),
#'        lty = 1, col = c(1,2,3,4), bty = "n")
#'
#' ### Reversal Aranda-Ordaz
#' curve(make.plink("raordaz")$linkinv(x, 0.5), xlim = c(-10, 10), ylim = c(0, 1),
#'       main = "Reversal Aranda-Ordaz link function", xlab = expression(eta), ylab = "Inverse")
#' curve(make.plink("raordaz")$linkinv(x, 1), add = TRUE, col = 2)
#' curve(make.plink("raordaz")$linkinv(x, 3), add = TRUE, col = 3)
#' curve(make.plink("raordaz")$linkinv(x, 6), add = TRUE, col = 4)
#' legend("bottomright", legend = c(expression(lambda == 0.5),
#'                                  expression(lambda == 1),
#'                                  expression(lambda == 3),
#'                                  expression(lambda == 6)),
#'        lty = 1, col = c(1,2,3,4), bty = "n")
#'
#' ### Power logit
#' curve(make.plink("plogit")$linkinv(x, 0.5), xlim = c(-10, 10), ylim = c(0, 1),
#'       main = "Power logit link function", xlab = expression(eta), ylab = "Inverse")
#' curve(make.plink("plogit")$linkinv(x, 1), add = TRUE, col = 2)
#' curve(make.plink("plogit")$linkinv(x, 3), add = TRUE, col = 3)
#' curve(make.plink("plogit")$linkinv(x, 6), add = TRUE, col = 4)
#' legend("bottomright", legend = c(expression(lambda == 0.5),
#'                                  expression(lambda == 1),
#'                                  expression(lambda == 3),
#'                                  expression(lambda == 6)),
#'        lty = 1, col = c(1,2,3,4), bty = "n")
#'
#'
#' ### Reversal power logit
#' curve(make.plink("rplogit")$linkinv(x, 0.5), xlim = c(-10, 10), ylim = c(0, 1),
#'       main = "Reversal power logit link function", xlab = expression(eta), ylab = "Inverse")
#' curve(make.plink("rplogit")$linkinv(x, 1), add = TRUE, col = 2)
#' curve(make.plink("rplogit")$linkinv(x, 3), add = TRUE, col = 3)
#' curve(make.plink("rplogit")$linkinv(x, 6), add = TRUE, col = 4)
#' legend("bottomright", legend = c(expression(lambda == 0.5),
#'                                  expression(lambda == 1),
#'                                  expression(lambda == 3),
#'                                  expression(lambda == 6)),
#'        lty = 1, col = c(1,2,3,4), bty = "n")
#'
#' @references
#' Bazán, J. L., Torres-Avilés, F., Suzuki, A. K., & Louzada, F. (2017). Power and reversal power
#'  links for binary regressions: An application for motor insurance policyholders.
#'   \emph{Applied Stochastic Models in Business and Industry}, \bold{33}(1), 22--34.
#'
#'
#'

 make.plink <- function(link)
 {
   switch(link,

          # Logit ------------------------------------------------------------------------------------------
          logit = {

            linkfun <- function(mu, lambda = NULL, ...) stats::qlogis(mu)

            linkinv <- function(eta, lambda = NULL) stats::plogis(eta)

            mu.eta <- function(eta, lambda = NULL) stats::dlogis(eta)

            mu2.eta2 <- function(eta, lambda = NULL) - stats::dlogis(eta) * (1 - exp(-eta)) / (1 + exp(-eta))

            rho <- function(eta, lambda = NULL) NULL

            rho2 <- function(eta, lambda = NULL) NULL

            rho.eta <- function(eta, lambda = NULL) NULL

            plink <- FALSE

            name <- "logit"

          },

          # Probit -----------------------------------------------------------------------------------------
          probit = {

            linkfun <- function(mu, lambda = NULL, ...) stats::qnorm(mu)

            linkinv <- function(eta, lambda = NULL) stats::pnorm(eta)

            mu.eta <- function(eta, lambda = NULL) stats::dnorm(eta)

            mu2.eta2 <- function(eta, lambda = NULL) -eta * stats::dnorm(eta)

            rho <- function(eta, lambda = NULL) NULL

            rho2 <- function(eta, lambda = NULL) NULL

            rho.eta <- function(eta, lambda = NULL) NULL

            plink <- FALSE

            name <- "probit"
          },

          # Cauchit ----------------------------------------------------------------------------------------
          cauchit = {

            linkfun <- function(mu, lambda = NULL, ...) stats::qcauchy(mu)

            linkinv <- function(eta, lambda = NULL) stats::pcauchy(eta)

            mu.eta <- function(eta, lambda = NULL) stats::dcauchy(eta)

            mu2.eta2 <- function(eta, lambda = NULL) - 2 * eta * stats::dcauchy(eta) / (1 + eta^2)

            rho <- function(eta, lambda = NULL) NULL

            rho2 <- function(eta, lambda = NULL) NULL

            rho.eta <- function(eta, lambda = NULL) NULL

            plink <- FALSE

            name <- "cauchit"
          },

          # Log-log ----------------------------------------------------------------------------------------
          loglog = {

            linkfun <- function(mu, lambda = NULL, ...) -log(-log(mu))

            linkinv <- function(eta, lambda = NULL) pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)

            mu.eta <- function(eta, lambda = NULL) exp(-eta - exp(-eta))

            mu2.eta2 <- function(eta, lambda = NULL) exp(-2 * eta - exp(-eta)) * (1 - exp(eta))

            rho <- function(eta, lambda = NULL) NULL

            rho2 <- function(eta, lambda = NULL) NULL

            rho.eta <- function(eta, lambda = NULL) NULL

            plink <- FALSE

            name <- "complement log-log"
          },

          # Complement log-log -----------------------------------------------------------------------------
          cloglog = {

            linkfun <- function(mu, lambda = NULL, ...) log(-log(1 - mu))

            linkinv <- function(eta, lambda = NULL) pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)

            mu.eta <- function(eta, lambda = NULL) exp(eta - exp(eta))

            mu2.eta2 <- function(eta, lambda = NULL)  exp(eta - exp(eta)) * (1 - exp(eta))

            rho <- function(eta, lambda = NULL) NULL

            rho2 <- function(eta, lambda = NULL) NULL

            rho.eta <- function(eta, lambda = NULL) NULL

            plink <- FALSE

            name <- "complement log-log"
          },

          # Identity ---------------------------------------------------------------------------------------
          identity = {

            linkfun <- function(mu, lambda = NULL, ...) mu

            linkinv <- function(eta, lambda = NULL) eta

            mu.eta <- function(eta, lambda = NULL) rep.int(1, length(eta))

            mu2.eta2 <- function(eta, lambda = NULL) rep.int(0, length(eta))

            rho <- function(eta, lambda = NULL) NULL

            rho2 <- function(eta, lambda = NULL) NULL

            rho.eta <- function(eta, lambda = NULL) NULL

            plink <- FALSE

            name <- "identity"
          },

          # Aranda-Ordaz ------------------------------------------------------------------------
          aordaz = {

            linkfun <- function(mu, lambda, ...) log( ((1 - mu)^(-lambda) - 1) / lambda )


            linkinv <- function(eta, lambda) 1 - (1 + lambda * exp(eta))^(-1 / lambda)

            mu.eta <- function(eta, lambda) exp(eta) * (1 + lambda * exp(eta))^(-1 - 1 / lambda)

            mu2.eta2 <- function(eta, lambda) exp(eta) * (1 - (1 + lambda) / (exp(-eta) + lambda)) /
              ((1 + lambda * exp(eta))^(1 + 1 / lambda))

            rho <- function(eta, lambda){

              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
              e <- pmin(pmax(exp(eta), -.Machine$double.xmax), .Machine$double.xmax)

              (1 / (exp(-eta) + lambda) - log1p(lambda * e) / lambda) /
                (lambda * (1 + lambda * e)^(1/lambda))
            }

            rho2 <- function(eta, lambda){

              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
              e <- pmin(pmax(exp(eta), -.Machine$double.xmax), .Machine$double.xmax)

              (2 * log(e * lambda + 1) / (lambda^3)  -
                  2 * e / (lambda^2 * (e * lambda + 1)) -
                  exp(2 * eta) / (lambda * (e * lambda + 1)^2) -
                  (log(e * lambda + 1) / (lambda^2) - e / (lambda * (e * lambda + 1)))^2) /
                ((e * lambda + 1)^(1 / lambda))
            }

            rho.eta <- function(eta, lambda){

              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
              e <- pmin(pmax(exp(eta), -.Machine$double.xmax), .Machine$double.xmax)

              e * (log1p(e * lambda) / lambda - (1 + lambda) / (exp(-eta) + lambda)) /
                (lambda * (1 + e * lambda)^(1 / lambda + 1))

            }

            plink <- TRUE

            name <- "Aranda-Ordaz"

          },

          # Power link functions --------------------------------------------------------------------

          ## Power logit ------------------------------------------------------------------------
          plogit = {

            linkfun <- function(mu, lambda, ...) stats::qlogis(mu^(1 / lambda))

            linkinv <- function(eta, lambda) stats::plogis(eta)^lambda

            mu.eta <- function(eta, lambda){

              g <- stats::plogis(eta)
              g. <- stats::dlogis(eta)
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- stats::plogis(eta)
              g. <- stats::dlogis(eta)
              g.. <- g. * (1 - exp(-eta)) / (1 + exp(-eta))
              lambda * ((lambda - 1) * g^(lambda - 2) * g.^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- stats::plogis(eta)
              g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- stats::plogis(eta)
              g^lambda  * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- stats::plogis(eta)
              g. <- stats::dlogis(eta)
              g. * g^(lambda  - 1) * (lambda  * log(g) + 1)

            }

            plink <- TRUE

            name <- "power logit"
          },

          # Power probit -----------------------------------------------------------------------------------
          pprobit = {

            linkfun <- function(mu, lambda, ...) stats::qnorm(mu^(1/lambda))

            linkinv <- function(eta, lambda) stats::pnorm(eta)^lambda

            mu.eta <- function(eta, lambda){

              g <- stats::pnorm(eta)
              g. <- stats::dnorm(eta)

              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- stats::pnorm(eta)
              g. <- stats::dnorm(eta)
              g.. <- -eta * g.

              lambda * ((lambda - 1) * g^(lambda - 2) * g.^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- stats::pnorm(eta)

              g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- stats::pnorm(eta)

              g^lambda  * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- stats::pnorm(eta)
              g. <- stats::dnorm(eta)

              g. * g^(lambda  - 1) * (lambda  * log(g) + 1)

            }

            plink <- TRUE

            name <- "power probit"
          },

          # Power cauchit -----------------------------------------------------------------------------------
          pcauchit = {

            linkfun <- function(mu, lambda, ...) stats::qcauchy(mu^(1 / lambda))

            linkinv <- function(eta, lambda) stats::pcauchy(eta)^lambda

            mu.eta <- function(eta, lambda){

              g <- stats::pcauchy(eta)
              g. <- stats::dcauchy(eta)

              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- stats::pcauchy(eta)
              g. <- stats::dcauchy(eta)
              g.. <- - 2 * eta * g. / (1 + eta^2)

              lambda * ((lambda - 1) * g^(lambda - 2) * g.^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- stats::pcauchy(eta)

              g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- stats::pcauchy(eta)

              g^lambda  * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- stats::pcauchy(eta)
              g. <- stats::dcauchy(eta)

              g. * g^(lambda  - 1) * (lambda  * log(g) + 1)

            }

            plink <- TRUE

            name <- "power cauchit"
          },

          # Power log-log -----------------------------------------------------------------------------------
          ploglog = {

            linkfun <- function(mu, lambda, ...) -log(-log(mu^(1 / lambda)))

            linkinv <- function(eta, lambda) pmax(pmin(exp(-exp(-eta))^lambda, 1 - .Machine$double.eps), .Machine$double.eps)

            mu.eta <- function(eta, lambda){

              g <- exp(-exp(-eta))
              g. <- exp(-exp(-eta) - eta)

              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- exp(-exp(-eta))
              g. <- exp(-exp(-eta) - eta)
              g.. <- exp(-2 * eta - exp(-eta)) * (1 - exp(eta))

              lambda * ((lambda - 1) * g^(lambda - 2) * g.^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- exp(-exp(-eta))

              g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- exp(-exp(-eta))

              g^lambda  * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- exp(-exp(-eta))
              g. <- exp(-exp(-eta) - eta)

              g. * g^(lambda  - 1) * (lambda  * log(g) + 1)

            }

            plink <- TRUE

            name <- "power loglog"
          },


          # Power complement log-log -----------------------------------------------------------------------------------
          pcloglog = {

            linkfun <- function(mu, lambda, ...) log(-log(1 - mu^(1 / lambda)))

            linkinv <- function(eta, lambda) pmax(pmin((-expm1(-exp(eta)))^lambda, 1 - .Machine$double.eps), .Machine$double.eps)

            mu.eta <- function(eta, lambda){

              g <- 1 - exp(-exp(eta))
              g. <- exp(eta - exp(eta))

              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- 1 - exp(-exp(eta))
              g. <- exp(eta - exp(eta))
              g.. <- g. * (1 - exp(eta))

              lambda * ((lambda - 1) * g^(lambda - 2) * g.^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- 1 - exp(-exp(eta))

              g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- 1 - exp(-exp(eta))

              g^lambda  * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- 1 - exp(-exp(eta))
              g. <- exp(eta - exp(eta))

              g. * g^(lambda  - 1) * (lambda  * log(g) + 1)

            }

            plink <- TRUE

            name <- "power complement loglog"
          },

          # Reversal power link functions -----------------------------------------------------------

          ## Reversal power logit ----
          rplogit = {

            linkfun <- function(mu, lambda = NULL, ...) -stats::qlogis((1 - mu)^(1/lambda))

            linkinv <- function(eta, lambda) 1 - stats::plogis(-eta)^lambda

            mu.eta <- function(eta, lambda){

              g <- stats::plogis(-eta)
              g. <- stats::dlogis(-eta)
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- stats::plogis(-eta)
              g. <- stats::dlogis(-eta)
              g.. <- g. * (1 - exp(eta)) / (1 + exp(eta))

              lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- stats::plogis(-eta)
              - g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- stats::plogis(-eta)

              - g^lambda * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- stats::plogis(-eta)
              g. <- stats::dlogis(-eta)
              g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.

            }

            plink <- TRUE

            name <- "reversal power logit"
          },

          ## Reversal power probit --------------------------------------------------------------------------
          rpprobit = {

            linkfun <- function(mu, lambda, ...) -stats::qnorm((1 - mu)^(1/lambda))

            linkinv <- function(eta, lambda) {
              1 - stats::pnorm(-eta)^lambda
            }

            mu.eta <- function(eta, lambda){

              g <- stats::pnorm(-eta)
              g. <- stats::dnorm(-eta)
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- stats::pnorm(-eta)
              g. <- stats::dnorm(-eta)
              g.. <- eta * g.

              lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- stats::pnorm(-eta)
              - g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- stats::pnorm(-eta)
              - g^lambda * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- stats::pnorm(-eta)
              g. <- stats::dnorm(-eta)
              g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.

            }

            plink <- TRUE

            name <- "reversal power probit"
          },

          ## Reversal power cauchit ------------------------------------------------------------------
          rpcauchit = {

            linkfun <- function(mu, lambda, ...) -stats::qcauchy((1 - mu)^(1 / lambda))

            linkinv <- function(eta, lambda) {
              1 - stats::pcauchy(-eta)^lambda
            }

            mu.eta <- function(eta, lambda){

              g <- stats::pcauchy(-eta)
              g. <- stats::dcauchy(-eta)
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- stats::pcauchy(-eta)
              g. <- stats::dcauchy(-eta)
              g.. <- 2 * eta * g. / (1 + eta^2)

              lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- stats::pcauchy(-eta)
              - g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- stats::pcauchy(-eta)

              - g^lambda * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- stats::pcauchy(-eta)
              g. <- stats::dcauchy(-eta)
              g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.

            }

            plink <- TRUE

            name <- "reversal power cauchit"
          },

          ## Reversal power log-log ------------------------------------------------------------------
          rploglog = {

            linkfun <- function(mu, lambda, ...) log(-log((1 - mu)^(1 / lambda)))

            linkinv <- function(eta, lambda){
              1 - pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^lambda
            }

            mu.eta <- function(eta, lambda){

              g <- exp(-exp(eta))
              g. <- exp(-exp(eta) + eta)
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- exp(-exp(eta))
              g. <- exp(-exp(eta) + eta)
              g.. <- exp(2 * eta - exp(eta)) * (1 - exp(-eta))

              lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- exp(-exp(eta))
              - g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- exp(-exp(eta))
              - g^lambda * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- exp(-exp(eta))
              g. <- exp(-exp(eta) + eta)
              g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.

            }

            plink <- TRUE

            name <- "reversal power log-log"
          },

          ## Reversal power complement log-log -------------------------------------------------------
          rpcloglog = {

            linkfun <- function(mu, lambda, ...) -log(-log(1 - (1 - mu)^(1 / lambda)))

            linkinv <- function(eta, lambda){
              1 - pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)^lambda
            }

            mu.eta <- function(eta, lambda){

              g <- 1 - exp(-exp(-eta))
              g. <- exp(-eta - exp(-eta))
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- 1 - exp(-exp(-eta))
              g. <- exp(-eta - exp(-eta))
              g.. <- g. * (1 - exp(-eta))

              lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- 1 - exp(-exp(-eta))
              - g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- 1 - exp(-exp(-eta))
              - g^lambda * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- 1 - exp(-exp(-eta))
              g. <- exp(-eta - exp(-eta))
              g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.

            }

            plink <- TRUE

            name <- "reversal power complement log-log"
          },

          # Reversal Aranda-Ordaz -------------------------------------------------------------------
          raordaz = {

            linkfun <- function(mu, lambda, ...) log(lambda / (mu^(-lambda) - 1))

            linkinv <- function(eta, lambda){
              e <- pmax(exp(-eta), .Machine$double.ep)
              (1 + lambda * e)^(-1 / lambda)
            }

            mu.eta <- function(eta, lambda){

              g <- 1 - (1 + lambda * exp(-eta))^(-1 / lambda)
              g. <- exp(-eta) * (1 + lambda * exp(-eta))^(-1 - 1 / lambda)
              lambda * g^(lambda - 1) * g.

            }

            mu2.eta2 <- function(eta, lambda) {

              g <- 1 - (1 + lambda * exp(-eta))^(-1 / lambda)
              g. <- exp(-eta) * (1 + lambda * exp(-eta))^(-1 - 1 / lambda)
              g.. <- exp(-eta) * (1 - (1 + lambda) / (exp(eta) + lambda)) /
                ((1 + lambda * exp(-eta))^(1 + 1 / lambda))

              lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)

            }

            rho <- function(eta, lambda){

              g <- 1 - (1 + lambda * exp(-eta))^(-1 / lambda)
              - g^lambda * log(g)

            }

            rho2 <- function(eta, lambda){

              g <- 1 - (1 + lambda * exp(-eta))^(-1 / lambda)
              - g^lambda * (log(g))^2

            }

            rho.eta <- function(eta, lambda){

              g <- 1 - (1 + lambda * exp(-eta))^(-1 / lambda)
              g. <- exp(-eta) * (1 + lambda * exp(-eta))^(-1 - 1 / lambda)
              g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.

            }

            plink <- TRUE

            name <- "reversal Aranda-Ordaz"
          },




          stop(gettextf("%s link not recognised", sQuote(link)), domain = NA))
   environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(mu2.eta2) <-
     environment(rho) <- environment(rho2) <- environment(rho.eta) <- asNamespace("stats")
   structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, mu2.eta2 = mu2.eta2,
                  rho = rho, rho2 = rho2, rho.eta = rho.eta, plink = plink, name = name), class = "link-glm")

}



# power_link_functions <- function(eta, lambda, g, g., g..) {
#
#   # First order derivatives
#   d1_eta <- lambda * g^(lambda - 1) * g.
#   d1_lambda <- g^lambda * log(g)
#
#   # Second order derivatives
#   d2_eta2  <- lambda * ((lambda - 1) * g^(lambda - 2) * g.^2 + g^(lambda - 1) * g..)
#   d2_lambda2  <- g^lambda  * (log(g))^2
#   d2_eta_lambda  <- g. * g^(lambda  - 1) * (lambda  * log(g) + 1)
#
#   list(linkinv = g^lambda,
#        mu.eta = d1_eta,
#        mu2.eta2 = d2_eta2,
#        rho = d1_lambda,
#        rho2 = d2_lambda2,
#        rho.eta = d2_eta_lambda)
# }

# compute_derivatives <- function(eta, lambda, g, g., g..) {
#
#   # OBS!!!
#   # g = g(-eta)
#   # g. = g.(-eta)
#   # g.. = g..(-eta)
#
#   # Derivadas de primeira ordem
#   mu.eta <- lambda * g^(lambda - 1) * g.
#   rho <- - g^lambda * log(g)
#
#   # Derivadas de segunda ordem
#   mu2.eta2 <- lambda * ((lambda - 1) * g^(lambda - 2) * (g.)^2 + g^(lambda - 1) * g..)
#   rho2 <- - g^lambda * (log(g))^2
#   rho.eta <- g^(lambda - 1) * g. + lambda * (lambda - 1) * g^(lambda - 2) * log(g) * g.
#
#   # Criando lista com as derivadas
#   #   list(mu.eta = mu.eta,
#   #        mu2.eta2 = mu2.eta2,
#   #        rho = rho,
#   #        rho2 = rho2,
#   #        rho.eta = rho.eta)
#
#   return(derivatives)
# }





# make.plinkOLD <- function(link)
# {
#   switch(link,
#
#          # Logit ------------------------------------------------------------------------------------------
#          logit = {
#
#            linkfun <- function(mu, lambda = NULL, ...) stats::qlogis(mu)
#
#            linkinv <- function(eta, lambda = NULL) {
#              thresh <- -stats::qlogis(.Machine$double.eps)
#              eta <- pmin(pmax(eta, -thresh), thresh)
#              stats::plogis(eta)
#            }
#
#            mu.eta <- function(eta, lambda = NULL) pmax(stats::dlogis(eta), .Machine$double.eps)
#
#            mu2.eta2 <- function(eta, lambda = NULL) {
#              - exp(eta) * expm1(eta) / (1 + exp(eta))^3
#            }
#
#            rho <- function(eta, lambda = NULL) NULL
#
#            rho2 <- function(eta, lambda = NULL) NULL
#
#            rho.eta <- function(eta, lambda = NULL) NULL
#
#            plink <- FALSE
#
#            name <- "logit"
#
#          },
#
#          # Probit -----------------------------------------------------------------------------------------
#          probit = {
#            linkfun <- function(mu, lambda = NULL, ...) stats::qnorm(mu)
#
#            linkinv <- function(eta, lambda = NULL) {
#              thresh <- -stats::qnorm(.Machine$double.eps)
#              eta <- pmin(pmax(eta, -thresh), thresh)
#              stats::pnorm(eta)
#            }
#
#            mu.eta <- function(eta, lambda = NULL) pmax(stats::dnorm(eta), .Machine$double.eps)
#
#            mu2.eta2 <- function(eta, lambda = NULL){
#              thresh <- -stats::qnorm(.Machine$double.eps)
#              eta <- pmin(pmax(eta, -thresh), thresh)
#              -eta * pmax(stats::dnorm(eta), .Machine$double.eps)
#            }
#
#            rho <- function(eta, lambda = NULL) NULL
#
#            rho2 <- function(eta, lambda = NULL) NULL
#
#            rho.eta <- function(eta, lambda = NULL) NULL
#
#            plink <- FALSE
#
#            name <- "probit"
#          },
#
#          # Cauchit ----------------------------------------------------------------------------------------
#          cauchit = {
#
#            linkfun <- function(mu, lambda = NULL, ...) stats::qcauchy(mu)
#
#            linkinv <- function(eta, lambda = NULL) {
#              thresh <- -stats::qcauchy(.Machine$double.eps)
#              eta <- pmin(pmax(eta, -thresh), thresh)
#              stats::pcauchy(eta)
#            }
#
#            mu.eta <- function(eta, lambda = NULL) pmax(stats::dcauchy(eta), .Machine$double.eps)
#
#            mu2.eta2 <- function(eta, lambda = NULL){
#              -2 * eta * pmax(stats::dcauchy(eta), .Machine$double.eps) / (1 + eta^2)
#            }
#
#            rho <- function(eta, lambda = NULL) NULL
#
#            rho2 <- function(eta, lambda = NULL) NULL
#
#            rho.eta <- function(eta, lambda = NULL) NULL
#
#            plink <- FALSE
#
#            name <- "cauchit"
#          },
#
#          # Log-log ----------------------------------------------------------------------------------------
#          loglog = {
#
#            linkfun <- function(mu, lambda = NULL, ...) -log(-log(mu))
#
#            linkinv <- function(eta, lambda = NULL) pmax(pmin(exp(-exp(-eta)),
#                                                              1 - .Machine$double.eps), .Machine$double.eps)
#            mu.eta <- function(eta, lambda = NULL) {
#              exp(-eta) * exp(-exp(-eta))
#            }
#
#            mu2.eta2 <- function(eta, lambda = NULL){
#              exp(-eta) * exp(-exp(-eta)) * expm1(-eta)
#            }
#
#            rho <- function(eta, lambda = NULL) NULL
#
#            rho2 <- function(eta, lambda = NULL) NULL
#
#            rho.eta <- function(eta, lambda = NULL) NULL
#
#            plink <- FALSE
#
#            name <- "complement log-log"
#          },
#
#          # Complement log-log -----------------------------------------------------------------------------
#          cloglog = {
#
#            linkfun <- function(mu, lambda = NULL, ...) log(-log(1 - mu))
#
#            linkinv <- function(eta, lambda = NULL) pmax(pmin(-expm1(-exp(eta)),
#                                                              1 - .Machine$double.eps), .Machine$double.eps)
#            mu.eta <- function(eta, lambda = NULL) {
#              exp(eta) * exp(-exp(eta))
#            }
#
#            mu2.eta2 <- function(eta, lambda = NULL){
#              eta <- pmin(eta, 700)
#              -exp(eta) * exp(-exp(eta)) * expm1(eta)
#            }
#
#            rho <- function(eta, lambda = NULL) NULL
#
#            rho2 <- function(eta, lambda = NULL) NULL
#
#            rho.eta <- function(eta, lambda = NULL) NULL
#
#            plink <- FALSE
#
#            name <- "complement log-log"
#          },
#
#          # Identity ---------------------------------------------------------------------------------------
#          identity = {
#
#            linkfun <- function(mu, lambda = NULL, ...) mu
#
#            linkinv <- function(eta, lambda = NULL) eta
#
#            mu.eta <- function(eta, lambda = NULL) rep.int(1, length(eta))
#
#            mu2.eta2 <- function(eta, lambda = NULL) rep.int(0, length(eta))
#
#            rho <- function(eta, lambda = NULL) NULL
#
#            rho2 <- function(eta, lambda = NULL) NULL
#
#            rho.eta <- function(eta, lambda = NULL) NULL
#
#            plink <- FALSE
#
#            name <- "identity"
#          },
#
#          # Asymmetric Aranda-Ordaz ------------------------------------------------------------------------
#          aordaz = {
#
#            linkfun <- function(mu, lambda, ...) {
#              log( ((1 - mu)^(-lambda) - 1) / lambda )
#            }
#
#            linkinv <- function(eta, lambda) {
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#
#              pmax(pmin((1 - (1 + lambda * exp(eta))^(-1 / lambda)),
#                        1 - .Machine$double.eps), .Machine$double.eps)
#            }
#
#            mu.eta <- function(eta, lambda) {
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              exp(eta) * (1 + lambda * exp(eta))^(-(1 + (1 / lambda)))
#
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#
#              exp(eta) * (1 + lambda * exp(eta))^(-(1 + (1/lambda))) *
#                (1 - (1 + lambda)/(exp(-eta) + lambda))
#            }
#
#            rho <- function(eta, lambda){
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              e <- pmin(pmax(exp(eta), -.Machine$double.xmax), .Machine$double.xmax)
#
#              (1 / (exp(-eta) + lambda) - log1p(lambda * e) / lambda) /
#                (lambda * (1 + lambda * e)^(1/lambda))
#            }
#
#            rho2 <- function(eta, lambda){
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              e <- pmin(pmax(exp(eta), -.Machine$double.xmax), .Machine$double.xmax)
#
#              (2 * log(e * lambda + 1) / (lambda^3)  -
#                  2 * e / (lambda^2 * (e * lambda + 1)) -
#                  exp(2 * eta) / (lambda * (e * lambda + 1)^2) -
#                  (log(e * lambda + 1) / (lambda^2) - e / (lambda * (e * lambda + 1)))^2) /
#                ((e * lambda + 1)^(1 / lambda))
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              e <- pmin(pmax(exp(eta), -.Machine$double.xmax), .Machine$double.xmax)
#
#              e * (log1p(e * lambda) / lambda - (1 + lambda) / (exp(-eta) + lambda)) /
#                (lambda * (1 + e * lambda)^(1 / lambda + 1))
#
#            }
#
#            plink <- TRUE
#
#            name <- "assymetric Aranda-Ordaz"
#
#          },
#
#
#          # Symmetric Aranda-Ordaz (parametric) ----
#          saordaz = {
#
#            linkfun <- function(mu, lambda, ...) {
#              (2/lambda)*(mu^lambda - (1 - mu)^lambda)/(mu^lambda + (1 - mu)^lambda)
#            }
#
#            linkinv <- function(eta, lambda) {
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#
#              etastar <- lambda * eta / 2
#
#              pmax(pmin(((1 + etastar)^(1 / lambda)) /
#                          (((1 - etastar)^(1 / lambda)) + ((1 + etastar)^(1 / lambda))),
#                        1 - .Machine$double.eps), .Machine$double.eps)
#
#            }
#
#            mu.eta <- function(eta, lambda) {
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              el <- lambda * eta
#
#              4 * ((4 - el^2)^(1 / lambda - 1)) / (((2 - el)^(1 / lambda) + (2 + el)^(1 / lambda))^2)
#
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              el <- lambda * eta
#
#              (8 * ((4 - el^2)^(1 / lambda - 1)) /
#                  (((2 - el)^(1 / lambda) + (2 + el)^(1 / lambda))^2)) *
#                (el * (lambda - 1) / (4 - el^2) -
#                   ((2 - el)^(1 / lambda - 1) + (2 + el)^(1 / lambda - 1)) /
#                   (((2 - el)^(1 / lambda) + (2 + el)^(1 / lambda))^3))
#            }
#
#            rho <- function(eta, lambda){
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              el <- lambda * eta
#
#              2 * ((4 - el^2)^(1 / lambda - 1)) *
#                ((el^2 - 4)*atanh(el/2) + 2 * el) /
#                ((lambda^2)*(((2 - el)^(1/lambda) + (2 + el)^(1/lambda))^2))
#            }
#
#            rho2 <- function(eta, lambda){
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              el <- lambda * eta
#
#              (2^(2 - 1/lambda))*((1 - 0.25*(el^2))^(1/lambda)) *
#                (4*eta*(lambda^2)*(((2-el)^(1/lambda))*(el^2 + eta - 2) +
#                                     ((2+el)^(1/lambda))*(el^2 - eta - 2)) +
#                   (el^2 - 4)*atanh(0.5*el)*
#                   (((2 - el)^(1/lambda))*((el^2 - 4)*atanh(0.5*el) +
#                                             lambda*(el^2 + 4*eta -4)) +
#                      ((2 + el)^(1/lambda))*(lambda*(el^2) + (4 - el^2)*atanh(0.5*el) -
#                                               4*lambda*(eta+1))))/
#                ((lambda^4)*((el^2 - 4)^2)*(((1 + 0.5*el)^(1/lambda) +
#                                               (1 - 0.5*el)^(1/lambda))^3))
#
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              eta <- pmin(pmax(eta, -.Machine$double.xmax), .Machine$double.xmax)
#              el <- lambda * eta
#
#              ret <- (2^(3 - 3/lambda)) * ((4 - el^2)^(1/lambda - 2)) *
#                (((2 + el)^(1/lambda))*((4 - el^2)*atanh(0.5*el) + el*(eta*(lambda^2) - 2)) +
#                   ((2 - el)^(1/lambda))*((el^2 - 4)*atanh(0.5*el) + el*(eta*(lambda^2) + 2)))/
#                ((lambda^2) * (((1 + 0.5*el)^(1/lambda) + (1 - 0.5*el)^(1/lambda))^3))
#            }
#
#            plink <- TRUE
#
#            name <- "symmetric Aranda-Ordaz"
#
#          },
#
#          # Power logit ------------------------------------------------------------------------------------
#          plogit = {
#
#            linkfun <- function(mu, lambda, ...) stats::qlogis(mu^(1/lambda))
#
#            linkinv <- function(eta, lambda = NULL) {
#              stats::plogis(eta)^lambda
#            }
#
#            mu.eta <- function(eta, lambda){
#              lambda * (stats::plogis(eta)^(lambda - 1)) * pmax(stats::dlogis(eta), .Machine$double.eps)
#            }
#
#            mu2.eta2 <- function(eta, lambda){
#
#              lambda * (stats::plogis(eta)^(lambda - 1)) *
#                (- exp(eta) * expm1(eta) / (1 + exp(eta))^3 + ((lambda - 1) / (stats::plogis(eta))) *
#                   (pmax(stats::dlogis(eta), .Machine$double.eps))^2)
#            }
#
#            rho <- function(eta, lambda){
#              (stats::plogis(eta)^lambda) * log(stats::plogis(eta))
#            }
#
#            rho2 <- function(eta, lambda){
#              (stats::plogis(eta)^lambda) * (log(stats::plogis(eta))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              (stats::plogis(eta)^(lambda-1)) * (1 + lambda * log(stats::plogis(eta))) *
#                pmax(stats::dlogis(eta), .Machine$double.eps)
#            }
#
#            plink <- TRUE
#
#            name <- "power logit"
#          },
#
#
#          # Power Probit -----------------------------------------------------------------------------------
#          pprobit = {
#
#            linkfun <- function(mu, lambda, ...) stats::qnorm(mu^(1/lambda))
#
#            linkinv <- function(eta, lambda) {
#              (stats::pnorm(eta))^lambda
#            }
#
#            mu.eta <- function(eta, lambda){
#              lambda * (stats::pnorm(eta)^(lambda - 1)) * pmax(stats::dnorm(eta), .Machine$double.eps)
#            }
#
#            mu2.eta2 <- function(eta, lambda){
#              lambda * (stats::pnorm(eta)^(lambda - 1)) *
#                (-eta * pmax(stats::dnorm(eta), .Machine$double.eps) +
#                   ((lambda - 1) / stats::pnorm(eta)) * (pmax(stats::dnorm(eta), .Machine$double.eps)^2))
#            }
#
#            rho <- function(eta, lambda){
#              (stats::pnorm(eta)^lambda) * log(stats::pnorm(eta))
#            }
#
#            rho2 <- function(eta, lambda){
#              (stats::pnorm(eta)^lambda) * (log(stats::pnorm(eta))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#              (stats::pnorm(eta)^(lambda-1)) *
#                pmax(stats::dnorm(eta), .Machine$double.eps) *
#                (1 + lambda*log(stats::pnorm(eta)))
#            }
#
#            plink <- TRUE
#
#            name <- "power probit"
#          },
#
#          # Power cauchit ----------------------------------------------------------------------------------
#          pcauchit = {
#
#            linkfun <- function(mu, lambda, ...) stats::qcauchy(mu^(1/lambda))
#
#            linkinv <- function(eta, lambda) {
#              (stats::pcauchy(eta)^lambda)
#            }
#
#            mu.eta <- function(eta, lambda){
#              lambda*(stats::pcauchy(eta)^(lambda-1))*
#                (pmax(stats::dcauchy(eta), .Machine$double.eps))
#            }
#
#            mu2.eta2 <- function(eta, lambda){
#              lambda * (stats::pcauchy(eta)^(lambda - 1)) *
#                (-2 * eta * pmax(stats::dcauchy(eta), .Machine$double.eps) / (1 + eta^2) +
#                   ((lambda - 1) / stats::pcauchy(eta)) *
#                   (pmax(stats::dcauchy(eta), .Machine$double.eps))^2)
#            }
#
#            rho <- function(eta, lambda){
#              (stats::pcauchy(eta)^lambda) * log(stats::pcauchy(eta))
#            }
#
#            rho2 <- function(eta, lambda){
#              (stats::pcauchy(eta)^lambda) * (log(stats::pcauchy(eta))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#              (stats::pcauchy(eta)^(lambda-1)) *
#                pmax(stats::dcauchy(eta), .Machine$double.eps) *
#                (1 + lambda*log(stats::pcauchy(eta)))
#            }
#
#            plink <- TRUE
#
#            name <- "power cauchit"
#          },
#
#
#          # Power log-log ----------------------------------------------------------------------------------
#          ploglog = {
#
#            linkfun <- function(mu, lambda, ...) -log(-log(mu^(1/lambda)))
#
#            linkinv <- function(eta, lambda){
#              pmax(pmin(exp(-exp(-eta))^lambda, 1 - .Machine$double.eps), .Machine$double.eps)
#            }
#
#            mu.eta <- function(eta, lambda) {
#
#              h_inv <- pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              lambda * (h_inv^(lambda - 1)) * pmax(exp(-eta) * exp(-exp(-eta)), .Machine$double.eps)
#
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#
#              h_inv <- pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              lambda * (h_inv^(lambda - 1)) * (exp(-eta) * exp(-exp(-eta)) * expm1(-eta) -
#                                                 ((1 - lambda) / h_inv) * (pmax(exp(-eta) * exp(-exp(-eta)), .Machine$double.eps)^2))
#            }
#
#            rho <- function(eta, lambda){
#              h_inv <- pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              (h_inv^lambda) * log(h_inv)
#            }
#
#            rho2 <- function(eta, lambda){
#              h_inv <- pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              (h_inv^lambda) * (log(h_inv)^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#              h_inv <- pmax(pmin(exp(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              (h_inv^(lambda - 1)) * pmax(exp(-eta) * exp(-exp(-eta)), .Machine$double.eps) *
#                (1 + lambda * log(h_inv))
#            }
#
#            plink <- TRUE
#
#            name <- "power log-log"
#          },
#
#
#          # Power complement log-log -----------------------------------------------------------------------
#          pcloglog = {
#
#            linkfun <- function(mu, lambda, ...) log(-log(1 - mu^(1/lambda)))
#
#            linkinv <- function(eta, lambda){
#              pmax(pmin(-expm1(-exp(eta)),
#                        1 - .Machine$double.eps),
#                   .Machine$double.eps)^lambda
#            }
#
#            mu.eta <- function(eta, lambda) {
#              lambda*(pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^(lambda-1)) *
#                pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
#
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#
#              e <- exp(eta)
#
#              lambda * (pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^(lambda-1))*
#                (exp(eta - e) * (1 - e) - ((1 - lambda) / pmax(pmin(-expm1(-exp(eta)),
#                                                                    1 - .Machine$double.eps), .Machine$double.eps)) *
#                   (pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)^2))
#            }
#
#            rho <- function(eta, lambda){
#              (pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^lambda)*
#                log(pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps))
#            }
#
#            rho2 <- function(eta, lambda){
#              (pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^lambda) *
#                (log(pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#              (pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^(lambda-1))*
#                pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps) *
#                (1 + lambda * log(pmax(pmin(-expm1(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)))
#            }
#
#            plink <- TRUE
#
#            name <- "power complement log-log"
#          },
#
#          # Reversal power logit ----
#          rplogit = {
#
#            linkfun <- function(mu, lambda = NULL, ...) -stats::qlogis((1 - mu)^(1/lambda))
#
#            linkinv <- function(eta, lambda) 1 - stats::plogis(-eta)^lambda
#
#            mu.eta <- function(eta, lambda){
#
#              lambda * (stats::plogis(-eta)^(lambda - 1)) *
#                pmax(stats::dlogis(-eta), .Machine$double.eps)
#
#            }
#
#            mu2.eta2 <- function(eta, lambda){
#
#              e <- exp(eta)
#
#              lambda * (stats::plogis(-eta)^(lambda - 1)) *
#                (exp(-eta) * expm1(-eta) / (1 + exp(-eta))^3 +
#                   ((1 - lambda) / (stats::plogis(-eta))) *
#                   ( pmax(stats::dlogis(-eta), .Machine$double.eps)^2))
#            }
#
#            rho <- function(eta, lambda){
#
#              - (stats::plogis(-eta)^lambda) *
#                log((stats::plogis(-eta)))
#
#            }
#
#            rho2 <- function(eta, lambda){
#
#              - (stats::plogis(-eta)^lambda) *
#                (log((stats::plogis(-eta)))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              (stats::plogis(-eta)^(lambda - 1)) *
#                pmax(stats::dlogis(-eta), .Machine$double.eps) * (1 + lambda * log((stats::plogis(-eta))))
#            }
#
#            plink <- TRUE
#
#            name <- "reversal power logit"
#          },
#
#          # Reversal power probit --------------------------------------------------------------------------
#          rpprobit = {
#
#            linkfun <- function(mu, lambda, ...) -stats::qnorm((1 - mu)^(1/lambda))
#
#            linkinv <- function(eta, lambda) {
#              1 - stats::pnorm(-eta)^lambda
#            }
#
#            mu.eta <- function(eta, lambda){
#
#              lambda * (stats::pnorm(-eta)^(lambda - 1)) *
#                pmax(stats::dnorm(-eta), .Machine$double.eps)
#            }
#
#            mu2.eta2 <- function(eta, lambda){
#
#              h_inv2.eta2 <- eta * pmax(stats::dnorm(-eta), .Machine$double.eps)
#
#              lambda * (stats::pnorm(-eta)^(lambda - 1)) *
#                ((pmax(stats::dnorm(-eta), .Machine$double.eps)^2) *
#                   (1 - lambda) / stats::pnorm(-eta)  - h_inv2.eta2)
#
#            }
#
#            rho <- function(eta, lambda){
#
#              - (stats::pnorm(-eta)^lambda) * log(stats::pnorm(-eta))
#            }
#
#            rho2 <- function(eta, lambda){
#
#              - (stats::pnorm(-eta)^lambda) * (log(stats::pnorm(-eta))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              (stats::pnorm(-eta)^(lambda - 1)) *
#                pmax(stats::dnorm(-eta), .Machine$double.eps) *
#                (1 + lambda * log(stats::pnorm(-eta)))
#            }
#
#            plink <- TRUE
#
#            name <- "reversal power probit"
#          },
#
#          # Reversal power cauchit ------------------------------------------------------------------
#          rpcauchit = {
#
#            linkfun <- function(mu, lambda, ...) -stats::qcauchy((1 - mu)^(1 / lambda))
#
#            linkinv <- function(eta, lambda) {
#              1 - stats::pcauchy(-eta)^lambda
#            }
#
#            mu.eta <- function(eta, lambda){
#              lambda * (stats::pcauchy(-eta)^(lambda - 1)) * (pmax(stats::dcauchy(-eta), .Machine$double.eps))
#            }
#
#            mu2.eta2 <- function(eta, lambda){
#              h_inv2.eta2 <- 2 * eta * pmax(stats::dcauchy(-eta), .Machine$double.eps) / (1 + eta^2)
#
#              lambda * (stats::pcauchy(-eta)^(lambda - 1)) *
#                (- h_inv2.eta2 + ((1 - lambda) / stats::pcauchy(-eta)) *
#                   (pmax(stats::dcauchy(-eta), .Machine$double.eps))^2)
#            }
#
#            rho <- function(eta, lambda){
#              - (stats::pcauchy(-eta)^lambda) * log(stats::pcauchy(-eta))
#            }
#
#            rho2 <- function(eta, lambda){
#              - (stats::pcauchy(-eta)^lambda) * (log(stats::pcauchy(-eta))^2)
#            }
#
#            rho.eta <- function(eta, lambda){
#              (stats::pcauchy(-eta)^(lambda - 1)) *
#                pmax(stats::dcauchy(-eta), .Machine$double.eps)*
#                (1 + lambda * log(stats::pcauchy(-eta)))
#            }
#
#            plink <- TRUE
#
#            name <- "reversal power cauchit"
#          },
#
#          # Reversal power log-log ------------------------------------------------------------------
#          rploglog = {
#
#            linkfun <- function(mu, lambda, ...) log(-log((1 - mu)^(1 / lambda)))
#
#            linkinv <- function(eta, lambda){
#              1 - pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^lambda
#            }
#
#            mu.eta <- function(eta, lambda) {
#              lambda * (pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)^(lambda - 1)) *
#                pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#
#              h_inv <- pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              h_inv2.eta2 <- exp(eta) * exp(-exp(eta)) * expm1(eta)
#
#              lambda * h_inv^(lambda - 1) *
#                ((pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps)^2) * (1 - lambda) / h_inv -
#                   h_inv2.eta2)
#            }
#
#            rho <- function(eta, lambda){
#
#              h_inv <- pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#
#              - h_inv^lambda * log(h_inv)
#
#            }
#
#            rho2 <- function(eta, lambda){
#
#              h_inv <- pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#
#              - h_inv^lambda * log(h_inv)^2
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              h_inv <- pmax(pmin(exp(-exp(eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#
#              h_inv^(lambda - 1) *
#                pmax(exp(eta) * exp(-exp(eta)), .Machine$double.eps) * (1 + lambda * log(h_inv))
#            }
#
#            plink <- TRUE
#
#            name <- "reversal power log-log"
#          },
#
#          # Reversal power complement log-log -------------------------------------------------------
#          rpcloglog = {
#
#            linkfun <- function(mu, lambda, ...) -log(-log(1 - (1 - mu)^(1 / lambda)))
#
#            linkinv <- function(eta, lambda){
#              1 - pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)^lambda
#            }
#
#            mu.eta <- function(eta, lambda) {
#
#              h_inv <- pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#
#              lambda * h_inv^(lambda-1) *
#                pmax(exp(-eta) * exp(-exp(-eta)), .Machine$double.eps)
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#
#              e <- exp(-eta)
#
#              h_inv <- pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              h_inv2.eta2 <- -exp(-eta) * exp(-exp(-eta)) * expm1(-eta)
#
#              lambda * h_inv^(lambda - 1) *
#                (- h_inv2.eta2 +
#                   ((1 - lambda) / h_inv) *
#                   (pmax(exp(-eta) * exp(-exp(-eta)), .Machine$double.eps)^2))
#            }
#
#            rho <- function(eta, lambda){
#              h_inv <- pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              - h_inv^lambda * log(h_inv)
#            }
#
#            rho2 <- function(eta, lambda){
#              h_inv <- pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#              - h_inv^lambda * log(h_inv)^2
#            }
#
#            rho.eta <- function(eta, lambda){
#
#              h_inv <- pmax(pmin(-expm1(-exp(-eta)), 1 - .Machine$double.eps), .Machine$double.eps)
#
#              h_inv^(lambda - 1) *
#                pmax(exp(-eta) * exp(-exp(-eta)), .Machine$double.eps) *
#                (1 + lambda * log(h_inv))
#            }
#
#            plink <- TRUE
#
#            name <- "reversal power complement log-log"
#          },
#
#          # Reversal Aranda-Ordaz -------------------------------------------------------------------
#          raordaz = {
#
#            linkfun <- function(mu, lambda, ...) log(lambda / (mu^(-lambda) - 1))
#
#            linkinv <- function(eta, lambda){
#              e <- pmax(exp(-eta), .Machine$double.ep)
#              (1 + lambda * e)^(-1 / lambda)
#            }
#
#            mu.eta <- function(eta, lambda) {
#              e <- pmax(exp(-eta), .Machine$double.ep)
#              e*(1 + lambda*e)^(-1/lambda - 1)
#            }
#
#            mu2.eta2 <- function(eta, lambda) {
#              e <- pmax(exp(-eta), .Machine$double.ep)
#              e2 <- pmax(exp(eta), .Machine$double.ep)
#              e*(((1+lambda)/(e2 + lambda)) - 1)/
#                ((1 + lambda*e)^(1/lambda + 1))
#            }
#
#            rho <- function(eta, lambda){
#              e <- pmax(exp(-eta), .Machine$double.ep)
#              e2 <- pmax(exp(eta), .Machine$double.ep)
#              (((1 + lambda*e)^(-1/lambda))/lambda)*
#                (log(1 + lambda*e)/lambda - 1/(lambda + e2))
#            }
#
#            rho2 <- function(eta, lambda){
#              e <- pmax(exp(-eta), .Machine$double.ep)
#              e2 <- pmax(exp(-2*eta), .Machine$double.ep)
#              ((lambda*e + 1)^(-1/lambda))*
#                ((log(1 + lambda*e)/(lambda^2) - e/(lambda*(1 + lambda*e)))^2 -
#                   2*log(1 + lambda*e)/(lambda^3) + 2*e/((lambda^2)*(lambda*e + 1)) +
#                   e2/(lambda*((1 + lambda*e)^2)))
#            }
#
#            rho.eta <- function(eta, lambda){
#              e <- pmax(exp(-eta), .Machine$double.ep)
#              e2 <- pmax(exp(eta), .Machine$double.ep)
#
#              e*(log(1 + lambda*e)/lambda - (1+lambda)/(e2 + lambda))/
#                (lambda*((1 + lambda*e)^(1 + 1/lambda)))
#            }
#
#            plink <- TRUE
#
#            name <- "reversal Aranda-Ordaz"
#          },
#
#          stop(gettextf("%s link not recognised", sQuote(link)), domain = NA))
#   environment(linkfun) <- environment(linkinv) <- environment(mu.eta) <- environment(mu2.eta2) <-
#     environment(rho) <- environment(rho2) <- environment(rho.eta) <- asNamespace("stats")
#   structure(list(linkfun = linkfun, linkinv = linkinv, mu.eta = mu.eta, mu2.eta2 = mu2.eta2,
#                  rho = rho, rho2 = rho2, rho.eta = rho.eta, plink = plink, name = name), class = "link-glm")
#
# }
