#' @name ugamma
#'
#' @title The Unit Gamma Distribution
#'
#' @description Density, distribution function, quantile function, and random generation for the
#'     unit gamma distribution parameterized in terms of mean (\code{mu}) and dispersion \code{sigma}
#'     parameters.
#'
#' @param x,q vector of quantiles.
#' @param p vector of probabilities.
#' @param n number of random values to return.
#' @param mu vector of means, taking values on (0, 1).
#' @param sigma vector of dispersion parameters, taking values on (0, 1).
#' @param lower.tail logical; if \code{TRUE} (default),
#'     probabilities are \code{P(X <= x)}, otherwise, \code{P(X > x)}.
#' @param log.p logical; if TRUE, probabilities \code{p} are given as
#'  \code{log(p)}.
#'
#' @return \code{dugamma} returns the probability function, \code{pugamma}
#' gives the distribution function, \code{qugamma} gives the quantile function,
#' and \code{rugamma} generates random observations.
#'
#' @details A continuous random variable \eqn{Y} is said to follow a unit gamma distribution with mean
#' \eqn{\mu} and dispersion parameter \eqn{\sigma} if its probability density function is given by
#' \deqn{
#' f(y; \mu,\sigma) = \dfrac{y^{d(\mu,\sigma) - 1} d(\mu,\sigma)^{1/\sigma^2-1}}{\Gamma\left(\dfrac{1}{\sigma^2}-1\right)}\left\{\log \left(\frac{1}{y}\right) \right\}^{1/\sigma^2-2}, \quad y \in (0, 1),
#' }
#' where
#' \deqn{
#'  d(\mu,\sigma)=\frac{\mu^{\sigma^2/(1-\sigma^2)}}{1-\mu^{\sigma^2/(1-\sigma^2)}},
#' }
#' for \eqn{\mu, \sigma \in (0, 1)}.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @examples
#' ## Probability density function for some combinations of
#' ## the parameter values
#' curve(dugamma(x, 0.25, 0.5), col = 1, ylim = c(0, 4), ylab = "Density")
#' curve(dugamma(x, 0.3, 0.5), col = 2, add = TRUE)
#' curve(dugamma(x, 0.5, 0.5), col = 3, add = TRUE)
#' curve(dugamma(x, 0.6, 0.5), col = 4, add = TRUE)
#' curve(dugamma(x, 0.73, 0.5), col = 6, add = TRUE)
#' legend("topleft", c(expression(mu == 0.25~","~ sigma==0.5),
#'                     expression(mu == 0.30~","~ sigma==0.5),
#'                     expression(mu == 0.50~","~ sigma==0.5)),
#'        lty = 1, col = 1:3, bty = "n")
#'
#' legend("top", c(expression(mu == 0.60~","~ sigma==0.5),
#'                 expression(mu == 0.73~","~ sigma==0.5)),
#'        lty = 1, col = c(4, 6), bty = "n")
#'
#'
#' ## Random generation
#' y <- rugamma(1000, 0.25, 0.5)
#'
#' hist(y, prob = TRUE, col = "white")
#' curve(dugamma(x, 0.25, 0.5), col = "blue", add = TRUE, lwd = 2)
#'
#' plot(ecdf(y), col = "grey")
#' curve(pugamma(x, 0.25, 0.5), col = "blue", add = TRUE)
#'
#' plot(ppoints(1000), quantile(y, probs = ppoints(1000)),
#'      xlab = "p", ylab = expression(p-"Quantile"), pch = 16, col = "grey")
#' curve(qugamma(x, 0.25, 0.5), col = "blue", add = TRUE)
NULL

#' @rdname ugamma
#' @export
dugamma <- function(x, mu, sigma, log.p = FALSE){

  if (any(mu <= 0 | mu >= 1))
    warning("The mean parameter must be on (0, 1)")
  if (any(sigma <= 0 | sigma >= 1))
    warning("The dispersion parameter must be on (0, 1)")

  if (is.vector(x))
    x <- matrix(x, ncol = length(x))

  n <- dim(x)[1]
  d <- dim(x)[2]

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)

  if(n%%dim(mu)[1] != 0)
    stop("x and mean have non-conforming size")
  if(n%%dim(sigma)[1] != 0)
    stop("x and sigma have non-conforming size")

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))

  pmf <- matrix(-Inf, n, d)

  # NaN indexation
  pmf[which(mu <= 0 | mu >= 1 |
              sigma <= 0 | sigma >= 1, arr.ind = TRUE)] <- NaN

  # Positive density index
  id <- which(0 < x & x < 1 & !is.nan(pmf), arr.ind = TRUE)

  mu_star <- mu[id]^(sigma[id]^2 / (1 - sigma[id]^2))
  dp <- mu_star / (1 - mu_star)

  pmf[id] <- (1 / sigma[id]^2 - 1) * log(dp) + (dp - 1) * log(x[id]) +
    (1 / sigma[id]^2 - 2) * log(- log(x[id])) - lgamma(1 / sigma[id]^2 - 1)

  if (n == 1) pmf <- as.numeric(pmf)
  if (!log.p) pmf <- exp(pmf)

  pmf
}

#' @rdname ugamma
#' @export
pugamma <- function(q, mu, sigma, lower.tail = TRUE){

  if (any((sigma <= 0) | (sigma >= 1)))
    warning("The dispersion parameter must be on the unit interval (0, 1).")
  if (any((mu <= 0) | (mu >= 1)))
    warning("The mean parameter must be on the unit interval (0, 1).")

  if (is.vector(q))
    q <- matrix(q, ncol = length(q))

  n <- dim(q)[1]
  d <- dim(q)[2]

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)

  if(n%%dim(mu)[1] != 0)
    stop("q and mean have non-conforming size")
  if(n%%dim(sigma)[1] != 0)
    stop("q and sigma have non-conforming size")

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))

  cdf <- matrix(0, n, d)

  # NaN indexation
  cdf[which(mu <= 0 | mu >= 1 |
              sigma <= 0 | sigma >= 1, arr.ind = TRUE)] <- NaN

  # Positive density index
  id <- which(0 < q & !is.nan(cdf), arr.ind = TRUE)

  mu_star <- mu[id]^(sigma[id]^2 / (1 - sigma[id]^2))
  dp <- mu_star / (1 - mu_star)

  cdf[id] <- stats::pgamma(-log(q[id]), 1 / sigma[id]^2 - 1, dp, lower.tail = FALSE)

  if (!lower.tail) cdf <- 1 - cdf
  if (n == 1L) cdf <- as.numeric(cdf)

  cdf
}

#' @rdname ugamma
#' @export
qugamma <- function(p, mu, sigma, lower.tail = TRUE){

  if (any((sigma <= 0) | (sigma >= 1)))
    warning("The dispersion parameter must be on the unit interval (0, 1).")
  if (any((mu <= 0) | (mu >= 1)))
    warning("The mean parameter must be on the unit interval (0, 1).")

  if (is.vector(p))
    p <- matrix(p, ncol = length(p))

  n <- dim(p)[1]
  d <- dim(p)[2]

  mu <- matrix(mu, ncol = d)
  sigma <- matrix(sigma, ncol = d)

  if(n%%dim(mu)[1] != 0)
    stop("x and mean have non-conforming size")
  if(n%%dim(sigma)[1] != 0)
    stop("x and sigma have non-conforming size")

  mu <- do.call(rbind, replicate(n/dim(mu)[1], mu, simplify = FALSE))
  sigma <- do.call(rbind, replicate(n/dim(sigma)[1], sigma, simplify = FALSE))

  qtl <- matrix(0, n, d)

  # NaN indexation
  qtl[which(mu <= 0 | mu >= 1 | sigma <= 0 | sigma >= 1 | p < 0 | p > 1, arr.ind = TRUE)] <- NaN

  id <- which(0 < p & p <= 1 & !is.nan(qtl), arr.ind = TRUE)

  if(!lower.tail) p <- 1 - p

  mu_star <- mu[id]^(sigma[id]^2 / (1 - sigma[id]^2))
  dp <- mu_star / (1 - mu_star)

  qtl[id] <- exp(-stats::qgamma(1 - p[id], 1 / sigma[id]^2 - 1, dp))

  if (n == 1) qtl <- as.numeric(qtl)

  qtl
}

#' @rdname ugamma
#' @export
rugamma <- function(n, mu, sigma){
  #qugamma(stats::runif(n), mu, sigma)

  alpha <- 1/sigma^2 - 1
  beta <- mu^(1 / alpha) / (1 - mu^(1 / alpha))
  exp(-stats::rgamma(n, alpha, beta))

}
