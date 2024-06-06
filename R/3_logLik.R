## Log-likelihood ----
ll_ug <- function(par, y, X, Z, link = "aordaz", sigma.link = "aordaz")
{

  k <- ncol(X); l <- ncol(Z); n <- length(y)

  beta <- as.vector(par[1:k])
  gamma <- as.vector(par[1:l + k])
  lambda1 <- NULL
  lambda2 <- NULL

  nlambdas <- make.plink(link)$plink + make.plink(sigma.link)$plink
  if (make.plink(link)$plink)
    lambda1 <-  as.numeric(par[k + l + 1])
  if (make.plink(sigma.link)$plink)
    lambda2 <- as.numeric(par[k + l + nlambdas])

  mu <- make.plink(link)$linkinv(X%*%beta, lambda1)
  sigma <- make.plink(sigma.link)$linkinv(Z%*%gamma, lambda2)

  if (any(mu <= 0) | any(mu >= 1) | any(sigma <= 0) | any(sigma >= 1) |
      any(lambda1 < 0) | any(lambda2 < 0)){
    NaN
  }else{
    ll <- suppressWarnings(dugamma(y, mu, sigma, log.p = TRUE))

    if (any(!is.finite(ll)))
      NaN
    else
      sum(ll)
  }
}


## Score function ----
U_ug <- function(par, y, X, Z, link = "aordaz", sigma.link = "aordaz")
{
  k <- ncol(X); l <- ncol(Z); n <- length(y)

  beta <- as.vector(par[1:k])
  gamma <- as.vector(par[1:l + k])
  lambda1 <- NULL
  lambda2 <- NULL

  g <- make.plink

  nlambdas <- g(link)$plink + g(sigma.link)$plink
  if (g(link)$plink)
    lambda1 <-  as.numeric(par[k + l + 1])
  if (g(sigma.link)$plink)
    lambda2 <- as.numeric(par[k + l + nlambdas])

  mu <- g(link)$linkinv(X%*%beta, lambda1)
  sigma <- g(sigma.link)$linkinv(Z%*%gamma, lambda2)

  phi <- as.vector(((1-(sigma^2))/(sigma^2)))
  mu_star <- mu^(1/phi)
  d <- (mu_star)/(1-mu_star)

  T1 <- diag(as.vector(g(link)$mu.eta(X%*%beta, lambda1)))
  T2 <- diag(as.vector(g(sigma.link)$mu.eta(Z%*%gamma, lambda2)))
  a1 <- (1 - mu_star - (sigma^2 * mu_star * log(y)) /
           (sigma^2 - 1))/(mu * (1-mu_star)^2)
  a2 <- (2 * log(mu) *
           (mu_star * (sigma^2 + sigma^2 * log(y) - 1) + 1 - sigma^2)) /
    ((1-sigma^2)^2 * sigma * (1-mu_star)^2) +
    2 * (digamma(phi) - log(-log(y)) - log(d)) / sigma^3

  Ub <- t(X) %*% T1 %*% a1
  Ug <- t(Z) %*% T2 %*% a2

  if (g(link)$plink == TRUE){
    rho1 <- g(link)$rho(X%*%beta, lambda1)
    Ul1 <- sum(rho1 * a1)
  }else{
    Ul1 <- NULL
  }

  if (g(sigma.link)$plink == TRUE){
    rho2 <- g(sigma.link)$rho(Z%*%gamma, lambda2)
    Ul2 <- sum(rho2 * a2)
  }else{
    Ul2 <- NULL
  }

  c(Ub, Ug, Ul1, Ul2)
}

## Fisher information matrix ----
K_ug <- function(par, X, Z, link = "aordaz", sigma.link = "aordaz", inverse = FALSE)
{
  k <- ncol(X); l <- ncol(Z)

  beta <- as.vector(par[1:k])
  gamma <- as.vector(par[1:l + k])
  lambda1 <- NULL
  lambda2 <- NULL

  g <- make.plink

  nlambdas <- g(link)$plink + g(sigma.link)$plink
  if (g(link)$plink)
    lambda1 <-  as.numeric(par[k + l + 1])
  if (g(sigma.link)$plink)
    lambda2 <- as.numeric(par[k + l + nlambdas])

  mu <- g(link)$linkinv(X%*%beta, lambda1)
  sigma <- g(sigma.link)$linkinv(Z%*%gamma, lambda2)

  phi <- as.vector(((1-(sigma^2))/(sigma^2)))
  mu_star <- mu^(1/phi)

  T1 <- diag(as.vector(g(link)$mu.eta(X%*%beta, lambda1)))
  T2 <- diag(as.vector(g(sigma.link)$mu.eta(Z%*%gamma, lambda2)))
  a <- 1/(mu*((1-mu_star)^2))
  b <- 2*mu*log(mu)/(sigma*(1-sigma^2))

  W <- diag(c( (a^2*sigma^2*(1-mu_star)^2/(1-sigma^2)) *((g(link)$mu.eta(X%*%beta, (lambda1)) )^2) ))
  C <- diag(c( (a/(1-sigma^2))*( 2*(1-mu_star)/sigma + b*(sigma^2)/mu ) ))
  V <- diag(c( a^2*sigma^2*(1-mu_star)^2/(1-sigma^2) ))
  Q <- diag(c( (4/(sigma^6))*(trigamma(phi)) +
                 a^2*b^2*sigma^2*((1-mu_star)^2)/(1-sigma^2)+
                 4*a*b*(1-mu_star)/(sigma*(1-sigma^2)) ))

  Kbb <- t(X)%*%W%*%X
  Kbg <- t(X)%*%C%*%T1%*%T2%*%Z
  Kgb <- t(Kbg)
  Kgg <- t(Z)%*%Q%*%T2%*%T2%*%Z

  if (g(link)$plink){
    eta1 <- X%*%beta
    rho1 <- g(link)$rho(X%*%beta, lambda1)
    Kbl1 <- t(X)%*%V%*%T1%*%rho1
    Kl1b <- t(Kbl1)
    Kgl1 <- t(Z)%*%C%*%T2%*%rho1
    Kl1g <- t(Kgl1)
    Kl1l1 <- t(rho1)%*%V%*%rho1
    ind_l1 <- 1
  }else{
    Kbl1 <- Kl1b <- Kgl1 <- Kl1g <-Kl1l1 <- NULL
    ind_l1 <- 0
  }

  if (g(sigma.link)$plink){
    eta2 <- Z%*%gamma
    rho2 <- g(sigma.link)$rho(Z%*%gamma, lambda2)
    Kbl2 <- t(X)%*%C%*%T1%*%rho2
    Kl2b <- t(Kbl2)
    Kgl2 <- t(Z)%*%Q%*%T2%*%rho2
    Kl2g <- t(Kgl2)
    Kl2l2 <- t(rho2)%*%Q%*%rho2
    ind_l2 <- 1
  }else{
    Kbl2 <- Kl2b <- Kgl2 <- Kl2g <-Kl2l2 <- NULL
    ind_l2 <- 0
  }

  if (ind_l1 == 1 && ind_l2 == 1){
    Kl1l2 <- t(rho1)%*%C%*%rho2
    Kl2l1 <- t(Kl1l2)
  }else{
    Kl1l2 <- Kl2l1 <- NULL
  }

  K <- rbind(cbind(Kbb,Kbg,Kbl1,Kbl2),
             cbind(Kgb,Kgg,Kgl1,Kgl2),
             cbind(Kl1b,Kl1g,Kl1l1,Kl1l2),
             cbind(Kl2b,Kl2g,Kl2l1,Kl2l2))

  if(!inverse) K else chol2inv(chol(K))

}



## Hessian matrix ----
J_ug <- function(par, y, X, Z, link = "aordaz", sigma.link = "aordaz", inverse = FALSE)
{
  k <- ncol(X); l <- ncol(Z); n <- length(y)

  beta <- as.vector(par[1:k])
  gamma <- as.vector(par[1:l + k])
  lambda1 <- NULL
  lambda2 <- NULL

  g <- make.plink

  nlambdas <- g(link)$plink + g(sigma.link)$plink
  if (g(link)$plink)
    lambda1 <-  as.numeric(par[k + l + 1])
  if (g(sigma.link)$plink)
    lambda2 <- as.numeric(par[k + l + nlambdas])

  mu <- g(link)$linkinv(X%*%beta, lambda1)
  sigma <- g(sigma.link)$linkinv(Z%*%gamma, lambda2)

  # Necessary quantities
  phi <- as.vector(((1-(sigma^2))/(sigma^2)))
  mu_star <- mu^(1/phi)
  y_star <- 1+mu_star*log(y)*(sigma^2)/(1-sigma^2)
  d <- (mu_star)/(1-mu_star)
  mu_dag <- digamma(phi) - log(d)
  y_dag <- log(-log(y))

  T1 <- diag(as.vector(g(link)$mu.eta(X%*%beta, lambda1)))
  T1_prime <- diag(as.vector(g(link)$mu2.eta2(X%*%beta, lambda1)))
  T2 <- diag(as.vector(g(sigma.link)$mu.eta(Z%*%gamma, lambda2)))
  T2_prime <- diag(as.vector(g(sigma.link)$mu2.eta2(Z%*%gamma, lambda2)))

  a <- 1/(mu*((1-mu_star)^2))
  b <- 2*mu*log(mu)/(sigma*(1-sigma^2))

  V_star <- diag(c( (a^2 * (1 - mu_star) / (1 - sigma^2)) *
                      (mu_star * (1 - mu_star) -
                         (1 - y_star)*(2*sigma^2 + mu_star - 1)) -
                      a^2 * (1 - mu_star)^2   ))

  Q_star <- diag(c(   6 * (y_dag - mu_dag) / (sigma^4) -
                        4 * trigamma(phi) / (sigma^6) -
                        (a^2) * (b^2) * (sigma^2) * (1 - y_star - mu_star + mu_star^2 -
                                                       mu_star^3 + y_star * (mu_star^2)) / (1 - sigma^2) -
                        a * b * (4 + 3 * mu_star * (sigma^2) - 3 * y_star * (sigma^2) -
                                   3 * mu_star - y_star)/(sigma * (1 - sigma^2)) ))

  C_star <- diag(c(  -a * (2*(1 - y_star)/sigma + b * (sigma^2) *
                             (1 - y_star - y_star * mu_star + mu_star^2)/
                             (mu * (1-mu_star)))/(1 - sigma^2)   ))

  Dlmu <- diag(c(a * (y_star - mu_star)))
  Dlsigma <- diag(c(a*b*(y_star - mu_star) - 2*(y_dag - mu_dag)/(sigma^3)))

  Jbb <- t(X)%*%(V_star%*%(T1^2) + Dlmu%*%T1_prime)%*%X
  Jbg <- t(X)%*%C_star%*%T1%*%T2%*%Z
  Jgb <- t(Jbg)
  Jgg <- t(Z)%*%(Q_star%*%(T2^2) + Dlsigma%*%T2_prime)%*%Z

  if (g(link)$plink){
    rho1 <- g(link)$rho(X%*%beta, lambda1)
    rho1_prime <- g(link)$rho2(X%*%beta, lambda1)
    rho1_eta <- g(link)$rho.eta(X%*%beta, lambda1)

    Jbl1 <- t(X)%*%(V_star%*%T1%*%rho1 + Dlmu%*%rho1_eta)
    Jl1b <- t(Jbl1)
    Jgl1 <- t(Z)%*%C_star%*%T2%*%rho1
    Jl1g <- t(Jgl1)
    Jl1l1 <- t(rho1)%*%V_star%*%rho1 + matrix(1, ncol = n)%*%Dlmu%*%rho1_prime
    ind_l1 <- 1
  }else{
    Jbl1 <- Jl1b <- Jgl1 <- Jl1g <-Jl1l1 <- NULL
    ind_l1 <- 0
  }

  if (g(sigma.link)$plink){
    rho2 <- g(sigma.link)$rho(Z%*%gamma, lambda2)
    rho2_prime <- g(sigma.link)$rho2(Z%*%gamma, lambda2)
    rho2_eta <- g(sigma.link)$rho.eta(Z%*%gamma, lambda2)

    Jbl2 <- t(X)%*%C_star%*%T1%*%rho2
    Jl2b <- t(Jbl2)
    Jgl2 <- t(Z)%*%(Q_star%*%T2%*%rho2 + Dlsigma%*%rho2_eta)
    Jl2g <- t(Jgl2)
    Jl2l2 <- t(rho2)%*%Q_star%*%rho2 + matrix(1, ncol = n)%*%Dlsigma%*%rho2_prime
    ind_l2 <- 1
  }else{

    Jbl2 <- Jl2b <- Jgl2 <- Jl2g <- Jl2l2 <- NULL
    ind_l2 <- 0
  }

  if (ind_l1 == 1 && ind_l2 == 1){
    Jl1l2 <- t(rho1)%*%C_star%*%rho2
    Jl2l1 <- t(Jl1l2)
  }else{
    Jl1l2 <- Jl2l1 <- NULL
  }

  J <- rbind(cbind(Jbb, Jbg, Jbl1, Jbl2),
             cbind(Jgb, Jgg, Jgl1, Jgl2),
             cbind(Jl1b,Jl1g, Jl1l1, Jl1l2),
             cbind(Jl2b,Jl2g, Jl2l1, Jl2l2))

  if(!inverse) J else chol2inv(chol(J))

}
