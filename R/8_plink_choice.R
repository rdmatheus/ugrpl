#' Choosing the Link Function in the Unit Gamma Regression
#'
#' @param object a "\code{ugrpl}" object.
#' @param link,sigma.link a vector of characters specifying which link functions should be
#'     considered for the mean and dispersion parameter submodels, respectively. By default,
#'     possible combinations between the link functions \code{"logit"}, \code{"probit"},
#'     \code{"plogit"}, \code{"pprobit"}, \code{"rplogit"}, \code{"rpprobit"}, \code{"aordaz"},
#'     and \code{"raordaz"} are considered. See \code{\link{make.plink}}.
#' @param criterion criteria for choosing the "best" combination. AIC, BIC, or Upsilon statistics
#'     can be used.
#' @param control a list of control parameters passed as arguments for the \code{optim} function
#'     specified via \code{\link{ug_control}}.
#' @param ... arguments passed to \code{\link{ug_control}}.
#'
#' @return The \code{plink.choice} function returns a list where each component is a unit gamma
#'     regression fit with a combination of the links provided by the \code{link} and \code{sigma.link}
#'     arguments.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
#'
#' @export
#'
plink_choice <- function(object, link, sigma.link,
                         criterion = c("aic", "bic", "upsilon"),
                         control = ug_control(...), ...){

  criterion <- match.arg(criterion, c("aic", "bic", "upsilon"))

  n <- object$nobs

  object_link <- object$link
  object_sigma.link <- object$sigma.link

  if (missing(link)) link <- c("logit", "probit", "plogit", "pprobit",
                               "rplogit", "rpprobit", "aordaz", "raordaz")

  if (missing(sigma.link)) sigma.link <- c("logit", "probit", "plogit", "pprobit",
                                           "rplogit", "rpprobit", "aordaz", "raordaz")

  FIT <- list()
  l <- 1L
  AIC <- BIC <- Upsilon <- matrix(NA, length(link), length(sigma.link))
  pb <- utils::txtProgressBar(max = length(link), style = 3)
  for (i in link){

    m <- 1L

    for (j in sigma.link) {

      if (i != object_link | j != object_sigma.link) {
        FIT[[i]][[j]] <- suppressWarnings(try(stats::update(object, link = i, sigma.link = j, control = control), silent = TRUE))
      } else {
        FIT[[i]][[j]] <- object
      }


      if (!any(grepl("Error", FIT[[i]][[j]]))){

        res <- stats::residuals(FIT[[i]][[j]], type = "quantile")

        AIC[l, m] <- stats::AIC(FIT[[i]][[j]])
        BIC[l, m] <- stats::AIC(FIT[[i]][[j]], k = log(n))
        Upsilon[l, m] <- mean(abs(sort(res) - EnvStats::evNormOrdStats(n = n)))
      }

      m <- m + 1


    }

    l <- l + 1
    utils::setTxtProgressBar(pb, l)
  }

  cat("\n")

  rownames(AIC) <- rownames(BIC) <- rownames(Upsilon) <- link
  colnames(AIC) <- colnames(BIC) <- colnames(Upsilon) <- sigma.link

  FIT$AIC <- AIC
  FIT$BIC <- BIC
  FIT$Upsilon <- Upsilon
  FIT$criterion <- criterion
  FIT$link <- link
  FIT$sigma.link <- sigma.link

  class(FIT) <- "plink_choice"
  FIT

}


#' @name plink_choice-methods
#' @title Methods for \code{"plink_choice"} objects
#' @param x an object of class \code{"plink_choice"}.
#' @param ... further arguments passed to or from other methods.
#'
#' @author Rodrigo M. R. de Medeiros <\email{rodrigo.matheus@live.com}>
NULL

# Print
#' @rdname plink_choice-methods
#' @export
print.plink_choice <- function(x, ...) {

  cat("\nLink function selection for the mean and dispersion submodels")

  AIC <- x$AIC
  BIC <- x$BIC
  Upsilon <- x$Upsilon
  criterion <- x$criterion

  if (criterion == "aic") {
    id <- which(AIC == min(AIC, na.rm = TRUE), arr.ind = TRUE)
    cat("\n\nAIC:\n")
    print(AIC)
  } else if (criterion == "bic") {
    id <- which(BIC == min(BIC, na.rm = TRUE), arr.ind = TRUE)
    cat("\n\nBIC:\n")
    print(BIC)
  } else if (criterion == "upsilon") {
    id <- which(Upsilon == min(Upsilon, na.rm = TRUE), arr.ind = TRUE)
    cat("\n\nUpsilon statistic:\n")
    print(Upsilon)
  }

  cat("\nBest link for the mean submodel:", x$link[id[1]])
  cat("\nBest link for the dispersion submodel:", x$sigma.link[id[2]], "\n")

  invisible(x)
}
