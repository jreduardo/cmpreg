#' @name cmpreg-methods
#' @title Method for 'cmpreg' objects
#' @param object an object of class \code{cmpreg}.
#' @param x and object of class \code{cmpreg}.
#' @param digits minimal number of _significant_ digits, see
#'   \code{\link[base]{print.default}}.
#' @param ... Currently not used.
#' @return .
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
NULL

#-----------------------------------------------------------------------
#' @rdname cmpreg-methods
# Print method
print.cmpreg <- function(x, digits = max(3L, getOption("digits") - 3L)) {
  cat("\nCOM-Poisson regression models", sep = "")
  cat("\nCall:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Mean coefficients:", "\n", sep = "")
  print.default(format(x$mean_coefficients, digits = digits),
                print.gap = 2, quote = FALSE)
  cat("\n")
  if (x$dformula == ~ 1) {
    cat("Dispersion coefficient: log(nu) = ",
        format(x$disp_coefficients, digits = digits),
        "\n", sep = "")
  } else {
    cat("Dispersion coefficients:", "\n", sep = "")
    print.default(format(x$disp_coefficients, digits = digits),
                  print.gap = 2, quote = FALSE)
  }
  cat("\n")
  cat("Residual degrees of freedom: ", x$df.residual,
      "\n", sep = "")
  cat("Minus twice the log-likelihood: ", -2 * x$loglik,
      "\n", sep = "")
  invisible(x)
}

#-----------------------------------------------------------------------
#' @rdname cmpreg-methods
# Get the log-likelihood
logLik.cmpreg <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  ll <- object$loglik
  attr(ll, "df") <- object$nobs - object$df.residual
  attr(ll, "nobs") <- object$nobs
  class(ll) <- "logLik"
  return(ll)
}

#-----------------------------------------------------------------------
#' @rdname cmpreg-methods
#' @param what a character indicating which parameter coefficient is
#'   required, parameters for the \code{"mean"} or for the \code{"mean"}
#'   model.
# Get the parameter estimates
coef.cmpreg <- function(object,
                        what = c("all", "mean", "dispersion"), ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  type <- match.arg(type)
  out <- switch(type,
                "all"        = object$coefficients,
                "mean"       = object$mean_coefficients,
                "dispersion" = object$disp_coefficients)
  return(out)
}

#-----------------------------------------------------------------------
#' @rdname cmpreg-methods
# Get the variance-covariance matrix
vcov.cmpreg <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  return(object$vcov)
}
