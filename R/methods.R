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
# Print method
#' @rdname cmpreg-methods
#'
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
# Get the log-likelihood
#' @rdname cmpreg-methods
#'
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
# Get the parameter estimates
#' @rdname cmpreg-methods
#' @param what a character indicating which parameter coefficients are
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"dispersion"} model. If \code{"all"}, a list with
#'   coefficients for the \code{mean} and for the \code{dispersion}
#'   model is returned.
#' @export
#'
coef.cmpreg <- function(object,
                        what = c("mean", "dispersion", "all"), ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  what <- match.arg(what)
  out <- switch(what,
                "all"        = list(
                  mean       = object$mean_coefficients,
                  dispersion = object$disp_coefficients),
                "mean"       = object$mean_coefficients,
                "dispersion" = object$disp_coefficients)
  return(out)
}

#-----------------------------------------------------------------------
# Get the variance-covariance matrix
#' @rdname cmpreg-methods
#'
vcov.cmpreg <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  return(object$vcov)
}

#-----------------------------------------------------------------------
# Get the design matrices
#' @rdname cmpreg-methods
#'
model.matrix.cmpreg <- function(object, ...) {
  if (!missing(...))
    warning("extra arguments discarded")
  list(X = object$data$X, Z = object$data$Z)
}
