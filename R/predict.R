#-----------------------------------------------------------------------
#' @title Predict method for COM-Poisson models
#' @param object a fitted object of class from \code{"glm"}.
#' @param newdata optionally, a data frame in which to look for
#'   variables with which to predict. If omitted, the fitted linear
#'   predictors are used.
#' @param newmatrices optionally, a list with named design matrices
#'   (\code{"X"} for mean model and \code{"Z"} for dispersion model)
#'   used to predict. If omitted, the fitted linear predictors are used.
#' @param what a character indicating which parameter coefficient is
#'   required, parameters for the \code{"mean"} or for the
#'   \code{"dispersion"} model.
#' @param type the type of prediction required. The default
#'   \code{"link"} is on the scale of the linear predictors; the
#'   alternative \code{"response"} is on the scale of the response
#'   variable.
#' @param se.fit logical switch indicating if standard errors are
#'   required. Default is \code{FALSE}.
#' @param augment_data logical indicating if \code{newdata} should be
#'   augmented with the predict values and, possibly, standard
#'   errors. Default is \code{FALSE}.
#' @param ... currently not used.
#' @return a tibble with fitted values (code{fit}) and, possibly,
#'   standard erros (\code{ste}). If \code{augment_data}, the result
#'   will be the \code{newdata} with new columns \code{fit}, fitted
#'   values, and code{ste}, standard errors.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
predict.cmpreg <- function(object,
                           newdata,
                           newmatrices = NULL,
                           what = c("mean", "dispersion", "all"),
                           type = c("link", "response"),
                           se.fit = FALSE,
                           augment_data = FALSE,
                           ...) {
  #------------------------------------------
  type <- match.arg(type)
  what <- match.arg(what)
  missingdata <- missing(newdata)
  missingmats <- is.null(newmatrices)
  #------------------------------------------
  Vcov <- object$vcov
  indb <- grep("beta", rownames(Vcov))
  indg <- grep("gama", rownames(Vcov))
  Vbeta <- Vcov[indb, indb, drop = FALSE]
  Vgama <- Vcov[indg, indg, drop = FALSE]
  Vbega <- Vcov[indb, indg, drop = FALSE]
  #------------------------------------------
  if (!missingdata & !missingmats) {
    stop("Use only 'newdata' or 'newmatrices'.")
  }
  if (missingdata & missingmats) {
    mats <- switch(what,
                   "mean"         = list(X = object$data$X),
                   "dispersion"   = list(Z = object$data$Z),
                   "all"          = list(X = object$data$X,
                                         Z = object$data$Z))
  }
  if (!missingmats & missingdata) {
    mats <- switch(what,
                   "mean"         = list(X = newmatrices$X),
                   "dispersion"   = list(Z = newmatrices$Z),
                   "all"          = list(X = newmatrices$X,
                                         Z = newmatrices$Z))
  }
  if (missingmats & !missingdata) {
    mats <- switch(
      what,
      "mean"         =
        list(X = model.matrix(object$formula[-2], newdata)),
      "dispersion"   =
        list(Z = model.matrix(object$dformula,    newdata)),
      "all"          =
        list(X = model.matrix(object$formula[-2], newdata),
             Z = model.matrix(object$dformula,    newdata)))
  }
  #------------------------------------------
  if (what == "mean") {
    X <- as.matrix(mats$X)
    Vcond <- Vbeta - tcrossprod(Vbega %*% solve(Vgama), Vbega)
    beta <- object$mean_coefficients
    out <- data.frame(fit = c(X %*% beta))
    if (se.fit) {
      ste <- sqrt(diag(tcrossprod(X %*% Vcond, X)))
      out <- cbind(out, ste = ste)
    }
  }
  #------------------------------------------
  if (what == "dispersion") {
    Z <- as.matrix(mats$Z)
    Vcond <- Vgama - crossprod(Vbega, solve(Vbeta) %*% Vbega)
    gama <- object$disp_coefficients
    out <- data.frame(fit = c(Z %*% gama))
    if (se.fit) {
      ste <- sqrt(diag(tcrossprod(Z %*% Vcond, Z)))
      out <- cbind(out, ste = ste)
    }
  }
  #--------------------------------------------
  if (what == "all") {
    Z <- as.matrix(mats$Z)
    X <- as.matrix(mats$X)
    Vcondbeta <- Vbeta - tcrossprod(Vbega %*% solve(Vgama), Vbega)
    Vcondgama <- Vgama - crossprod(Vbega, solve(Vbeta) %*% Vbega)
    beta <- object$mean_coefficients
    gama <- object$disp_coefficients
    out <- data.frame(what = rep(c("mean", "dispersion"),
                                 c(nrow(X), nrow(Z))),
                      fit = c(X %*% beta, Z %*% gama))
    if (se.fit) {
      ste_mean <- sqrt(diag(tcrossprod(X %*% Vcondbeta, X)))
      ste_disp <- sqrt(diag(tcrossprod(Z %*% Vcondgama, Z)))
      out <- cbind(out, ste = c(ste_mean, ste_disp))
    }
  }
  if (type == "response") {
    out$fit <- exp(out$fit)
    if (se.fit) {
      out$ste <- out$fit * out$ste
    }
  }
  if (augment_data & missingdata)
    warning("Needs 'newdata' for augment data.")
  if (augment_data & !missingdata) {
    if(what == "all") newdata <- rbind(newdata, newdata)
    out <- cbind(newdata, out)
    rownames(out) <- NULL
  }
  class(out) <- c("tbl_df", "tbl", "data.frame")
  return(out)
}


#-----------------------------------------------------------------------
# Fitted method
#' @rdname cmpreg-methods
#' @export
fitted.cmpreg <- function(object,
                          what = c("mean", "dispersion", "all"),
                          ...) {
  what <- match.arg(what)
  switch(what,
         "all"        = list(mean       = object$fitted_mean,
                             dispersion = object$fitted_disp),
         "mean"       = object$fitted_mean,
         "dispersion" = object$fitted_disp)
}
