#-----------------------------------------------------------------------
#' @title Summary of the COM-Poisson models (individual t-tests)
#' @param object an object of class \code{cmpreg}, a result of call
#'   \code{\link{cmp}(...)}.
#' @param correlation logical; if \code{TRUE}, the correlation matrix of
#'   the estimated parameters is returned and printed. Default is
#'   \code{FALSE}.
#' @param ... Currently not used.
#' @return an object of class \code{"summary.cmpreg"}, a list with
#'   components.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats cov2cor pnorm printCoefmat
#' @export
#'
summary.cmpreg <- function(object, correlation = FALSE, ...) {
  if (is.null(object$vcov))
    stop(paste("Refit the model with `hessian=TRUE` to compute",
               "the standard errors."))
  #------------------------------------------
  p <- ncol(object$data$X)
  q <- ncol(object$data$Z)
  #------------------------------------------
  estimates  <- object$coefficients
  stderror   <- sqrt(diag(object$vcov))
  zvalue     <- estimates/stderror
  pvalue     <- 2 * pnorm(-abs(zvalue))
  ctableall  <- cbind("Estimate"    = estimates,
                      "Std. Error"  = stderror,
                      "Z value"     = zvalue,
                      "Pr(>|z|)"    = pvalue)
  ctable <- list(mean       = ctableall[1:p, ,drop = FALSE],
                 dispersion = ctableall[1:q + p, ,drop = FALSE])
  ctable <- mapply(`rownames<-`,
                   ctable,
                   lapply(object$data[-3], colnames))
  #------------------------------------------
  if (correlation) {
    ctable$correlation <- cov2cor(object$vcov)
  }
  #--------------------------------------------
  out <- list(coeftable   = ctable,
              loglik      = object$loglik,
              df.residual = object$df.residual,
              correlation = correlation,
              call        = object$call)
  class(out) <- "summary.cmpreg"
  return(out)
}

#-----------------------------------------------------------------------
# Print method for summary COM-Poisson models
#' @rdname cmpreg-methods
#'
print.summary.cmpreg <- function(x,
                                 digits = max(3L, getOption("digits") - 3L),
                                 ...) {
  cat("\nIndividual Wald-tests for COM-Poisson regression models",
      sep = "")
  cat("\nCall:  ",
      paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")
  cat("Mean coefficients:", "\n", sep = "")
  printCoefmat(x$coeftable$mean,
               digits = digits,
               has.Pvalue = TRUE)
  cat("\n")
  cat("Dispersion coefficients:", "\n", sep = "")
  printCoefmat(x$coeftable$dispersion,
               digits = digits,
               has.Pvalue = TRUE)
  cat("\n")
  if (x$correlation) {
    cat("Correlation of coefficients:", "\n", sep = "")
    corr <- x$coeftable$correlation
    corr <- format(round(corr, 2L), nsmall = 2L,
                   digits = digits)
    corr[!lower.tri(corr)] <- ""
    print(corr[-1, -ncol(corr), drop = FALSE], quote = FALSE)
    cat("\n")
  }
  cat("Residual degrees of freedom: ", x$df.residual,
      "\n", sep = "")
  cat("Minus twice the log-likelihood: ", -2 * x$loglik,
      "\n", sep = "")
  invisible(x)
}

#-----------------------------------------------------------------------
#' @title Likelihood ratio tests for nested COM-Poisson models
#' @param object an object of class \code{cmpreg}, a result of call
#'   \code{\link{cmp}(...)}.
#' @param heading logical; if \code{TRUE}, the model formulas are
#'   printed as heading. Default is \code{TRUE}.
#' @param ... Currently not used.
#' @return an object of class \code{"anova"}, a table with compenents
#'   for LR test.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats pchisq
#' @export
#'
anova.cmpreg <- function(object, ..., heading = TRUE) {
  dots <- list(...)
  iscmp <- class(object) == "cmpreg"
  #------------------------------------------
  # Organize the list of models
  if (iscmp & length(dots) == 0)
    stop("Currently, this method only compare two (or more) models.")
  if (!iscmp & length(dots) == 0) mlist <- object
  if (!iscmp & length(dots) > 0)  mlist <- c(object, dots)
  if (iscmp & length(dots)  > 0)  mlist <- c(list(object), dots)
  # print(str(mlist))
  #------------------------------------------
  # Test wheter models are comparable
  #   TODO: check wheter models are nested.
  if (any(!vapply(mlist, inherits, what = "cmpreg", 0L)))
    stop("Not all objects are of class \"cmpreg\"")
  obs <- vapply(mlist, function(x) x$nobs, 0L)
  if (diff(range(obs)) != 0) {
    warning("The models were fitted for different sample sizes.")
  }
  forms_mean <- vapply(mlist, function(x) deparse(x$formula[-2]), "")
  forms_disp <- vapply(mlist, function(x) deparse(x$dformula), "")
  same_mean <- all_identical(forms_mean)
  same_disp <- all_identical(forms_disp)
  if (same_mean & same_disp) {
    warning("The models are equal.")
    heading <- ""
  }
  #--------------------------------------------
  # Compute stats
  rds <- vapply(mlist, function(x) x$df.residual, 0L)
  lls <- vapply(mlist, function(x) x$loglik, 0)
  aic <- -2 * lls + 2 * (obs - rds)
  bic <- -2 * lls + log(obs) * (obs - rds)
  lrs <- c(NA, abs(diff(-2 * lls)))
  dfs <- c(NA, abs(diff(rds)))
  pvs <- pchisq(q = lrs, df = dfs, lower.tail = FALSE)
  tab <- cbind("AIC"        = aic,
               "BIC"        = bic,
               "Resid.df"  = rds,
               "loglik"     = lls,
               "Chisq.df"   = dfs,
               "Chisq"      = lrs,
               "Pr(>Chisq)" = pvs)
  rownames(tab) <- sprintf("Model %i", seq(mlist))
  #--------------------------------------------
  # Heading
  if (heading) {
    forms_mean <- gsub("~", " ~ ", forms_mean)
    forms_disp <- gsub("~", " ~ ", forms_disp)
    blank <- paste(rep(" ", 7 + nchar(max(seq(mlist)))), collapse = "")
    if (same_mean & !same_disp) {
      heading <- paste(sprintf("Model %s: log(nu)%s",
                               seq(mlist), forms_disp),
                       collapse = "\n")
    }
    if (!same_mean & same_disp) {
      heading <- paste(sprintf("Model %s: log(mu)%s",
                               seq(mlist), forms_mean),
                       collapse = "\n")
    }
    if (!same_mean & !same_disp) {
      heading <- paste(sprintf("Model %i: log(mu)%s\n%s log(nu)%s",
                               seq(mlist),
                               forms_mean,
                               rep(blank, length(mlist)),
                               forms_disp),
                       collapse = "\n")
    }
  } else {
    heading <- NULL
  }
  attr(tab, "heading") <- paste0("\nLikelihood ratio test for ",
                                 "COM-Poisson regression models",
                                 "\n\n",  heading, "\n")
  class(tab) <- "anova"
  return(tab)
}

#-----------------------------------------------------------------------
#' @title Likelihood ratio test for equidispersion
#' @param object an object (or a list of objects) of class
#'   \code{"cmpreg"}.
#' @param ... possibly, others fitted models (objects of class
#'   \code{"cmpreg"}).
#' @param heading logical; if \code{TRUE}, the model formulas are
#'   printed as heading. Default is \code{TRUE}.
#' @details Note that, when more than one model is passed the models
#'   \strong{are not compared among themselves}. Here, each model is
#'   compared to Poisson and them organized in a table.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @export
#'
equitest <- function (object, ..., heading = TRUE) {
  UseMethod("equitest", object)
}

#' @export
equitest.cmpreg <- function(object, ..., heading = TRUE) {
  dots <- list(...)
  iscmp <- class(object) == "cmpreg"
  #------------------------------------------
  # Organize the list of models
  if (!iscmp &  is.null(dots)) mlist <- object
  if ( iscmp &  is.null(dots)) mlist <- list(object)
  if (!iscmp & !is.null(dots)) mlist <- c(object, dots)
  if ( iscmp & !is.null(dots)) mlist <- c(list(object), dots)
  if (any(!vapply(mlist, inherits, what = "cmpreg", 0L)))
    stop("Not all objects are of class \"cmpreg\"")
  #--------------------------------------------
  lpo <- vapply(mlist, function(x) x$poissonfit$loglik, 0)
  npo <- vapply(mlist, function(x) length(x$poissonfit$coef), 0L)
  rds <- vapply(mlist, function(x) x$df.residual, 0L)
  lls <- vapply(mlist, function(x) x$loglik, 0)
  obs <- vapply(mlist, function(x) x$nobs, 0L)
  lrs <- 2 * (lls - lpo)
  dfs <- obs - rds - npo
  pvs <- pchisq(q = lrs, df = dfs, lower.tail = FALSE)
  tab <- cbind("Resid.df"   = rds,
               "Loglik"     = lls,
               "LRT_stat"   = lrs,
               "LRT_df"     = dfs,
               "Pr(>LRT_stat)" = pvs)
  rnames <- sprintf("Model %i", seq_along(mlist))
  rownames(tab) <- rnames
  #--------------------------------------------
  if (heading) {
    calls <- vapply(mlist, function(x)
      deparse(x$call, width.cutoff = 500), "")
    calls <- paste(paste0(rnames, ":"), calls)
    calls <- gsub(", ", ",\n             ", calls)
    calls <- paste(calls, collapse = "\n")
  } else {
    calls <- NULL
  }
  attr(tab, "heading") <- calls
  class(tab) <- "equitest.cmpreg"
  return(tab)
}

#' @export
equitest.list <- function(object, ..., heading = TRUE) {
  out <- equitest.cmpreg(object, ...)
  if (heading & !is.null(names(object))) {
    calls <- sprintf("Model %i: %s", seq_along(object), names(object))
    attr(out, "heading") <- paste(calls, collapse = "\n")
  }
  return(out)
}

print.equitest.cmpreg <- function(x,
                                  digits = max(getOption("digits") - 2L, 3L),
                                  signif.stars = getOption("show.signif.stars"),
                                  ...) {
  cat("\nLikelihood ratio test for equidispersion", "\n\n")
  if (!is.null(attr(x, "heading"))) {
    cat(attr(x, "heading"), "\n\n")
  }
  printCoefmat(x, cs.ind = NA, tst.ind = 2:3, zap.ind = 5L,
               digits = digits,
               signif.stars = signif.stars,
               ...)
  invisible(x)
}
