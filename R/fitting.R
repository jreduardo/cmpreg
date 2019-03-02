#' @useDynLib cmpreg, .registration = TRUE
#' @importFrom Rcpp sourceCpp
NULL

#-----------------------------------------------------------------------
#' @rdname llcmp-compute
#' @title Negative of the log-likelihood function from the
#'   reparametrized COM-Poisson model
#' @description Compute the negative of the logarithm of the likelihood
#'   function for a set of observations \code{y} from the reparametrized
#'   COM-Poisson model given the \eqn{\mu} and \eqn{\nu} parameters (see
#'   Details).
#' @details The log-likelihood function is given by
#'   \deqn{\ell(\beta,\nu) = \sum y \log\lambda - \nu\log y! -
#'   \log[Z(\lambda,\nu)],} where \deqn{\log\lambda = \nu \log\left(\mu
#'   - \frac{\nu-1}{2\nu}\right),} and \eqn{\log[Z(\lambda, \nu)]} is a
#'   normalizing constant computed in log space to avoid numerical
#'   issues (see \code{\link{compute_logz}}).
#' @param beta A vector of \eqn{\beta} parameter.
#' @param gama A vector of \eqn{\gamma} parameter.
#' @param X Design matrix related to the (approximate) mean parameter
#'   \eqn{\mu = \exp(X \beta)}.
#' @param Z Design matrix related to the dispersion parameter \eqn{\nu =
#'   \exp(Z \gamma)}.
#' @param y Vector of observed count data.
#' @param params A vector of the model parameters \code{params} =
#'   \code{c(beta, gama)}.
#' @return The computed log-likelihood function.
#' @author Eduardo Jr <edujrrib@gmail.com>
#'
llcmp_fixed <- function(beta, gama, X, Z, y) {
  # Map to CMP parameter space
  Xbeta <- X %*% beta
  Zgama <- Z %*% gama
  nu <- exp(Zgama)
  loglambda <- suppressWarnings(
    nu * log(exp(Xbeta) + (nu - 1) / (2 * nu))
  )
  # Get the normalizing constants
  logz <- compute_logz(loglambda, nu)
  # Calcula o negativo do log da função de verossimilhança
  ll <- sum(y * loglambda - nu * lfactorial(y) - logz)
  return(-ll)
}

#' @rdname llcmp-compute
#'
llcmp <- function(params, X, Z, y) {
  p <- ncol(X)
  q <- ncol(Z)
  llcmp_fixed(beta = params[1:p],
              gama = params[1:q + p],
              X = X, Z = Z, y = y)
}

#-----------------------------------------------------------------------
#' @title Maximize the COM-Poisson log-likelihood function
#' @param start Initial parameters
#' @param X Design matrix related to the (approximate) mean parameter
#'   \eqn{\mu = \exp(X \beta)}.
#' @param Z Design matrix related to the dispersion parameter \eqn{\nu =
#'   \exp(Z \gamma)}.
#' @param y Vector of observed count data.
#' @param method Argument passed to \code{\link[stats]{optim}(...)}.
#' @param lower Argument passed to \code{\link[stats]{optim}(...)}.
#' @param upper Argument passed to \code{\link[stats]{optim}(...)}.
#' @param hessian Argument passed to \code{\link[stats]{optim}(...)}.
#' @param control Argument passed to \code{\link[stats]{optim}(...)}.
#' @return A list with estimated parameters and hessian matrix
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats glm.fit poisson optim
#' @export
#'
cmp_fit <- function(start   = NULL,
                    X, Z, y,
                    method  = c("BFGS",
                                "Nelder-Mead",
                                "CG",
                                "L-BFGS-B",
                                "SANN"),
                    lower   = -Inf,
                    upper   = Inf,
                    hessian = TRUE,
                    control = list()) {
  #--------------------------------------------
  # Dimensions
  n <- length(y)
  p <- ncol(X)
  q <- ncol(Z)
  #-------------------------------------------
  # Initial values
  if (is.null(start)) {
    model <- glm.fit(x = X, y = y, family = poisson())
    start <- c("beta" = model$coefficients, rep(0, q))
    names(start)[1:q + p] <- paste0("gama.", colnames(Z))
  } else {
    if (is.null(names(start)))
      names(start)  <- c(paste0("beta.", colnames(X)),
                         paste0("gama.", colnames(Z)))
  }
  #------------------------------------------
  # Maximization
  method <- match.arg(method)
  out <- optim(par = start,
               fn = llcmp,
               method = method,
               lower = lower,
               upper = upper,
               hessian = hessian,
               control = control,
               X = X,
               Z = Z,
               y = y)
  return(out)
}

#-----------------------------------------------------------------------
#' @title Fitting COM-Poisson models with varying dispersion
#' @param formula An object of class "\code{\link[stats]{formula}}"
#'   describe the model for mean.
#' @param dformula An object of class "\code{\link[stats]{formula}}"
#'   describe the model for dispersion.
#' @param data An optional data frame containing the variables in the
#'   model. If not found in data, the variables are taken from
#'   \code{environment(formula)}.
#' @param ... Arguments to be used by \code{\link{cmp_fit}}.
#' @return An object of class \code{cmpreg}.
#' @author Eduardo Jr <edujrrib@gmail.com>
#' @importFrom stats model.frame model.matrix model.response
#' @export
#'
cmp <- function(formula, dformula = ~1, data, ...) {
  #--------------------------------------------
  if (missing(data))
    data <- environment(formula)
  #-------------------------------------------
  # Get matrices
  frame <- model.frame(formula, data)
  terms <- attr(frame, "terms")
  X <- model.matrix(terms, frame)
  Z <- model.matrix(dformula, data)
  y <- model.response(frame)
  if (dformula == ~ 1) colnames(Z) <- "log(nu)"
  #-------------------------------------------
  # Fit model
  details <- cmp_fit(X = X, Z = Z, y = y, ...)
  mean_coefficients <- details$par[1:ncol(X)]
  disp_coefficients <- details$par[1:ncol(Z) + ncol(X)]
  names(mean_coefficients) <- colnames(X)
  names(disp_coefficients) <- colnames(Z)
  vcov <- NULL
  if ("hessian" %in% names(details))
    vcov <- solve(details$hessian)
  #--------------------------------------------
  # Fitted values
  fitted_mean <- exp(X %*% mean_coefficients)
  fitted_disp <- exp(Z %*% disp_coefficients)
  #--------------------------------------------
  # Output
  out <- list(call = match.call(),
              formula = formula,
              dformula = dformula,
              nobs = length(y),
              df.residual = length(y) - length(details$par),
              details = details,
              loglik = -details$value,
              vcov = vcov,
              coefficients = details$par,
              mean_coefficients = mean_coefficients,
              disp_coefficients = disp_coefficients,
              fitted_mean = c(fitted_mean),
              fitted_disp = c(fitted_disp),
              data = list(X = X, Z = Z, y = y))
  class(out) <- "cmpreg"
  return(out)
}
