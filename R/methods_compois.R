#' @title Log-Verossimilhança para modelos Conway-Maxwell-Poisson
#' @name logLik.compois
#' @description Exibe a log-verossimilhança calculada para o modelo
#'     Conway-Maxwell-Poisson ajustado em objetos de classe
#'     \code{compois}.
#' @param object O objeto produzido por \code{\link{glm_cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @return Um objeto de classe \code{logLik}
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

logLik.compois <- function(object, ...) {
    ll <- object$logLik
    attr(ll, "df") <- object$df
    attr(ll, "nobs") <- object$nobs
    class(ll) <- "logLik"
    return(ll)
}


#' @title Exibição do Ajuste do Modelo Conway-Maxwell-Poisson
#' @name print.compois
#' @description Exibe o ajuste do modelo Conway-Maxwell-Poisson para
#'     objetos de classe \code{compois}
#' @param x O objeto produzido por \code{\link{glm_cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

print.compois <- function(x, ...) {
    cat("\nCall:\t")
    print(x$call)
    cat("\nCoefficients:", sep = "\n")
    print(x$betas)
    cat("\nDegrees of Freedom:", x$nobs - 1,
        "Total (i.e. Null); ", x$nobs - x$df, "Residual")
    cat("\nNull Deviance:\t")
    cat("\nResidual Deviance:\t", NULL)
    cat("\tAIC: ")
    cat(AIC(x), sep = "\n")
    invisible(x)
}
