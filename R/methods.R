#' @importFrom stats AIC binomial coef delete.response fitted glm
#'     glm.fit logLik model.frame model.matrix model.offset
#'     model.response optim pchisq plogis pnorm poisson predict
#'     printCoefmat qnorm resid residuals terms vcov
NULL

#' @title Log-Verossimilhança para modelos Conway-Maxwell-Poisson
#' @name logLik.compois
#' @description Exibe a log-verossimilhança calculada para o modelo
#'     Conway-Maxwell-Poisson ajustado em objetos de classe
#'     \code{compois}.
#' @param object O objeto produzido por \code{\link{cmp}}
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
#' @param x O objeto produzido por \code{\link{cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

print.compois <- function(x, ...) {
    cat("\nCall:\t")
    print(x$call)
    cat("\nCoefficients:", sep = "\n")
    print(x$betas)
    cat("\nDegrees of Freedom:", x$nobs - 2,
        "Total (i.e. Null); ", x$df, "Residual")
    cat("\nNull Deviance:", "TODO", "\t")
    cat("\nResidual Deviance:", "TODO", "\t", NULL)
    cat("\tAIC: ")
    cat(AIC(x), sep = "\n")
    invisible(x)
}

#' @title Matriz de Variância e Covariância do Modelo
#'     Conway-Maxwell-Poisson Ajustado
#' @name vcov.compois
#' @description Calcula a matriz de variância e covariância numérica de
#'     um modelo COM-Poisson
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @return Uma matriz nomeada, \code{matrix} com as variâncias e
#'     covariâncias entre os parâmetros estimados
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

vcov.compois <- function(object, ...) {
    return(solve(object$hessian))
}

#' @title Estimativas dos parâmetros
#' @name coef.compois
#' @description Extrai os coeficientes do modelo COM-Poisson ajustado
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @return Uma lista, \code{list} nomeada em que
#' \describe{
#' \item{betas}{É um vetor com o(s) coeficiente(s) estimados para o
#'     preditor linear.}
#' \item{nu}{É o valor estimado para o coeficiente \code{nu} do modelo
#'     COM-Poisson.}
#' }
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

coef.compois <- function(object, ...) {
    betas <- object$betas
    nu <- c("nu" = exp(object$phi))
    return(list(betas = betas, nu = nu))
}

#' @title Valores preditos pelo Modelo Conway-Maxwell-Poisson
#' @name predict.compois
#' @description Calcula os valores preditos pelo modelo COM-Poisson
#'     ajustado na escala da função de ligação, que neste caso não é
#'     expressamente definidas, mas é função dos parâmetros \eqn{\lambda}
#'     e \eqn{\nu} da distribuição, ou na escala da média da
#'     contagem. Intervalos de confiança para a média da contagem ou
#'     para os parâmetros \eqn{\lambda}'s estimados também estão
#'     implementados.
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param newdata Valores das covariáveis para predição
#' @param type A escala da predição. Assume valores \code{link} para a
#'     escala do parâmetro \eqn{\lambda} e \code{response} para escala
#'     da contagem
#' @param interval Para construção do intervalo de confiança. Assume
#'     valores \code{none} onde o intervalo não é calculado e
#'     \code{confidence} para construção do intervalo de confiança para
#'     a média da contagem sob o nível de confiança \code{level}.
#' @param level Nível de confiança para o intervalo predito
#' @param ... Argumentos adicionais para métodos S3
#' @return Um vetor, \code{vector} de valores preditos pelo modelo
#'     COM-Poisson, no caso de \code{interval = "none"} e uma matriz,
#'     \code{matrix} com os valores ajustados, inferiores e superiores
#'     (\emph{fit}, \emph{lwr} e \emph{upr} respectivamente) no caso de
#'     \code{interval = "confidence"}
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com} e Walmes
#'     M. Zeviani
#' @export

predict.compois <- function(object, newdata,
                            type = c("link", "response"),
                            interval = c("none", "confidence"),
                            level = 0.95, ...) {
    if (!missing(newdata)) {
        terms <- delete.response(object$terms)
	frame <- model.frame(terms, newdata)
        X <- model.matrix(terms, frame)
    } else {
        X <- object$data$X
    }
    lambdas <- c(X %*% object$betas)
    if (interval[1] == "none" && type[1] == "link" ) {
        pred <- lambdas
    }
    if (interval[1] == "none" && type[1] == "response" ) {
        pred <- sapply(lambdas, FUN = function(x) {
            compoisson::com.mean(exp(x), exp(object$phi))
        })
        names(pred) <- NULL
    }
    if (interval[1] == "confidence") {
        V <- vcov(object)
        Vc <- V[-1, -1] - V[-1, 1] %*% solve(V[1, 1]) %*% V[1, -1]
        U <- chol(Vc)
        se <- sqrt(apply(X %*% t(U),
                         MARGIN = 1, FUN = function(x) {
                             sum(x^2)
                         }))
        qn <- qnorm((1 - level)/2, lower.tail = FALSE) *
            c(fit = 0, lwr = -1, upr = 1)
        aux <- lambdas + outer(se, qn, FUN = "*")
        if (type[1] == "link") {
            pred <- aux
        }
        if (type[1] == "response") {
            ## Programar função que calcula a média. Para o caso de
            ## prever os dados utilizados para estimação pode-se evitar
            ## o calculo da constante de normalização novamente.
            pred <- apply(aux, MARGIN = 2, FUN = function(col) {
                sapply(col, FUN = function(x) {
                    compoisson::com.mean(exp(x), exp(object$phi))
                })
            })
        }
    }
    return(pred)
}

#' @title Valores Ajustados pelo Modelo Conway-Maxwell-Poisson
#' @name fitted.compois
#' @description Contagens ajustadas pelo Modelo Conway-Maxwell-Poisson
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @return Um vetor, \code{vector} de valores de contagem ajustados pelo
#'     modelo COM-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

fitted.compois <- function(object, ...) {
    fit <- predict(object, type = "response")
    return(fit)
}

#' @title Resíduos do Modelo Conway-Maxwell-Poisson Ajustado
#' @name residuals.compois
#' @description Calcula os resíduos do modelo ajustado
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param type Tipo de resíduo a ser calculado
#' @param ... Argumentos adicionais para métodos S3
#' @return Um vetor, \code{vector} de resíduos calculados sob o modelo
#'     COM-Poisson ajustado
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

residuals.compois <- function(object,
                              type = c("raw", "deviance", "pearson"),
                              ...) {
    if (type[1] == "raw") {
        return(object$data$y - fitted(object))
    } else {cat("\tTODO", sep = "\n")}
}


#' @title Resumo do Modelo Conway-Maxwell-Poisson Ajustado
#' @name summary.compois
#' @description Calcula as estatísticas resumo do modelo ajustado
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param ... Argumentos adicionais para métodos S3
#' @return Uma lista, \code{list} de componentes utilizados para resumir
#'     o ajuste do modelo. Objeto de classe \code{summary.compois}
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com} e Walmes
#'     M. Zeviani
#' @export

summary.compois <- function(object, ...) {
    call <- object$call
    betas <- object$betas
    phi <- object$phi
    est <- c(nu = exp(phi), betas)
    ##
    V <- vcov(object)
    sdt <- sqrt(diag(V))
    sdt[1] <- est[1] * sdt[1]
    zval <- est / sdt
    pval <- 2 * pnorm(-abs(zval))
    ##
    tab <- cbind("Estimate" = est, "Std. Error" = sdt,
                 "z value" = zval, "Pr(>|z|)" = pval)
    tabBetas <- tab[-1, ]
    tabNu <- cbind(t(tab[1, ]))
    ## Apenas para exibição do summary, alinhar tabelas
    nch <- max(nchar(rownames(tabBetas))) - 2
    rownames(tabNu) <- paste(paste(rep(" ", nch), collapse = ""), "nu")
    ##
    out <- list(
        call = call,
        residuals = summary(residuals(object))[-4],
        nu = est[1],
        tab = tab,
        tabBetas = tabBetas,
        tabNu = tabNu,
        df = object$df,
        nobs = object$nobs,
        aic = AIC(object),
        logLik = logLik(object)
    )
    class(out) <- "summary.compois"
    return(out)
}

#' @title Exibição do Resumo do Modelo Conway-Maxwell-Poisson Ajustado
#' @name print.summary.compois
#' @description Exibe as estatísticas resumo do modelo ajustado
#' @param x O objeto produzido por \code{\link{summary.compois}}
#' @param ... Argumentos adicionais para métodos S3
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

print.summary.compois <- function(x, ...) {
    cat("\nCall:\n")
    print(x$call)
    cat("\nDeviance Residuals:\n")
    print(x$residuals)
    cat("\nCoefficients:", sep = "\n")
    printCoefmat(x$tabBetas, P.values = TRUE, has.Pvalue = TRUE)
    cat("\nDispersion Coefficient:", sep = "\n")
    printCoefmat(x$tabNu, P.values = TRUE, has.Pvalue = TRUE,
                 signif.legend = FALSE)
    cat("\n(Dispersion parameter for COM-Poisson model is",
        paste0(round(x$nu, 3), ","), "Poisson occurs when it's 1)\n")
    cat("\n    Null deviance:", "TODO", "on", x$nobs - 2,
        "degrees of fredom")
    cat("\nResidual deviance:", "TODO", "on", x$df,
        "degrees of fredom")
    cat("\nAIC:", x$aic, "\t\tlogLik:", x$logLik, "\n")
    invisible(x)
}

#' @title Teste de Razão de Verossimilhanças
#' @name anova.compois
#' @description Realiza o teste de razão de verossimilhanças (TRV) para
#'     modelos COM-Poisson aninhados. Quando executado em um único
#'     modelo o TRV será feito com os modelos subjacentes ao informado,
#'     cujo serão removidos os efeitos preditores um a um. Caso
#'     executado com mais de um modelo o TRV será feito entre os modelos
#'     informados.
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param ... Demais modelos a serem submetidos ao teste de razão de
#'     verossimilhança
#' @return Uma lista, \code{list} de componentes utilizados para exibir
#'     o TRV. Objeto de classe \code{anova.compois}
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com} e Walmes
#'     M. Zeviani
#' @export

anova.compois <- function(object, ...){
    dots <- list(...)
    if (length(dots) == 0) {
        cat("\tTODO\t:Realizar ajustes sequencias pode demorar...",
            sep = "\n")
        ## Repartir o modelo em termos/parametros (via terms) e ajustar
        ## a sequencia de modelos acrescentando cada termo
        ##
        nmodels <- 1
        out <- list(nmodels = nmodels)
        attr(out, "type") <- "anodev"
    } else {
        modelsList <- list(object, ...)
        forms <- sapply(modelsList, function(x) deparse(x$form))
        nus <- sapply(modelsList, function(x) exp(x$phi))
        lls <- sapply(modelsList, function(x) logLik(x))
        nps <- sapply(modelsList, function(x) length(x$betas) + 1)
        cst <- 2 * diff(lls)
        pvs <- pchisq(abs(cst), df = abs(diff(nps)), lower.tail = FALSE)
        ##
        tab <- cbind(
            "nu" = nus,
            "Npar" = nps,
            "logLik" = lls,
            "LR stat" = c(NA, cst),
            "Npar diff" = c(NA, diff(nps)),
            "Pr(>Chisq)" = c(NA, pvs))
        ##
        nmodels <- length(modelsList)
        out <- list(tab = tab,
                    forms = forms,
                    nmodels = nmodels)
        attr(out, "type") <- "trv"
    }
    class(out) <- "anova.compois"
    return(out)
}

#' @title Exibição da ANOVA do Modelo Conway-Maxwell-Poisson
#' @name print.anova.compois
#' @description Exibe a análise de log-verossimilhança entre modelos
#'     COM-Poisson
#' @param x O objeto produzido por \code{\link{anova.compois}}
#' @param ... Argumentos adicionais para métodos S3
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

print.anova.compois <- function(x, ...) {
    if (attr(x, "type") == "anodev") {
        cat("Analysis of Deviance Table", sep = "\n")
    }
    if (attr(x, "type") == "trv") {
        cat("Likelihood Ratio Tests for COM-Poisson Models\n",
            sep = "\n")
        names <- paste0("model", 1:x$nmodels)
        names2 <- paste(names, ": ", x$forms)
        tab <- x$tab
        tab[,"nu"] <- round(tab[, "nu"], 2)
        rownames(tab) <- names
        ##
        sapply(names2, function(x) cat("  ", x, "\n"))
        cat("\n")
        printCoefmat(tab, P.values = TRUE, has.Pvalue = TRUE,
                     na.print = "", cs.ind = 2)
    }
    invisible(x)
}
