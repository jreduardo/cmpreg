#' @importFrom stats qnorm binomial coef delete.response formula glm
#'     glm.fit model.frame model.matrix model.offset model.response
#'     optim plogis poisson resid terms pchisq as.formula
NULL

#' @title Teste de Hipóteses para o parâmetro \eqn{phi}
#' @author Eduardo Junior, \email{edujrrib@gmail.com}
#' @description Realiza um teste de hipóteses bilateral para hipótese
#'     nula \eqn{H_0:\, \phi = 0}, ou seja, se o modelo COM-Poisson pode
#'     ser reduzido ao Poisson.
#' @param ... Uma sequência de modelos em que se realizará o teste.
#' @export

cmptest <- function(...) {
    cmp.list <- list(...)
    ##-------------------------------------------
    cls <- sapply(cmp.list, function(model) class(model))
    if (sum(!cls %in% "mle2") > 0 ) {
        stop(paste("Fun\\u00e7\\u00e3o dispon\\u00edvel apenas",
                   "para objectos de classe `mle2`"))
    }
    ##-------------------------------------------
    cll <- sapply(cmp.list, function(model)
        as.character(model@call.orig), simplify = FALSE)
    mdl <- sapply(cll, function(call) {
        grep(x = call, pattern = "\\b(llcmp|llhurdle|llmixed)\\b",
             value = TRUE)
    })
    cond <- sum(!mdl %in% c("llcmp", "llhurdle", "llmixed"))
    if (cond > 0) {
        stop(paste("Func\\u00e3o usada como `minuslogl`",
                   "n\\u00e3o reconhecida."))
    } else {
        cond <- sum(mdl %in% c("llhurdle", "llmixed"))
        if (cond > 0) {
            stop(paste("Teste para o par\\u00e2metro phi em",
                        paste(mdl, collapse = ", "),
                       "ainda n\\u00e3o implementado."))
        }
    }
    ##-------------------------------------------
    ## ## Fazer método print desta função
    ## msg1 <- paste0("  Teste de hip\\u00f3ts bilateral para o",
    ##                " par\\u00e2metro phi igual a zero\n")
    ## if (comment) cat(msg1, sep = "\n")
    ##-------------------------------------------
    llP <- sapply(cmp.list, function(model) {
        with(model@data, {
            fit <- glm.fit(X, y, family = poisson())
            loglik <- with(fit, rank - aic/2)
            loglik
        })
    })
    llC <- sapply(cmp.list, function(model) -model@min)
    ##-------------------------------------------
    trv <- 2 * (llC - llP)
    pvs <- pchisq(q = trv, df = 1, lower.tail = FALSE)
    ##-------------------------------------------
    phi <- sapply(cmp.list, function(model) model@fullcoef[1])
    ##-------------------------------------------
    out <- cbind(phi = phi, "P(>Chisq)" = pvs)
    rownames(out) <- rep("phi == 0", length(cmp.list))
    return(out)
}

#' @title Obtenção dos Efeitos Aleatórios no Modelo COM-Poisson Misto
#' @export
#' @author Eduardo Junior, \email{edujrrib@gmail.com}
#' @description Encontra os efeitos aleatórios que maximizam a
#'     log-verossimilhança do modelo  COM-Poisson de efeitos mistos.
#' @param object Uma sequência de modelos em que se realizará o teste.

mixedcmp.ranef <- function(object) {
    ##-------------------------------------------
    id <- as.character(object@data$form.Z[[1]][[3]])
    out <- with(object@data, {
        sapply(dados.id, function(dados) {
            ##-------------------------------------------
            mf <- model.frame(formula, data = dados)
            yi <- model.response(mf)
            Xi <- model.matrix(form.X, dados)
            Zi <- t(as.matrix(lme4::mkReTrms(form.Z, mf)$Zt))
            ##-------------------------------------------
            optim(0, llicmp, method = "BFGS",
                  control = list(fnscale = -1),
                  params = object@fullcoef,
                  y = yi,
                  X = Xi,
                  Z = Zi)$par
        })
    })
    return(out)
}

#' @title Obtenção Pontual e Intervalar dos Preditores Lineares
#' @author Walmes Zeviani, \email{walmes@ufpr.br}.
#' @description Função para obter o valores de \eqn{\eta = X\beta} que
#'    é o preditor da parte de locação do modelo de regressão,
#'    incluindo o intervalo de confiança para \eqn{\eta}, caso
#'    corretamente especificado pelo argumento \code{qn}.
#' @param V Matriz de variância e covariância das estimativas dos
#'    parâmetros da parte de locação do modelo, necessário para
#'    calcular a média.
#' @param X Matriz de funções lineares para obter \eqn{\hat{eta} = X
#'    \hat{\beta}}, em que \eqn{\beta} são os coeficientes da porção de
#'    locação.
#' @param b Estimativas dos parâmetros da parte de locação, ou seja,
#'    \eqn{\hat{\beta}}.
#' @param qn Quantis da distribuição normal padrão apropriados para um
#'    intervalo de confiança conhecida.
#' @return Um vetor se \code{length(qn) == 1} e uma matriz caso
#'    contrário.
cholV_eta <- function(V, X, b, qn) {
    eta <- c(X %*% b)
    if (length(qn) > 1) {
        U <- chol(V)
        stderr <- sqrt(apply(X %*% t(U),
                             MARGIN = 1,
                             FUN = function(x) {
                                 sum(x^2)
                             }))
        eta <- sweep(x = outer(stderr, qn, FUN = "*"),
                     MARGIN = 1,
                     STATS = eta,
                     FUN = "+")
    }
    return(eta)
}


#' @title Valores preditos pelo Modelo Conway-Maxwell-Poisson de efeitos
#'     mistos
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @description Calcula os valores preditos pelo modelo COM-Poisson
#'     ajustado na escala da função de ligação \eqn{\log(\lambda)}
#'     (default), ou na escala da média da contagem, que neste caso não
#'     é expressamente definida, mas é função dos parâmetros
#'     \eqn{\lambda} e \eqn{\nu} da distribuição.
#' @param object O objeto produzido por \code{\link{mixedcmp}}
#' @param newdata Valores das covariáveis para predição, em formato de
#'     \code{data.frame}.
#' @param type A escala da predição. Assume valores \code{link} para a
#'     escala logarítmica do parâmetro \eqn{\lambda} e \code{response}
#'     para escala da contagem.
#' @return Um vetor, \code{vector} de valores preditos pelo modelo
#'     COM-Poisson de efeitos mistos.

mixedcmp.predict <- function(object, newdata = NULL,
                             type = c("link", "response")) {
    ##-------------------------------------------
    type <- match.arg(type)
    params <- object@fullcoef
    ##-------------------------------------------
    id <- as.character(object@data$form.Z[[1]][[3]])
    nZ <- sum(grepl("lsigma", names(object@fullcoef)))
    nX <- length(params) - nZ - 1
    ##-------------------------------------------
    id <- as.character(object@data$form.Z[[1]][[3]])
    if (!is.null(newdata)) {
        dados.id <- split(newdata, newdata[, id])
    } else {
        dados.id <- object@data$dados.id
    }
    ##-------------------------------------------
    phi <- params[1]
    betas <- params[(2 + nZ):(nZ + nX + 1)]
    ranef <- mixedcmp.ranef(object)
    ##-------------------------------------------
    out <- vector("list", length = length(dados.id))
    for (i in 1:length(out)) {
        ##-------------------------------------------
        dados <- dados.id[[i]]
        ##-------------------------------------------
        Xi <- model.matrix(object@data$form.X[-2], dados)
        form.me <- as.formula(paste0("~0 + ", id))
        Zi <- model.matrix(form.me, dados)
        ##-------------------------------------------
        loglambdas <- Xi %*% betas + Zi %*% ranef
        out[[i]] <- switch(
            type,
            "link" = {
                loglambdas
            },
            "response" = {
                calc_mean_cmp(loglambda = loglambdas, phi = phi)
            }
        )
    }
    return(unlist(out))
}

#' @title Valores preditos pelo Modelo Conway-Maxwell-Poisson de efeitos
#'     fixos
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @description Calcula os valores preditos pelo modelo COM-Poisson
#'     ajustado na escala da função de ligação \eqn{\log(\lambda)}
#'     (default), ou na escala da média da contagem, que neste caso não
#'     é expressamente definidas, mas é função dos parâmetros
#'     \eqn{\lambda} e \eqn{\nu} da distribuição. Intervalos de
#'     confiança para a média da contagem ou para os parâmetros
#'     \eqn{\lambda}'s estimados também estão implementados.
#' @param object O objeto produzido por \code{\link{cmp}}
#' @param newdata Valores das covariáveis para predição, em formato de
#'     \code{data.frame} ou \code{matrix}, respeitando a nomenclatura
#'     das variáveis que compõem o modelo.
#' @param type A escala da predição. Assume valores \code{link} para a
#'     escala logarítmica do parâmetro \eqn{\lambda} e \code{response}
#'     para escala da contagem.
#' @param interval Para construção dos intervalos de confiança. Assume
#'     valores \code{none} onde o intervalo não é calculado e
#'     \code{confidence} para construção de intervalos de confiança para
#'     a média da contagem sob o nível de confiança \code{level}.
#' @param level Nível de confiança para o intervalo predito.
#' @return Um vetor, \code{vector} de valores preditos pelo modelo
#'     COM-Poisson, no caso de \code{interval = "none"} e uma matriz,
#'     \code{matrix} com os valores ajustados, inferiores e superiores
#'     (\emph{lwr}, \emph{fit} e \emph{upr} respectivamente) no caso de
#'     \code{interval = "confidence"}.

cmp.predict <- function(object, newdata,
                        type = c("link", "response"),
                        interval = c("none", "confidence"),
                        level = 0.95) {
    ##-------------------------------------------
    type <- match.arg(type)
    interval <- match.arg(interval)
    ##----------------------------------------
    if (!is.null(newdata)) {
        if (is.matrix(newdata)) {
            if (all(colnames(newdata) == names(object@fullcoef)[-1])) {
                X <- newdata
            } else {
                stop(paste("Nomes das colunas em `newdata` nao",
                           "bate com dos coeficientes."))
            }
        } else {
            if (is.data.frame(newdata)) {
                terms <- delete.response(object@data$terms)
                frame <- model.frame(terms, newdata)
                X <- model.matrix(terms, frame)
            } else {
                stop("`newdata` deve ser matriz ou data.frame.")
            }
        }
    } else {
        X <- object@data$X
    }
    ##-------------------------------------------
    qn <- -qnorm((1 - level[1])/2)
    qn <- switch(interval[1],
                 confidence = qn * c(lwr = -1, fit = 0, upr = 1),
                 none = qn * c(fit = 0))
    ##-------------------------------------------
    V <- object@vcov
    Vc <- V[-1, -1] - V[-1, 1] %*%
        solve(V[1, 1]) %*% V[1, -1]
    eta <- cholV_eta(Vc, X,
                     b = object@fullcoef[-1],
                     qn = qn)
    ##-------------------------------------------
    pred <- switch(type,
                   "link" = eta,
                   "response" = {
                       apply(as.matrix(eta),
                             MARGIN = 2,
                             FUN = calc_mean_cmp,
                             phi = object@fullcoef[1],
                             sumto = object@data$sumto)}
                   )
    pred <- cbind(pred)
    colnames(pred) <- names(qn)
    return(pred)
}

#' @title Valores preditos pelo Modelo Hurdle Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @description Calcula os valores preditos pelo modelo Hurdle
#'     COM-Poisson ajustado. São disponíveis os valores preditos na
#'     escala da contagem ou ainda as probabilidades preditas em um
#'     intervalo \code{at}. Adicionalmente disponibilizamos predições
#'     para os componentes zero e de contagem não nula. Nesses outros
#'     dois tipos de predição o valor do preditor linear \eqn{\eta} é
#'     escalonado para a média da distribuição considera, Binomial e
#'     COM-Poisson para as componentes zero e de contagem
#'     respectivamente.
#' @param object O objeto produzido por \code{\link{hurdlecmp}}
#' @param newdata Valores das covariáveis para predição, em formato de
#'     \code{data.frame} ou \code{matrix}, respeitando a nomenclatura
#'     das variáveis que compõem o modelo. Lembre-se de qque como agora
#'     temos dois processos atuantes, a definição dessa matriz deve
#'     contemplar os preditor para ambas as componentes.
#' @param type A escala da predição. Assume valores \code{response} para a
#'     média da contagem considerando o modelo combinado, \code{prob},
#'     para as probabilidades calculadas para \code{at} sob o modelo
#'     ajustado, \code{zero} para retornar o valor de \eqn{\hat{\pi}}, a
#'     probabilidade prevista pelo modelo Binomial para y>0 e
#'     \code{count} para a média da variável COM-Poisson desconsiderando
#'     a componente de contagens nulas.
#' @param at Intervalos de valores inteiros sob os quais calcular-se-á
#'     as probabilidades descritas pelo modelo. Esse argumento só tem
#'     efeito quando \code{type = "prob"}.
#' @return Os valores preditos pelo modelo, em formato de vetor para
#'     \code{type = c("link", "response", "count", "zero")} e de forma
#'     matricial para \code{type = "prob"}.

hurdle.predict <- function(object, newdata,
                           type = c("response", "count",
                                    "zero", "prob"),
                           at = NULL) {
    ##-------------------------------------------
    type <- match.arg(type)
    if (type == "prob" && is.null(at)) {
        at <- seq(min(object@data$y), max(object@data$y))
        message(paste0(
            "Intervalo para calculo das probabilidades",
            "nao definido. Probabilidades calculadas nos",
            "valores de y = ", at[1], " a y = ", rev(at)[1], "."))
    }
    ##-------------------------------------------
    coefnames <- names(object@fullcoef)
    indexc <- grepl(pattern = "^count.*", coefnames)
    ##----------------------------------------
    if (!is.null(newdata)) {
        if (is.matrix(newdata)) {
            if (all(colnames(newdata) == coefnames[-1])) {
                X <- newdata[, indexc]
                Z <- newdata[, !indexc]
            } else {
                stop(paste("Nomes das colunas em `newdata` nao",
                           "bate com dos coeficientes. Lembre-se de",
                           "colocar os prefixos `count.` e `zero.`"))
            }
        } else {
            if (is.data.frame(newdata)) {
                termsX <- delete.response(object@data$termsc)
                termsZ <- delete.response(object@data$termsz)
                X <- model.matrix(termsX, newdata)
                Z <- model.matrix(termsZ, newdata)
            } else {
                stop("`newdata` deve ser matriz ou data.frame.")
            }
        }
    } else {
        X <- object@data$Xc
        Z <- object@data$Xz
    }
    ##-------------------------------------------
    phi <- object@fullcoef[indexc][1]
    count.pars <- object@fullcoef[indexc][-1]
    zero.pars <- object@fullcoef[!indexc]
    ##-------------------------------------------
    etaX <- X %*% count.pars
    etaZ <- Z %*% zero.pars
    ##-------------------------------------------
    muz <- plogis(etaZ)
    muc <- calc_mean_cmp(etaX, phi = phi)
    p0.count <- dcmp(0, etaX, phi = phi)
    ##-------------------------------------------
    out <- switch(
        type,
        "response" = {
            c(muc * muz / (1 - p0.count))
        },
        "count" = {
            c(muc)
        },
        "zero" = {
            c(muz)
        },
        "prob" = {
            at <- unique(at)
            probs <- matrix(NA, nrow = length(muc), ncol = length(at))
            probs[, 1] <- (1 - muz)
            for (i in 2:length(at)) {
                ## probs[, i] <- (muz) * dcmp(at[i], etaX, phi = phi) *
                ## (1 - p0.count)^-1
                probs[, i] <- exp(
                    log(muz) +
                    dcmp(at[i], etaX, phi = phi, log = TRUE) -
                    log(1 - p0.count))
            }
            ## colnames(probs) <- paste0("y = ", at)
            ## rownames(probs) <- paste0("P(Y[", 1:length(muc), "] = y)")
            colnames(probs) <- at
            rownames(probs) <- 1:length(muc)
            probs
        }
    )
    return(out)
}

## @title Valores preditos pelo Modelo Conway-Maxwell-Poisson
## @author Eduardo E. R. Junior, \email{edujrrib@gmail.com} e Walmes
##    M. Zeviani
## @description Calcula os valores preditos pelo modelo COM-Poisson
##    ajustado na escala da função de ligação, que neste caso não é
##    expressamente definidas, mas é função dos parâmetros \eqn{\lambda}
##    e \eqn{\nu} da distribuição, ou na escala da média da
##    contagem. Intervalos de confiança para a média da contagem ou
##    para os parâmetros \eqn{\lambda}'s estimados também estão
##    implementados.
## @param object O objeto produzido por \code{\link{cmp}}
## @param newdata Valores das covariáveis para predição
## @param type A escala da predição. Assume valores \code{link} para a
##    escala do parâmetro \eqn{\lambda} e \code{response} para escala
##    da contagem
## @param interval Para construção do intervalo de confiança. Assume
##    valores \code{none} onde o intervalo não é calculado e
##    \code{confidence} para construção do intervalo de confiança para
##    a média da contagem sob o nível de confiança \code{level}.
## @param level Nível de confiança para o intervalo predito
## @param ... Argumentos adicionais para métodos S3
## @return Um vetor, \code{vector} de valores preditos pelo modelo
##    COM-Poisson, no caso de \code{interval = "none"} e uma matriz,
##    \code{matrix} com os valores ajustados, inferiores e superiores
##    (\emph{fit}, \emph{lwr} e \emph{upr} respectivamente) no caso de
##    \code{interval = "confidence"}

predict.mle2 <- function(object, newdata, ...) {
    ##-------------------------------------------
    cll <- as.character(object@call.orig)
    mdl <- grep(x = cll,
                pattern = "\\b(llcmp|llhurdle|llmixed)\\b",
                value = TRUE)
    if (!mdl %in% c("llcmp", "llhurdle", "llmixed")) {
        stop(paste("Func\\u00e3ao usada como `minuslogl`",
                   "n\\u00e3ao reconhecida."))
    }
    ##-------------------------------------------
    if (missing(newdata)) {
        newdata = NULL
    }
    ##-------------------------------------------
    out <- switch(mdl,
                  "llcmp" = cmp.predict(
                      object = object, newdata = newdata, ...),
                  "llhurdle" = hurdle.predict(
                      object = object, newdata = newdata, ...),
                  "llmixed" = mixedcmp.predict(
                      object = object, newdata = newdata, ...)
                  )
    return(out)
}


## @title Valores Ajustados pelo Modelo Conway-Maxwell-Poisson
## @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
## @usage
## fitted(object)
## @description Contagens ajustadas pelo Modelo Conway-Maxwell-Poisson
## @param object O objeto produzido por \code{\link{cmp}}
## @param ... Argumentos adicionais para métodos S3
## @return Um vetor, \code{vector} de valores de contagem ajustados pelo
##    modelo COM-Poisson
fitted.mle2 <- function(object, ...) {
    fit <- predict.mle2(object, type = "response")
    return(fit)
}

## @title Resíduos do Modelo Conway-Maxwell-Poisson Ajustado
## @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
## @description Calcula os resíduos do modelo ajustado
## @usage
## residuals(object,
##          type = c("raw", "pearson"))
## @param object O objeto produzido por \code{\link{cmp}}
## @param type Tipo de resíduo a ser calculado
## @param ... Argumentos adicionais para métodos S3
## @return Um vetor, \code{vector} de resíduos calculados sob o modelo
##    COM-Poisson ajustado
residuals.mle2 <- function(object,
                           type = c("raw", "pearson"),
                           ...) {
    type <- match.arg(type)
    raw <- object@data$y - fitted.mle2(object)
    switch(type,
           "raw" = {
               out <- raw
           },
           "pearson" = {
               betas <- object@fullcoef[-1]
               phi <- object@fullcoef[1]
               loglambda <- object@data$X %*% betas
               sumto <- object@data$sumto
               ##-------------------------------------------
               variancias <- calc_var_cmp(
                   loglambda = loglambda, phi = phi, sumto = sumto)
               ##-------------------------------------------
               residp <- raw/sqrt(variancias)
               out <- residp
           })
    return(out)
}

## @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
## @description Apenas para criar o método \code{predict} para objetos
##     retornados pela \code{\link[bbmle]{mle2}}.
methods::setMethod(f = "predict",
                   signature = "mle2",
                   definition = predict.mle2)

## @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
## @description Apenas para criar o método \code{residuals} para objetos
##     retornados pela \code{\link[bbmle]{mle2}}.
methods::setMethod(f = "residuals", signature = "mle2",
                   definition = residuals.mle2)
