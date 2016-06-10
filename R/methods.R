#' @importFrom stats qnorm binomial coef delete.response formula glm
#'     glm.fit model.frame model.matrix model.offset model.response
#'     optim plogis poisson resid terms
NULL

## @title Obtenção Pontual e Intervalar dos Preditores Lineares
## @author Walmes Zeviani, \email{walmes@@ufpr.br}.
## @description Função para obter o valores de \eqn{\eta = X\beta} que
##    é o preditor da parte de locação do modelo de regressão,
##    incluindo o intervalo de confiança para \eqn{\eta}, caso
##    corretamente especificado pelo argumento \code{qn}.
## @param V Matriz de variância e covariância das estimativas dos
##    parâmetros da parte de locação do modelo, necessário para
##    calcular a média.
## @param X Matriz de funções lineares para obter \eqn{\hat{eta} = X
##    \hat{\beta}}, em que \eqn{\beta} são os coeficientes da porção de
##    locação.
## @param b Estimativas dos parâmetros da parte de locação, ou seja,
##    \eqn{\hat{\beta}}.
## @param qn Quantis da distribuição normal padrão apropriados para um
##    intervalo de confiança conhecida.
## @return Um vetor se \code{length(qn) == 1} e uma matriz caso
##    contrário.
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

## @title Valores preditos pelo Modelo Conway-Maxwell-Poisson
## @author Eduardo E. R. Junior, \email{edujrrib@gmail.com} e Walmes
##    M. Zeviani
## @usage
## predict(object,
##        newdata,
##        type = c("link", "response"),
##        interval = c("none", "confidence"),
##        level = 0.95, ...)
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
predict.mle2 <- function(object, newdata,
                         type = c("link", "response"),
                         interval = c("none", "confidence"),
                         level = 0.95, ...) {
    ##-------------------------------------------
    cll <- as.character(object@call.orig)
    mdl <- grep(x = cll,
                pattern = "\\b(llcmp|llhurdle|llmixed)\\b",
                value = TRUE)
    if (!mdl %in% c("llcmp", "llhurdle", "llmixed")) {
        stop(paste("Func\\u00e3ao usada como `minuslogl`",
                   "n\\u00e3ao reconhecida."))
    } else {
        if (mdl %in% c("llhurdle", "llmixed")) {
            stop(paste("M\\u00e9todo de predi\\u00e7c\\u00e3ao",
                       "para ", mdl, "ainda n\\u00e3ao implementado."))
        }
    }
    ##----------------------------------------
    if (!missing(newdata)) {
        if (is.matrix(newdata)) {
            if (all(colnames(newdata) == names(coef(object)[-1]))) {
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
    interval <- match.arg(interval)
    qn <- -qnorm((1 - level[1])/2)
    qn <- switch(interval[1],
                 confidence = qn * c(lwr = -1, fit = 0, upr = 1),
                 none = qn * c(fit = 0))
    ##-------------------------------------------
    type <- match.arg(type)
    pred <-
        switch(mdl,
               "llcmp" = {
                   V <- object@vcov
                   Vc <- V[-1, -1] - V[-1, 1] %*%
                       solve(V[1, 1]) %*% V[1, -1]
                   eta <- cholV_eta(Vc, X,
                                    b = coef(object)[-1],
                                    qn = qn)
                   switch(type,
                          "link" = eta,
                          "response" = {
                              apply(as.matrix(eta),
                                    MARGIN = 2,
                                    FUN = calc_mean_cmp,
                                    phi = coef(object)[1],
                                    sumto = object@data$sumto)})
               })
    pred <- cbind(pred)
    colnames(pred) <- names(qn)
    return(pred)
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
