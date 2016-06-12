#' @title Log-Verossimilhança do Modelo Hurdle Conway-Maxwell-Poisson
#' @description Calcula a log-verossimilhança de um modelo Hurdle
#'     COM-Poisson considerando os dados e as estimativas dos parâmetros
#'     informadas. As distribuições consideradas para a porção de
#'     valores nulos (\eqn{y = 0}) e para a contagem truncada \eqn{y >
#'     0} são Binomial e COM-Poisson respectivamente.
#' @details A função de log-verossimilhança do modelo toma a forma:
#'     \deqn{llcmp_{hurdle} = \log(\pi_i) I_{y = 0} + ( \log(1 - \pi_i)+
#'     y \log(\lambda) - e^{\phi} \log(y!) - \log(Z) - \log(1 - 1/Z))
#'     I_{y > 0}}, onde \eqn{Z = \sum \frac{\lambda^i}{(i!)^\nu}}
#' @param params Um vetor de estimativas para os parâmetros da
#'     distribuição Conway-Maxwell-Poisson truncada em zero. O verto
#'     deve ser nomeado com prefixo \code{count}, para os parâmetros do
#'     modelo COM-Poisson e \code{zero} para os que se referem ao modelo
#'     Binomial associado a contagens 0.
#' @param y Um vetor de contagens, considerado como variável resposta
#' @param Xc A matriz de delineamento do modelo para contagens \eqn{y >
#'     0}
#' @param Xz A matriz de delineamento do modelo para contagens nulas
#'     \eqn{y = 0}
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo das constantes normalizadoras. Como padrão, para cálculo
#'     dessa constante faz-se uso de um processo iterativo, porém esse
#'     processo é demasiadamente demorado quando se ajusta os modelos
#'     via otimização numérica. Portanto indicar esse valor tornará o
#'     procedimento de estimação dos parâmetros mais veloz.
#' @return O valor da log-verossimilhança do modelo
#'     Hurdle Conway-Maxwell-Poisson com os parâmetros e dados informados
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @seealso \code{\link[tccPackage]{hurdlecmp}}

llhurdle <- function(params, y, Xc, Xz, sumto = NULL) {
    ##-------------------------------------------
    ## Separa os parametros das duas partes do modelo
    zeropars <- params[grep("zero", names(params))]
    countpars <- params[grep("count", names(params))]
    ##-------------------------------------------
    ## Para a porção em zero
    yz <- ifelse(y == 0, 0, 1)
    Xbz <- Xz %*% zeropars
    muz <- plogis(Xbz)
    llz <- sum(yz * log(muz) + (1 - yz) * log(1 - muz))
    ##-------------------------------------------
    ## Para contagem
    betas <- countpars[-1]
    phi <- countpars[1]
    yc <- y[y > 0]
    Xc <- Xc[y > 0, ]
    Xbc <- Xc %*% betas
    ##-------------------------------------------
    ## Obtendo a constante normatizadora Z.
    if (is.null(sumto)) {
        zs <- sapply(Xbc, function(loglambda)
            computez(loglambda, phi = phi, maxit = 1000))
    } else {
        i <- 0:sumto
        zs <- sapply(Xbc, function(loglambda)
            sum(exp(i * loglambda - exp(phi) * lfactorial(i))))
    }
    llc <- sum(yc * Xbc - exp(phi) * lfactorial(yc) - log(zs) -
               log(1 - 1/zs))
    ##-------------------------------------------
    ## Verossimilhança combinada
    ll <- llz + llc
    return(-ll)
}

#' @title Ajuste de um Modelo Hurdle Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @description Estima os parâmetros de um modelo Hurdle COM-Poisson sob
#'     a otimização da função de log-verossimilhança. A sintaxe
#'     assemelha-se com a função \code{\link[pscl]{hurdle}} (Hurdle
#'     Models for Count Data Regression).
#' @param formula Um objeto da classe \code{\link{formula}}. O preditor
#'     linear do modelo para contagens 0 e para contagens acima de zero
#'     deve ser separado pelo operador \code{|}, e.g. esspecificando
#'     \code{y ~ x1 | x2} o preditor do modelo Binomial (para contagens
#'     0) considerará apenas \code{x2} e o preditor para o modelo
#'     COM-Poisson somente o \code{x1}. Se não for especificado o
#'     operador \code{|} o mesmo proditor será considerado para ambas as
#'     partes do modelo.
#' @param data Um objeto de classe \code{data.frame}, cujo contém as
#'     variáveis descritas na \code{formula}
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo das constantes normalizadoras. Como padrão, para cálculo
#'     dessa constante faz-se uso de um processo iterativo, porém esse
#'     processo é demasiadamente demorado quando se ajusta os modelos
#'     via otimização numérica. Portanto indicar esse valor tornará o
#'     procedimento de estimação dos parâmetros mais veloz.
#' @param ... Argumentos opcionais do framework de maximização numérica
#'     \code{\link[bbmle]{mle2}}
#' @return Um objeto de classe \code{mle2}, retornado da função de
#'     \code{\link[bbmle]{mle2}}, usada para ajuste de modelos por
#'     máxima verossimilhança.
#' @importFrom bbmle parnames mle2
#' @export

hurdlecmp <- function(formula, data, sumto = NULL, ...) {
    ##-------------------------------------------
    ## Dividindo a formula para os zeros (ffz) e para as contagens (ffc)
    if(length(formula[[3]]) > 1 && formula[[3]][[1]] == as.name("|")) {
        ff <- formula
        formula[[3]][1] <- call("+")
        ffc <- . ~ .
        ffz <- ~ .
        ffc[[3]] <- ff[[3]][[2]]
        ffc[[2]] <- ff[[2]]
        ffz[[3]] <- ff[[3]][[3]]
        ffz[[2]] <- NULL
    } else {
        ffz <- ffc <- ff <- formula
        ffz[[2]] <- NULL
    }
    if(inherits(try(terms(ffz), silent = TRUE), "try-error")) {
        ffz <- eval(parse(
            text = sprintf(paste("%s -", deparse(ffc[[2]])),
                           deparse(ffz))))
    }
    ##-------------------------------------------
    ## Criando as matrizes dos modelos
    ## Para zero (y = 0)
    framez <- model.frame(ffz, data)
    termsz <- attr(framez, "terms")
    Xz <- model.matrix(termsz, framez)
    ## Para contagem (y > 0)
    framec <- model.frame(ffc, data)
    termsc <- attr(framec, "terms")
    Xc <- model.matrix(termsc, framec)
    ##-------------------------------------------
    offz <- model.offset(framez)
    offc <- model.offset(framec)
    if(!is.null(offc) || !is.null(offz)) {
        stop("Este modelo ainda nao suporta offset")
    }
    ##-------------------------------------------
    ## Resposta
    y <- model.response(framec)
    ##-------------------------------------------
    ## Start - parametros iniciais
    startz <- glm.fit(Xz, factor(y > 0), family = binomial())$coef
    startc <- c(phi = 0, glm.fit(Xc, y, family = poisson())$coef)
    start <- c(count = startc, zero = startz)
    ##-------------------------------------------
    ## Otimização via bbmle
    bbmle::parnames(llhurdle) <- names(start)
    model <- bbmle::mle2(llhurdle, start = start,
                         data = list(y = y, Xc = Xc, Xz = Xz,
                                     sumto = sumto),
                         vecpar = TRUE, ...)
    return(model)
}
