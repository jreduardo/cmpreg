## =====================================================================
## Estimação por Máxima verossimilhança via BBMLE
##                                                        Eduardo Junior
##                                                    edujrrib@gmail.com
##                                                            2016-05-13
## =====================================================================

## As funções abaixo têm cálculos idênticos aos já programados, contudo
## aqui a estrutura de input e output será voltada para utilização da
## função bbmle::mle2

##----------------------------------------------------------------------
## Log-verossimilhança do modelo COM-Poisson
##----------------------------------------------------------------------
llcmp <- function(params, y, X, sumto){
    ##-------------------------------------------
    betas <- params[-1]
    phi <- params[1]
    nu <- exp(phi)
    ##-------------------------------------------
    Xb <- X %*% betas
    ##-------------------------------------------
    ## Obtendo a constante normatizadora Z.
    i <- 0:sumto
    zs <- sapply(Xb, function(loglambda)
        sum(exp(i * loglambda - nu * lfactorial(i))))
    Z <- sum(log(zs))
    ##-------------------------------------------
    ll <- sum(y * Xb - nu * lfactorial(y)) - Z
    return(-ll)
}

##----------------------------------------------------------------------
## Densidade de probabilidade do modelo COM-Poisson
##----------------------------------------------------------------------
dcmp <- function (y, loglambda, phi, sumto, log = FALSE) {
    py <- sapply(y, function(yi) {
        -llcmp(c(phi, loglambda), y = yi, X = 1, sumto = sumto)
    })
    if(!log) py <- exp(py)
    return(py)
}

##----------------------------------------------------------------------
## Ajuste de modelos COM-Poisson, com layout similiar a glm
##----------------------------------------------------------------------
cmp <- function(formula, data, start = NULL, sumto = NULL, ...) {
    ##-------------------------------------------
    ## Constrói as matrizes do modelo
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    y <- model.response(frame)
    X <- model.matrix(terms, frame)
    ## off <- model.offset(frame)
    if (is.null(sumto)) sumto <- ceiling(max(y)^1.5)
    ##-------------------------------------------
    ## Define os chutes iniciais
    if (is.null(start)) {
        m0 <- glm.fit(x = X, y = y, family = poisson())
        start <- c("phi" = 0, m0$coefficients)
    }
    ##-------------------------------------------
    ## Nomeia os parâmetros da função para métodos bbmle
    bbmle::parnames(llcmp) <- names(start)
    model <- bbmle::mle2(llcmp, start = start,
                         data = list(y = y, X = X,
                                     sumto = sumto),
                         vecpar = TRUE, ...)
    return(model)
}
