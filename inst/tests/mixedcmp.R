## =====================================================================
## Poisson e COM-Poisson misto (sem covariaveis)
##                                                        Eduardo Junior
##                                                    edujrrib@gmail.com
##                                                            2016-05-13
## =====================================================================

##======================================================================
## Pacotes
library(latticeExtra)
library(bbmle)
library(lme4)

##======================================================================
## Funções de verossimilhanca e densidade COM-Poisson
source("cmp-bbmle.R")

##======================================================================
## Define a função para aproximação da Laplace da integral
##     Fonte: material MCIE pág. 141
##======================================================================
laplace <- function(funcao, otimizador, n.dim, ...) {
    log_integral <- -sqrt(.Machine$double.xmax)
    inicial <- rep(0, n.dim)
    temp <- try(optim(inicial, funcao, ..., method = otimizador, 
                      hessian = TRUE, control = list(fnscale = -1)))
    if (class(temp) != "try-error") {
        ## Já deixa na escala do log, pois é assim que utilzamos no
        ## processo de otimização
        log_integral <- temp$value + 0.5 * log(2*pi) - 0.5 *
            determinant(-temp$hessian)$modulus
    }
    return(log_integral)
}

##======================================================================
## Simulando os dados
##======================================================================
simula.pois <- function(n, r, b0, sig, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    b <- rnorm(n, 0, sig)
    lambda <- exp(b0 + b)
    y <- rpois(n * r, lambda)
    da <- data.frame(y = y, id = 1:n)
    return(da = da[order(da$id), ])
}
da <- simula.pois(10, 10, 2, 0.8, seed = 2016)
(xy <- xyplot(y ~ id, groups = id, data = da,
       jitter.x = TRUE, pch = 19, grid = TRUE))

##======================================================================
## Modelo misto Poisson
##======================================================================

##-------------------------------------------
## log-verossimilhança por bloco
Li <- function(b, params, y) {
    lambda <- exp(params[1] + b)
    desvio <- exp(params[2])
    out <- sum(dpois(y, lambda, log = TRUE)) +
        dnorm(b, 0, desvio, log = TRUE)
    return(out)
}

##-------------------------------------------
## Função de verossimilhança marginal
vero.pois <- function(params, data, log = TRUE) {
    ## Necessita que os dados tenham uma coluna de nome (id)
    ## que identifica os blocos similares
    dados.id <- split(data, data$id)
    ## Calcula as integrais
    ints <- lapply(dados.id, function(dados) {
        laplace(Li, otimizador = "BFGS", n.dim = 1,
                params = params, y = dados$y)
    })
    ## Monta a log-verossimilhanca (soma das integrais dos grupos)
    ll <- sum(unlist(ints))
    ## print(params)
    return(-ll)
}

## Ajuste do modelo
parnames(vero.pois) <- c("b0", "lsig")
model.pois <- mle2(vero.pois, start = c("b0" = 1, "lsig" = 1),
           data = list(data = da))
model.pois

## Via glmer
model.pois2 <- glmer(y ~ 1|id, data = da, family = poisson)
model.pois2

## Comparando os resultados da implementação vero.pois com a glmer
coef.pois <- c(coef(model.pois)[1], exp(coef(model.pois)[2]))
coef.pois2 <- c(fixef(model.pois2), sqrt(VarCorr(model.pois2)$id))

cbind("vero.pois" = coef.pois, "glmer" = coef.pois2)
cbind("vero.pois" = logLik(model.pois), "gmler" = logLik(model.pois2))

##----------------------------------------------------------------------
## Valores preditos

## Obtém os efeitos aleatórios
ranef <- sapply(split(da, da$id), function(dados) {
    optim(0, Li, method = "BFGS", control = list(fnscale = -1),
          params = coef(model.pois), y = dados$y)$par
})
densityplot(~ranef)
## cbind(ranef(model.pois2)$id, ranef)

## Obtém os valores preditos
eta <- coef(model.pois)[1] + ranef
mu <- exp(eta)
## cbind(unique(predict(model.pois2, type = "response")), mu)

(xy2 <- update(xy, alpha = 0.5) +
     layer(panel.points(x = seq_along(mu) - 0.1, y = mu, col = 1,
                        pch = 15, cex = 1.2)))

##======================================================================
## Modelo misto COM-Poisson
##======================================================================

##-------------------------------------------
## log-verossimilhança por bloco
Li <- function(b, params, y) {
    phi <- params[3]
    loglambda <- params[1] + b
    desvio <- exp(params[2])
    out <- sum(dcmp(y, loglambda, phi, sumto = 50, log = TRUE)) +
        dnorm(b, 0, desvio, log = TRUE)
    return(out)
}

##-------------------------------------------
## log-verossimilhança marginal
vero.cmp <- function(params, data, log = TRUE) {
    ## Necessita que os dados tenham uma coluna de nome (id)
    ## que identifica os blocos similares
    dados.id <- split(data, data$id)
    ## Calcula as integrais
    ints <- lapply(dados.id, function(dados) {
        laplace(Li, otimizador = "BFGS", n.dim = 1,
                params = params, y = dados$y)
    })
    ## Monta a log-verossimilhanca (soma das integrais dos grupos)
    ll <- sum(unlist(ints))
    ## print(params)
    return(-ll)
}

## Ajuste do modelo
parnames(vero.cmp) <- c("b0", "esig", "phi")
model.cmp <- mle2(vero.cmp, start = c("b0" = 1, "esig" = 1, "phi" = 0),
                  data = list(data = da))
model.cmp

## Comparando os resultados com o modelo Poisson
coef.cmp <- coef(model.cmp)
coef.cmp[2] <- exp(coef.cmp[2])

cbind("vero.cmp" = coef.cmp, "vero.pois" = c(coef.pois, NA))
cbind("vero.cmp" = logLik(model.cmp), "vero.pois" = logLik(model.pois))

##----------------------------------------------------------------------
## Valores preditos

## Obtém os efeitos aleatórios
ranef <- sapply(split(da, da$id), function(dados) {
    optim(0, Li, method = "BFGS", control = list(fnscale = -1),
          params = coef(model.cmp), y = dados$y)$par
})
densityplot(~ranef)
## cbind(ranef(model.pois2)$id, ranef)

## Obtém os valores preditos
eta <- coef(model.cmp)[1] + ranef
mu2 <- sapply(eta, function(loglambda) {
    x <- 0:50
    px <- dcmp2(x, loglambda = loglambda, phi = coef.cmp[3],
                sumto = 50, log = FALSE)
    sum(x*px)    
})

update(xy2, alpha = 0.5) +
    layer(panel.points(x = seq_along(mu2) + 0.1, y = mu2, col = 1,
                       pch = 17, cex = 1.2))
## cbind(mu, mu2)

