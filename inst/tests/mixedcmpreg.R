## =====================================================================
## Poisson e COM-Poisson misto (com covariaveis)
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
Li <- function(b, params, y, X, Z) {
    ##-------------------------------------------
    nZ <- ncol(Z)
    nX <- ncol(X)
    lsigmas <- params[1:nZ]
    betas <- params[(nZ + 1):(nX + nZ)]
    ##-------------------------------------------
    lambda <- exp(X %*% betas + Z %*% b)
    ## Aqui fazemos efeitos aleatórios independentes, assim poder-se-ia
    ## somar dnorms porém para deixar de mais fácil generalização
    ## utiliza-se a dmvnorm (normal multivariada)
    sigma <- diag(nZ)
    diag(sigma) <- exp(lsigmas)^2
    out <- sum(dpois(y, lambda, log = TRUE)) +
        mvtnorm::dmvnorm(b, rep(0, nZ), sigma, log = TRUE)
    return(out)
}

##-------------------------------------------
## Função de verossimilhança marginal
vero.pois <- function(params, formula, data, log = TRUE) {
    ## Necessita que os dados tenham uma coluna de nome (id)
    ## que identifica os blocos similares
    ##-------------------------------------------
    ## Separa os dados nos blocos 
    dados.id <- split(data, data$id)
    ##-------------------------------------------
    ## Constrói as fórmulas para os efeitos fixos e aleatórios
    form.X <- lme4::nobars(formula)
    form.Z <- lme4::findbars(formula)
    ## Calcula as integrais
    ints <- lapply(dados.id, function(dados) {
        mf <- model.frame(lme4::subbars(formula), data = dados)
        Xi <- model.matrix(form.X, dados)
        Zi <- t(as.matrix(lme4::mkReTrms(form.Z, mf)$Zt))
        laplace(Li, otimizador = "BFGS", n.dim = 1,
                params = params, y = dados$y, X = Xi, Z = Zi)
    })
    ## Monta a log-verossimilhanca (soma das integrais dos grupos)
    ll <- sum(unlist(ints))
    ## print(params)
    return(-ll)
}

## Ajuste do modelo
parnames(vero.pois) <- c("lsig", "b0")
model.pois <- mle2(vero.pois, start = c("lsig" = 1, "b0" = 1),
                   data = list(data = da, formula = y ~ 1 | id))
model.pois

## Via glmer
model.pois2 <- glmer(y ~ 1|id, data = da, family = poisson)
model.pois2

## Comparando os resultados da implementação vero.pois com a glmer
coef.pois <- c(exp(coef(model.pois)[1]), coef(model.pois)[2])
coef.pois2 <- c(sqrt(VarCorr(model.pois2)$id), fixef(model.pois2))

cbind("vero.pois" = coef.pois, "glmer" = coef.pois2)
cbind("vero.pois" = logLik(model.pois), "gmler" = logLik(model.pois2))

##----------------------------------------------------------------------
## Valores preditos

## Obtém os efeitos aleatórios
ranef <- sapply(split(da, da$id), function(dados) {
    formula <- y ~ 1|id
    form.X <- lme4::nobars(formula)
    form.Z <- lme4::findbars(formula)
    mf <- model.frame(lme4::subbars(formula), data = dados)
    Xi <- model.matrix(form.X, dados)
    Zi <- t(as.matrix(lme4::mkReTrms(form.Z, mf)$Zt))
    optim(0, Li, method = "BFGS", control = list(fnscale = -1),
          params = coef(model.pois), y = dados$y, X = Xi, Z = Zi)$par
})

densityplot(~ranef)
## cbind(ranef(model.pois2)$id, ranef)

## Obtém os valores preditos
eta <- coef(model.pois)[2] + ranef
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
Li <- function(b, params, y, X, Z) {
    ##-------------------------------------------
    ## Modelo sem correlação entre os efeitos aleatórios
    ## **params**: vetor que representa os parametros do modelo. Deve
    ##             ser especificado em ordem:
    ##     - `phi`: parametro de precisão da distribuição COM-Poisson
    ##     - `lsigma`: parametros relacionados a matrix Z dos efeitos
    ##        aleatorios
    ##     - `betas`: parametros relacionados a matrix X dos efeitos
    ##        fixos
    ##-------------------------------------------
    nZ <- ncol(Z)
    nX <- ncol(X)
    phi <- params[1]
    lsigmas <- params[2:(1 + nZ)]
    betas <- params[(2 + nZ):(nZ + nX + 1)]
    ##-------------------------------------------
    loglambda <- X %*% betas + Z %*% b
    sigma <- diag(nZ)
    diag(sigma) <- exp(lsigmas)^2
    ##-------------------------------------------
    out <- sum(dcmp(y, loglambda, phi, sumto = 50, log = TRUE)) +
        mvtnorm::dmvnorm(b, rep(0, nZ), sigma, log = TRUE)
    return(out)
}

##-------------------------------------------
## log-verossimilhança marginal
vero.cmp <- function(params, formula, data, log = TRUE) {
    ## Necessita que os dados tenham uma coluna de nome (id)
    ## que identifica os blocos similares
    ##-------------------------------------------
    ## Separa os dados nos blocos 
    dados.id <- split(data, data$id)
    ##-------------------------------------------
    ## Constrói as fórmulas para os efeitos fixos e aleatórios
    form.X <- lme4::nobars(formula)
    form.Z <- lme4::findbars(formula)
    ## Calcula as integrais
    ints <- lapply(dados.id, function(dados) {
        mf <- model.frame(lme4::subbars(formula), data = dados)
        Xi <- model.matrix(form.X, dados)
        Zi <- t(as.matrix(lme4::mkReTrms(form.Z, mf)$Zt))
        laplace(Li, otimizador = "BFGS", n.dim = ncol(Zi),
                params = params, y = dados$y, X = Xi, Z = Zi)
    })
    ## Monta a log-verossimilhanca (soma das integrais dos grupos)
    ll <- sum(unlist(ints))
    ## print(params)
    return(-ll)
}

## Implementação coerente com o caso particular Poisson !
vero.cmp(c(0, coef(model.pois)), formula = y ~ 1|id, data = da)

## Ajuste do modelo
parnames(vero.cmp) <- c("phi", "lsig", "b0")
model.cmp <- mle2(vero.cmp, start = c("phi" = 0, "lsig" = 0, "b0" = 1),
                  data = list(data = da, formula = y ~ 1|id))
model.cmp

## Comparando os resultados com o modelo Poisson
coef.cmp <- coef(model.cmp)
coef.cmp[2] <- exp(coef.cmp[2])

cbind("vero.cmp" = coef.cmp, "vero.pois" = c(NA, coef.pois))
cbind("vero.cmp" = logLik(model.cmp), "vero.pois" = logLik(model.pois))

##----------------------------------------------------------------------
## Valores preditos

## Obtém os efeitos aleatórios
ranef <- sapply(split(da, da$id), function(dados) {
    formula <- y ~ 1|id
    form.X <- lme4::nobars(formula)
    form.Z <- lme4::findbars(formula)
    mf <- model.frame(lme4::subbars(formula), data = dados)
    Xi <- model.matrix(form.X, dados)
    Zi <- t(as.matrix(lme4::mkReTrms(form.Z, mf)$Zt))
    optim(0, Li, method = "BFGS", control = list(fnscale = -1),
          params = coef(model.cmp), y = dados$y, X = Xi, Z = Zi)$par
})
densityplot(~ranef)
## cbind(ranef(model.pois2)$id, ranef)

## Obtém os valores preditos
eta <- coef(model.cmp)[3] + ranef
mu2 <- sapply(eta, function(loglambda) {
    x <- 0:50
    px <- dcmp(x, loglambda = loglambda, phi = coef.cmp[1],
                sumto = 50, log = FALSE)
    sum(x*px)    
})

## Valores preditos sem o efeito aleatório
eta <- coef(model.cmp)[3]
mu3 <- sapply(eta, function(loglambda) {
    x <- 0:50
    px <- dcmp(x, loglambda = loglambda, phi = coef.cmp[1],
                sumto = 50, log = FALSE)
    sum(x*px)    
})

update(xy2, alpha = 0.5) +
    layer(panel.points(x = seq_along(mu2) + 0.1, y = mu2, col = 1,
                       pch = 17, cex = 1.2)) +
    layer(panel.abline(h = mu3))

## cbind(mu, mu2)
