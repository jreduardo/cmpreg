#' @title Integração Numérica por meio do método de Laplace
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com},
#'     implementação original em Ribeiro Jr et al. (2012).
#' @description Aproxima a integral de uma função pelo método de
#'     Laplace, ou seja, \deqn{\int f(x) e^{-x^2} dx =
#'     (2\pi)^{n/2}|{Q}|''(\hat{x})^{-1/2}e^{Q(\hat{x})}}.
#' @param funcao Uma função cujo qual se tem interesse na obtenção do
#'     valor de sua integral.
#' @param n.dim Número de dimensões a ser integral avaliada.
#' @param log Argumento lógico se a escala retorna é logarítmica ou
#'     não. Por default \code{log = TRUE}.
#' @param ... Argumentos adicionais passados a função
#'     \code{\link{optim}}.

laplace <- function(funcao, n.dim, log = TRUE, ...) {
    log_integral <- -sqrt(.Machine$double.xmax)
    inicial <- rep(0, n.dim)
    temp <- try(optim(inicial, funcao, ...,
        hessian = TRUE, control = list(fnscale = -1)))
    if (class(temp) != "try-error") {
        ## Já deixa na escala do log, pois é assim que utilizamos no
        ## processo de otimização
        log_integral <- temp$value + 0.5 * log(2 * pi) - 0.5 *
            determinant(-temp$hessian)$modulus
    }
    out <- log_integral
    if (!log)
        out <- exp(log_integral)
    return(out)
}

#' @title Log-verossimilhança por grupo do modelo COM-Poisson misto
#' @description Avalia a log-verossimilhança de cada grupo formado pelas
#'     categóricas da variável considerada como aleatória.
#' @param b Efeito aleatório nosiderada para avaliação da
#'     log-verossimilhança.
#' @param params Vetor de parâmetros do modelo. Deve seguir a ordenação
#'     dada abaixo: \itemize{
#'
#' \item{\code{phi}}{Parâmetro extra da distribuição COM-Poisson;}
#'
#' \item{\code{lsigma}}{Parâmetros relacionados a matriz Z dos efeitos;
#'     aleatórios. Representam as \eqn{\log(\sqrt(\textrm{variâncias}))}
#'     na diagonal da matriz de variâncias e covariâncias \eqn{Sigma} da
#'     Normal q-variada considerada para os efeitos aleatórios;}
#'
#' \item{\code{betas}}{Parâmetros relacionados aos efeitos fixos do
#'     modelo.}
#'
#' }
#' @param y Um vetor de contagens do grupo considerado.
#' @param X A matriz de delineamento dos efeitos fixos do modelo no
#'     grupo considerado.
#' @param Z A matriz de delineamento dos efeitos aleatórios do modelo no
#'     grupo considerado.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.

llicmp <- function(b, params, y, X, Z, sumto = NULL) {
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
    out <- sum(dcmp(y, loglambda, phi, sumto = sumto, log = TRUE)) +
        mvtnorm::dmvnorm(b, rep(0, nZ), sigma, log = TRUE)
    return(out)
}

#' @title Log-verossimilhança Marginal do Modelo COM-Poisson Misto
#' @description Avalia a log-verossimilhaça marginal de um modelo
#'     COM-Poisson misto. Aqui as log-verossimilhanças de cada grupo são
#'     integradas numericamente em relação aos efeitos aleatórios.
#' @param params Vetor de parâmetros do modelo. Deve seguir a ordenação
#'     dada abaixo: \itemize{
#'
#' \item{\code{phi}}{Parâmetro extra da distribuição COM-Poisson;}
#'
#' \item{\code{lsigma}}{Parâmetros relacionados a matriz Z dos efeitos;
#'     aleatórios. Representam as \eqn{\log(\sqrt(\textrm{variâncias}))}
#'     na diagonal da matriz de variâncias e covariâncias \eqn{Sigma} da
#'     Normal q-variada considerada para os efeitos aleatórios;}
#'
#' \item{\code{betas}}{Parâmetros relacionados aos efeitos fixos do
#'     modelo.}
#'
#' }
#' @param form Um objeto da classe \code{\link{formula}}. Essa fórmula
#'     substitue os "pipes", da fármula genérica para modelos mistos,
#'     pelo símbolo de adição, veja \code{\link[lme4]{subbars}}.
#' @param form.X Um objeto da classe \code{\link{formula}}. Essa fórmula
#'     deve representar somente os efeitos fixos do modelo, veja
#'     \code{\link[lme4]{nobars}}.
#' @param form.Z Um objeto da classe \code{\link{formula}}. Essa fórmula
#'     deve representar somente os efeitos aleatórios do modelo, veja
#'     \code{\link[lme4]{findbars}}.
#' @param dados.id Um objeto de classe \code{\link{list}}. Esa lista
#'     deve ser uma divisão dos dados, onde cada slot contém a
#'     participação comum a um grupo identificado por níveis
#'     considerados como efeitos aleatórios.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.

llmixedcmp <- function(params, form, form.X, form.Z, dados.id,
                       sumto = NULL) {
    ##-------------------------------------------
    ## Constrói as fórmulas para os efeitos fixos e aleatórios
    ## Calcula as integrais
    ints <- lapply(dados.id, function(dados) {
        mf <- model.frame(form, data = dados)
        y <- model.response(mf)
        Xi <- model.matrix(form.X, dados)
        Zi <- t(as.matrix(lme4::mkReTrms(form.Z, mf)$Zt))
        laplace(llicmp, method = "BFGS", n.dim = ncol(Zi),
                params = params, y = y, X = Xi, Z = Zi, sumto = sumto)
    })
    ## Monta a log-verossimilhanca (soma das integrais dos grupos)
    ll <- sum(unlist(ints))
    ## print(params)
    return(-ll)
}


#' @title Ajuste de um Modelo Misto Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @description Estima os parâmetros de um modelo misto COM-Poisson sob
#'     a otimização da função de log-verossimilhança. A sintaxe
#'     assemelha-se com a função \code{\link[lme4]{lmer}}, ou
#'     \code{\link[lme4]{glmer}} (Generalized Linear Mixed-Effects
#'     Models).
#' @param formula Um objeto da classe \code{\link{formula}}. A
#'     especificação dos efeitos aleatórios é feita neste argumento,
#'     e.g. \code{y ~ x + (x | id)} especifica que o modelo será um
#'     COM-Poisson de intercepto e inclinação aleatórios. Os efeitos
#'     considerados são independentes e normalmente distribuídos.
#' @param data Um objeto de classe \code{data.frame}, cujo contém as
#'     variáveis descritas na \code{formula}
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @param ... Argumentos opcionais do framework de maximização numérica
#'     \code{\link[bbmle]{mle2}}.
#' @return Um objeto de classe \code{mle2}, retornado da função de
#'     \code{\link[bbmle]{mle2}}, usada para ajuste de modelos por
#'     máxima verossimilhança.
#' @importFrom bbmle mle2 parnames

mixedcmp <- function(formula, data, sumto = NULL, ...) {
    ##-------------------------------------------
    form <- lme4::subbars(formula)
    form.X <- lme4::nobars(formula)
    form.Z <- lme4::findbars(formula)
    if (length(form.Z) > 1) {
        msg <- paste("Efeitos aleat\\u00f3rios para mais",
                     "de um grupo ainda n\\u00e3o implementado")
        stop(msg)
    }
    ##-------------------------------------------
    ## Separa os dados nos blocos
    id <- as.character(lme4::findbars(formula)[[1]][[3]])
    dados.id <- split(data, data[, id])
    ##-------------------------------------------
    ## Parâmetros iniciais
    nZ <- sum(!as.character(form.Z[[1]]) %in% c("1", "|"))
    m0 <- glm(form.X, data = data, family = poisson)
    phi.init <- -log(sum(resid(m0, type = "pearson")^2)/m0$df.residual)
    lsigma.init <- rep(0, nZ)
    beta.init <- coef(m0)
    ##-------------------------------------------
    names(phi.init) <- "phi"
    names(lsigma.init) <- paste0("lsigma", 0:(nZ-1))
    start <- c(phi.init, lsigma.init, beta.init)
    ##------------------------------------------
    bbmle::parnames(llmixedcmp) <- names(start)
    dataL <- list(form = form, form.X = form.X, form.Z = form.Z,
                  dados.id = dados.id, sumto = sumto)
    model <- bbmle::mle2(llmixedcmp, start = start, data = dataL)
    return(model)
}
