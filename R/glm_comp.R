#' @title Log-Verossimilhança do Modelo Conway-Maxwell-Poisson
#' @description Calcula a log-verossimilhança de um modelo COM-Poisson
#'     considerando os dados e as estimativas dos parâmetros informadas.
#' @details A função de log-verossimilhança toma a forma: \deqn{-Z - y *
#'     \lambda - \nu \log{y!}}, onde \eqn{Z = \sum
#'     \frac{\lambda^i}{(i!)^\nu}}
#' @param betas Um vetor de estimativas para os parâmetros de regressão
#'     da distribuição Conway-Maxwell-Poisson.
#' @param phi Um valor estimado para o parâmetro de dispersão considerado
#'     na distrbuição Conway-Maxwell-Poisson \eqn{\phi = \log{\nu}}
#' @param y Um vetor de contagens, considerado como variável resposta
#' @param X A matriz de delineamento do modelo
#' @param offset Um vetor de valores a serem adicionados ao preditor
#'     linear
#' @return O valor da log-verossimilhança do modelo
#'     Conway-Maxwell-Poisson com os parâmetros e dados informados
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @seealso \code{\link[tccPackage]{glm_cmp}}

ll_cmp <- function(betas, phi, y, X, offset = NULL){
    nu <- exp(phi)
    Xb <- X %*% betas
    if (!is.null(offset)) 
        Xb <- Xb + offset
    kernel <- sum(y * Xb - nu * lfactorial(y))
    ## Obtendo a constante normatizadora Z.
    ## WARNING: Verificar a qtde de termos para a soma infinita
    i <- 1:150
    zs <- sapply(Xb, function(lam)
        sum(exp(i * lam - nu * lfactorial(i))))
    Z <- sum(log(zs + 1))
    ll <- kernel - Z
    return(ll)
}

#' @title Estimação do Modelo Conway-Maxwell-Poisson
#' @description Estima os parâmetros de um modelo COM-Poisson sob a
#'     otimização da função de log-verossimilhança. A sintaxe
#'     assemelha-se com a função \code{\link{glm}} (Generalized Linear
#'     Models).
#' @param formula Um objeto da classe \code{\link{formula}}. Se
#'     necessária a inclusão de \emph{offset} deve-se indicá-lo como
#'     \code{\link{offset}}
#' @param data Um objeto de classe \code{data.frame}, cujo contém as
#'     variáveis descritas na \code{formula}
#' @param ... Argumentos opcionais do framework de maximização numérica
#'     \code{\link{optim}}
#' @return Uma lista de componentes do ajuste. Objeto de classe
#'     \code{mcglm} cujo funções métodos foram implementadas.
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @examples
#'
#' #-------------------------------------------
#' # Conjunto 1
#' data(soyaBeans)
#' help(soyaBeans, h = "html")
#' str(soyaBeans)
#'
#' da <- transform(soyaBeans, K = factor(potassio), A = factor(agua))
#' (model <- glm_cmp(nv ~ bloco + K * A, data = da))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 2
#' data(whiteFly)
#' help(whiteFly, h = "html")
#' str(whiteFly)
#'
#' da <- droplevels(subset(whiteFly, grepl("BRS", x = cult)))
#' (model <- glm_cmp(ntot ~ bloco + cult * dias, data = da))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 3
#' data(eggs)
#' help(eggs, h = "html")
#' str(eggs)
#'
#' da <- aggregate(ovos ~ periodo + box + luz + gaiola,
#'                 data = eggs, FUN = sum)
#' da <- transform(da, off = 10 * 14)
#' (model <- glm_cmp(ovos ~ offset(log(off)) + periodo + box + luz,
#'                   data = da))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 4
#' data(cottonBolls)
#' help(cottonBolls, h = "html")
#' str(cottonBolls)
#'
#' (model <- glm_cmp(ncap ~ est:(des + I(des^2)), data = cottonBolls))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 5
#' data(cottonBolls2)
#' help(cottonBolls2, h = "html")
#' str(cottonBolls2)
#'
#' da <- transform(cottonBolls2, dexp = dexp - mean(range(dexp)))
#' (model <- glm_cmp(ncapu ~ dexp + I(dexp^2), data = da))
#' logLik(model)
#'
#' @export
glm_cmp <- function(formula, data, ...) {
    ##-------------------------------------------
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    ##
    y <- model.response(frame)
    X <- model.matrix(terms, frame)
    off <- model.offset(frame)
    ##
    poissonModel <- glm(formula, data = data,
                        family = poisson)
    betas <- poissonModel$coefficients
    phi <- 0
    ##
    fn <- function(params) {
        - ll_cmp(params[-1], params[1], y = y, X = X, offset = off)
    }
    opt <- optim(c(phi, betas), fn = fn, method = "BFGS",
                 hessian = TRUE, ...)
    ##-------------------------------------------
    fit <- list(
        call = match.call(),
        data = list(y = y, X = X, offset = off),
        nobs = length(y),
        df = length(opt$par),
        phi = opt$par[1],
        betas = opt$par[-1],
        logLik = -opt$value,
        niter = opt$count,
        convergence = opt$convergence,
        hessian = opt$hessian
    )
    class(fit) <- "compois"
    return(fit)
}
