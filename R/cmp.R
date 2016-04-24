#' @title Avaliação da Constante Normalizadora
#' @description Calcula o valor da constante de normalização do modelo
#'     COM-Poisson definida por: \deqn{Z = \sum
#'     \frac{\lambda^i}{(i!)^\nu}}. Para melhoria dos métodos de
#'     estimação o parâmetro \eqn{\nu} é tomado como \eqn{\exp{\phi}} e
#'     todas as inferências são tomadas a partir de \eqn{\phi}
#' @param lambda Parâmetro \eqn{\lambda} da função de distribuição de
#'     probabilidades COM-Poisson
#' @param phi Parâmetro \eqn{\phi = \log{\nu}} da função de distribuição
#'     de probabilidades COM-Poisson
#' @param tol Critério de parada do algoritmo, representa o valor
#'     tolerado para a diferença do valor de \eqn{Z(\lambda, \phi)}
#'     entre duas iterações. O valor padrão é 1e-3
#' @param maxit Número máximo de iterações a serem realizadas pelo
#'     algoritmo. Se este número for atingido e o critério de tolerância
#'     não for atendido, uma mensagem de aviso será exibida 
#' @return O valor da constante de normalização, \eqn{Z(\lambda,
#'     \phi)} da distribuição COM-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @examples
#'
#' ## Verificação da função
#' data(cottonBolls)
#' formula <- ncap ~ est:(des + I(des^2))
#' data <- cottonBolls
#' 
#' frame <- model.frame(formula, data)
#' terms <- attr(frame, "terms")
#' y <- model.response(frame)
#' X <- model.matrix(terms, frame)
#' 
#' betas <- coef(glm(formula, data = data, family = poisson))
#' phi <- 0
#' Xb <- X %*% betas
#' kernel <- sum(y * Xb - exp(phi) * lfactorial(y))
#' 
#' ## --- Utilizando a soma até 150
#' i <- 0:150
#' zs <- sapply(Xb, function(lam)
#'     sum(exp(i * lam - exp(phi) * lfactorial(i))))
#' (Z <- sum(log(zs)))
#' (ll <- kernel - Z)
#' 
#' ## --- Utilizando a computez, com critério de parada
#' zs1 <- sapply(exp(Xb), tccPackage:::computez, phi = phi)
#' (Z1 <- sum(log(zs1)))
#' (ll1 <- kernel - Z1)
#'

computez <- function(lambda, phi, tol = 1e-3, maxit = 1e3) {
    nu <- exp(phi)
    z0 <- 1000
    z <- 0
    i <- 2
    while (abs(z - z0) > tol && i < maxit) {
        z0 <- z
        z <- z0 + exp(i * log(lambda) - nu * lfactorial(i))
        i <- i + 1
    }
    if (abs(z - z0) > tol && i == maxit) {
        diff <- abs(z - z0)
        aviso <- sprintf(paste0(
            "  Valor da constante de normaliza\\u00e7\\u00e3o",
            " n\\u00e3o convergiu!\\n",
            "    Estimada em %.3f com %i itera\\u00e7\\u00f5es\\n",
            "    Diferen\\u00e7a entre 2 \\u00faltimas",
            " itera\\u00e7\\u00f5es %.5f"),
            z, i, diff)
        warning(aviso)
    }
    return(z + lambda + 1)
}

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
#' @seealso \code{\link[tccPackage]{cmp}}

llcmp <- function(betas, phi, y, X, offset = NULL){
    nu <- exp(phi)
    if (is.null(offset)) offset <- 0
    Xb <- X %*% betas + offset
    kernel <- sum(y * Xb - nu * lfactorial(y))
    ## Obtendo a constante normatizadora Z.
    ## WARNING: Verificar a qtde de termos para a soma infinita
    i <- 1:150
    zs <- sapply(Xb, function(lam)
        sum(exp(i * lam - nu * lfactorial(i))))
    Z <- sum(log(zs + 1))
    ##-------------------------------------------
    ## Ainda não funciona bem, o algoritmo de otimização leva a valores
    ## absurdos de phi, inviabilizando o cálculo de Z
    ## zs <- sapply(exp(Xb), computez, phi = phi)
    ## Z <- sum(log(zs))
    ##-------------------------------------------
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
#'     \code{compois} cujo funções métodos foram implementadas.
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @examples
#' 
#' \dontrun{
#' #-------------------------------------------
#' # Conjunto 1
#' data(soyaBeans)
#' help(soyaBeans, h = "html")
#' str(soyaBeans)
#'
#' da <- transform(soyaBeans, K = factor(potassio), A = factor(agua))
#' (model <- cmp(nv ~ bloco + K * A, data = da))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 2
#' data(whiteFly)
#' help(whiteFly, h = "html")
#' str(whiteFly)
#'
#' da <- droplevels(subset(whiteFly, grepl("BRS", x = cult)))
#' (model <- cmp(ntot ~ bloco + cult * dias, data = da))
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
#' (model <- cmp(ovos ~ offset(log(off)) + periodo + box + luz,
#'                   data = da))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 4
#' data(cottonBolls)
#' help(cottonBolls, h = "html")
#' str(cottonBolls)
#'
#' (model <- cmp(ncap ~ est:(des + I(des^2)), data = cottonBolls))
#' logLik(model)
#'
#' #-------------------------------------------
#' # Conjunto 5
#' data(cottonBolls2)
#' help(cottonBolls2, h = "html")
#' str(cottonBolls2)
#'
#' da <- transform(cottonBolls2, dexp = dexp - mean(range(dexp)))
#' (model <- cmp(ncapu ~ dexp + I(dexp^2), data = da))
#' logLik(model)
#' }
#' @export
#'
cmp <- function(formula, data, ...) {
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
        - llcmp(params[-1], params[1], y = y, X = X, offset = off)
    }
    opt <- optim(c(phi, betas), fn = fn, method = "BFGS",
                 hessian = TRUE, ...)
    ##-------------------------------------------
    fit <- list(
        call = match.call(),
        form = formula,
        terms = terms,
        data = list(y = y, X = X, offset = off),
        nobs = length(y),
        df = length(y) - length(opt$par),
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
