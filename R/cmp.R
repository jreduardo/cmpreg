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
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @return O valor da log-verossimilhança do modelo
#'     Conway-Maxwell-Poisson com os parâmetros e dados informados
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @seealso \code{\link[tccPackage]{cmp}}

llcmp <- function(betas, phi, y, X, sumto = NULL){
    nu <- exp(phi)
    Xb <- X %*% betas
    kernel <- sum(y * Xb - nu * lfactorial(y))
    ## Obtendo a constante normatizadora Z.
    ## WARNING: Verificar a qtde de termos para a soma infinita
    if (is.null(sumto)) sumto <- max(y)^1.2
    i <- 1:sumto
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

#' @title Probabilidades do Modelo Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @description Calcula as probabilidades para uma variável aleatória
#'     distribuída conforme modelo COM-Poisson.
#'
#' \deqn{p(y,\lambda,\nu) =
#'     \frac{\lambda^y}{(y!)^\nu Z(\lambda, \nu)}
#' }
#'
#' em que \eqn{Z(\lambda, \nu)} é a constante de normalização definida
#'     por \eqn{\sum_{j=0}^{\infty} \frac{\lambda^j}{(j!)^\nu}}.  Nesta
#'     implementação o número de incrementos considerados para cálculo
#'     dessa constante é definido por \code{sumto}. \eqn{\lambda > 0} e
#'     \eqn{\nu \geq 0} são os parâmetros da distribuição.
#'
#' @param y Valor da variável de contagem.
#' @param lambda Valor do parâmetro \eqn{\lambda} da distribuição
#'     COM-Poisson.
#' @param nu Valor do parâmetro \eqn{\nu} da distribuição COM-Poisson.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @examples
#' dpois(5, lambda = 5)
#' dcmp(5, lambda = 5, nu = 1, sumto = 20)
#'
#' probs <- data.frame(y = 0:30)
#' within(probs, {
#'     py0 <- dpois(y, lambda = 15)
#'     py1 <- dcmp(y, lambda = 15, nu = 1, sumto = 50)
#'     py2 <- dcmp(y, lambda = 915, nu = 2.5, sumto = 50)
#'     py3 <- dcmp(y, lambda = 2.2, nu = 0.3, sumto = 50)
#'     plot(py0 ~ y, type = "h",
#'          ylim = c(0, max(c(py0, py2, py3))),
#'          ylab = expression(Pr(Y == y)))
#'     points(y + 0.1, py1, type = "h", col = 2)
#'     points(y - 0.3, py2, type = "h", col = 3)
#'     points(y + 0.3, py3, type = "h", col = 4)
#'     legend("topleft", bty = "n",
#'            col = c(1:4), lty = 1,
#'            legend = expression(
#'                Poisson(lambda == 15),
#'                CMP(lambda == 15, nu == 1),
#'                CMP(lambda == 915, nu == 2.5),
#'                CMP(lambda == 2.2, nu == 0.3)))
#' })

dcmp <- Vectorize(
    FUN = function(y, loglambda, phi, sumto, log = FALSE) {
        py <- sapply(y, function(yi) {
            llcmp(betas = loglambda, phi = phi, y = yi, X = 1,
                  sumto = sumto)
        })
        if (!log)
            py <- exp(py)
        return(py)
    }, vectorize.args = c("y", "loglambda", "phi"))

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
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @param ... Argumentos opcionais do framework de maximização numérica
#'     \code{\link{optim}}
#' @return Uma lista de componentes do ajuste. Objeto de classe
#'     \code{compois} cujo funções métodos foram implementadas.
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

cmp <- function(formula, data, sumto = NULL, ...) {
    ##-------------------------------------------
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    ##
    y <- model.response(frame)
    X <- model.matrix(terms, frame)
    off <- model.offset(frame)
    if(!is.null(model.offset(frame))) {
        stop("Este modelo ainda nao suporta offset")
    }
    ##
    m0 <- glm(formula, data = data, family = poisson)
    betas.init <- m0$coefficients
    phi.init <- -log(sum(resid(m0, type = "pearson")^2)/m0$df.residual)
    ##
    fn <- function(params) {
        - llcmp(params[-1], params[1], y = y, X = X, sumto = sumto)
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
