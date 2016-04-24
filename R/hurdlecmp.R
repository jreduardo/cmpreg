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
#' @param offc Um vetor de valores a serem adicionados ao preditor
#'     linear da contagem \eqn{y > 0}
#' @param offz Um vetor de valores a serem adicionados ao preditor
#'     linear da contagem \eqn{y = 0}
#' @return O valor da log-verossimilhança do modelo
#'     Hurdle Conway-Maxwell-Poisson com os parâmetros e dados informados
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @seealso \code{\link[tccPackage]{hurdlecmp}}

llhurdle <- function(params, y, Xc, Xz, offc = NULL, offz = NULL) {
    zeropars <- params[grep("zero", names(params))]
    countpars <- params[grep("count", names(params))]
    ##-------------------------------------------
    ## Para a porção em zero
    yz <- ifelse(y == 0, 0, 1)
        if (is.null(offz)) offz <- 0
    Xbz <- Xz %*% zeropars + offz
    muz <- plogis(Xbz)
    llz <- sum(yz * log(muz) + (1 - yz) * log(1 - muz))
    ##-------------------------------------------
    ## Para contagem
    betas <- countpars[-1]
    phi <- countpars[1]
    yc <- y[y > 0]
    if (is.null(offc)) {
        offc <- 0
    } else
        offz <- offz[y > 0, ]
    Xc <- Xc[y > 0, ]
    Xbc <- Xc %*% betas + offc
    kernel <- sum(yc * Xbc - exp(phi) * lfactorial(yc))
    ## Obtendo a constante normatizadora Z.
    ## WARNING: Verificar a qtde de termos para a soma infinita
    i <- 0:300
    zs <- sapply(Xbc, function(lam)
        sum(exp(i * lam - exp(phi) * lfactorial(i))))
    Z <- sum(log(zs))
    llc <- kernel - Z - sum(log(1 - 1/zs))
    ##-------------------------------------------
    ## Verossimilhança combinada
    ll <- llz + llc
    return(ll)
}

#' @title Estimação do Modelo Hurdle Conway-Maxwell-Poisson
#' @description Estima os parâmetros de um modelo Hurdle COM-Poisson sob
#'     a otimização da função de log-verossimilhança. A sintaxe
#'     assemelha-se com a função \code{\link[pscl]{hurdle}} (Hurdle
#'     Models for Count Data Regression).
#' @param formula Um objeto da classe \code{\link{formula}}. Se
#'     necessária a inclusão de \emph{offset} deve-se indicá-lo como
#'     \code{\link{offset}}. O preditor linear do modelo para contagens
#'     0 e para contagens acima de zero deve ser separado pelo operador
#'     \code{|}, por exemplo se for descrito \code{y ~ x1 | x2} o
#'     preditor do modelo Binomial (para contagens 0) considerará apenas
#'     \code{x2} e o preditor para o modelo COM-Poisson somente o
#'     \code{x1}. Se não for especificado o operador \code{|} o mesmo
#'     proditor será considerado para ambas as porções do modelo.
#' @param data Um objeto de classe \code{data.frame}, cujo contém as
#'     variáveis descritas na \code{formula}
#' @param ... Argumentos opcionais do framework de maximização numérica
#'     \code{\link{optim}}
#' @return Uma lista de componentes do ajuste. Objeto de classe
#'     \code{hurdlecmp} cujo funções métodos foram implementadas.
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export

hurdlecmp <- function(formula, data, ...) {
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
    offz <- model.offset(framez)
    ## Para contagem (y > 0)
    framec <- model.frame(ffc, data)
    termsc <- attr(framec, "terms")
    Xc <- model.matrix(termsc, framec)
    offc <- model.offset(framec)
    ## Resposta
    y <- model.response(framec)
    ##-------------------------------------------
    ## Funcao de log-verossimilhanca
    fn <- function(params) {
         - llhurdle(params = params, y = y, Xc = Xc, Xz = Xz,
                 offz = offz, offc = offc)
    }
    ##-------------------------------------------
    ## Start - parametros iniciais
    startz <- glm.fit(Xz, factor(y > 0), family = binomial(),
                      offset = offz)$coef
    startc <- c(phi = 0,
                glm.fit(Xc, y, family = poisson(), offset = offc)$coef)
    ##-------------------------------------------
    ## Otimização
    opt <- optim(par = c(zero = startz, count = startc), fn = fn,
                 method = "BFGS", hessian = TRUE, ...)
    zeropars <- opt$par[grep("zero", names(opt$par))]
    countpars <- opt$par[grep("count", names(opt$par))]
    ##-------------------------------------------
    fit <- list(
        call = match.call(),
        form = formula,
        terms = list(zero = termsz, count = termsc),
        data = list(y = y, Xz = Xz, Xc = Xc, offz = offz, offc = offc),
        nobs = length(y),
        df = length(y) - length(opt$par),
        phi = opt$par["count.phi"],
        betas = countpars[-1],
        alphas = zeropars,
        logLik = -opt$value,
        niter = opt$counts,
        convergence = opt$convergence,
        hessian = opt$hessian
    )
    class(fit) <- "hurdlecmp"
    return(fit)
}

