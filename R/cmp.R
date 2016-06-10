#' @title Avaliação da Constante Normalizadora
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @description Calcula o valor da constante de normalização do modelo
#'     COM-Poisson definida por: \deqn{Z = \sum
#'     \frac{\lambda^i}{(i!)^\nu}}. Para melhoria dos métodos de
#'     estimação o parâmetro \eqn{\nu} é tomado como \eqn{\exp{\phi}} e
#'     todas as inferências são tomadas a partir de \eqn{\phi}
#' @param loglambda Valor de \eqn{\log(\lambda)}, sendo \eqn{\lambda} o
#'     parâmetro da função de distribuição de probabilidades COM-Poisson
#' @param phi Parâmetro \eqn{\phi = \log{\nu}} da função de distribuição
#'     de probabilidades COM-Poisson
#' @param tol Critério de parada do algoritmo, representa o valor
#'     tolerado para a diferença do valor de \eqn{Z(\lambda, \phi)}
#'     entre duas iterações. O valor padrão é 1e-3
#' @param maxit Número máximo de iterações a serem realizadas pelo
#'     algoritmo. Se este número for atingido e o critério de tolerância
#'     não for atendido, uma mensagem de aviso será exibida
#' @param incremento Número de incrementos da soma a serem considerados
#'     a cada iteração. Padrão definido como 10, ou seja, a cada
#'     iteração 10 incrementos são calculados.
#' @return O valor da constante de normalização, \eqn{Z(\lambda, \phi)}
#'     da distribuição COM-Poisson

computez <- function(loglambda, phi, tol = 1e-2, maxit = 500,
                     incremento = 10) {
    ##-------------------------------------------
    nu <- exp(phi)
    ##-------------------------------------------
    zg <- vector("list", maxit)
    t <- incremento
    i <- 1:t
    j <- 1
    ##
    zg[[j]] <- exp(i * loglambda - nu * lfactorial(i))
    ##
    while (abs(zg[[j]][t-1]) > tol && j < maxit) {
        i <- (i[t] + 1):(i[t] + t)
        j = j + 1
        zg[[j]] <- exp(i * loglambda - nu * lfactorial(i))
    }
    z <- unlist(zg)
    return(sum(z)+1)
}

#' @title Log-Verossimilhança do Modelo Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @description Calcula a log-verossimilhança de um modelo de regressão
#'     para o parâmetro \eqn{\lambda} considerando as respostas de
#'     contagem (y), condicionadas as suas covariáveis (X), distribuídas
#'     conforme modelo COM-Poisson.
#' @details A função de log-verossimilhança da COM-Poisson, na
#'     parametrização de modelo de regresssão é:
#'
#' \deqn{\ell(\beta, \nu, y) =
#'     \sum_{i=1}^{n} y_i \log(\lambda_i) - \nu \sum_{i=1}^{n}\log(y!)
#'     - \sum_{i=1}^{n} \log(Z(\lambda_i, \nu))}
#'
#' em que (i) \eqn{\lambda_i = \exp(X_i \beta)}, no modelo de regressão
#'     COM-Poisson um preditor linear é ligado à \eqn{\lambda} por meio
#'     da função de ligação log. Note que não estamos modelando
#'     diretamente a média, assim as estimativas dos parâmetros
#'     \eqn{\beta} não tem a mesma interpretação dos modelos Poisson,
#'     por exemplo. Contudo, os sinais desses parâmetros indicam efeitos
#'     de acréscimo ou descréscimo nas contagens médias.
#' (ii) \eqn{\nu} é o parâmetro de dispersão que indica equi, sub ou
#'     superdispersão das contagens y. Sendo \eqn{nu = 1} o caso de
#'     equidispersão, \eqn{0 \leq \nu < 1} superdispersão e \eqn{\nu >
#'     1} subdispersão. Vale ressaltar que a variância \eqn{V(Y)} não
#'     tem expressão fechada e não é definada unicamente por \eqn{\nu}.
#' (iii) \eqn{Z(\lambda_i, \nu)} é a constante de normalização definida
#'     por \deqn{\sum_{j=0}^{\infty} \frac{\lambda_i^j}{(j!)^\nu}}. Note
#'     que são cálculadas n constantes Z. Nesta implementação o número
#'     de incrementos considerados para cálculo dessas constantes é
#'     definido por \code{sumto}, o mesmo número de incrementos é
#'     considerado para o cálculo de todas as contantes. Uma verificação
#'     pós ajuste da escolha de \code{sumto} pode ser realizada a partir
#'     de \code{\link[MRDCr]{convergencez}}.
#'
#' Nesta parametrização o modelo COM-Poisson tem como casos particulares
#'     os modelos Poisson quando \eqn{\nu = 1}, Bernoulli quando
#'     \eqn{\nu \rightarrow \infty} (ou o modelo logístico considerando
#'     modelos de regressão) e Geométrico quando \eqn{\nu = 0} e
#'     \eqn{\lambda < 1}.
#'
#' Para que não seja necessário restringir o algoritmo de maximização da
#'     log-verossimilhança, a função foi implementada reparametrizando o
#'     parâmetro \eqn{\nu} para \eqn{\log(\phi)}. Assim o parâmetro
#'     estimado será \eqn{\phi} que tem suporte nos reais, assim como o
#'     vetor \eqn{\beta}.
#' @param params Um vetor de parâmetros do modelo COM-Poisson. O
#'     primeiro elemento desse vetor deve ser o parâmetro de dispersão
#'     do modelo, \eqn{\phi}, os restantes são os parâmetros
#'     \eqn{\beta}'s associados ao preditor linear em \eqn{\lambda}.
#' @param y Um vetor com variável dependente do modelo, resposta do tipo
#'     contagem.
#' @param X A matriz de delineamento correspondente ao modelo linear
#'     ligado à \eqn{\lambda} pela função de ligação log. A matriz do
#'     modelo pode ser construída com a função
#'     \code{\link[stats]{model.matrix}}.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo das constantes normalizadoras. Como padrão, para cálculo
#'     dessa constante faz-se uso de um processo iterativo, porém esse
#'     processo é demasiadamente demorado quando se ajusta os modelos
#'     via otimização numérica. Portanto indicar esse valor tornará o
#'     procedimento de estimação dos parâmetros mais veloz.
#' @return O negativo da log-verossimilhança do modelo
#'     Conway-Maxwell-Poisson com os parâmetros e dados informados.
#' @seealso \code{\link[bbmle]{mle2}}

llcmp <- function(params, y, X, sumto = NULL){
    ##-------------------------------------------
    betas <- params[-1]
    phi <- params[1]
    nu <- exp(phi)
    ##-------------------------------------------
    Xb <- X %*% betas
    ##-------------------------------------------
    ## Obtendo a constante normatizadora Z.
    if (is.null(sumto)) {
        zs <- sapply(Xb, function(loglambda)
            computez(loglambda, phi = phi, maxit = 1000))
    } else {
        i <- 0:sumto
        zs <- sapply(Xb, function(loglambda)
            sum(exp(i * loglambda - nu * lfactorial(i))))
    }
    ##-------------------------------------------
    ll <- sum(y * Xb - nu * lfactorial(y) - log(zs))
    return(-ll)
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
#' @param loglambda Valor do parâmetro tranformado \eqn{\log(\lambda)}
#'     da distribuição COM-Poisson.
#' @param phi Valor do parâmetro \eqn{\nu}, sob a reparametrização
#'     \eqn{\phi = \log(\nu)} da distribuição COM-Poisson.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @param log Argumento lógico qque indica se a densidade da COM-Poisson
#'     será retornada na escala logarítmica, seu valor padrão é
#'     \code{FALSE}.
#' @examples
#' dpois(5, lambda = 5)
#' dcmp(5, loglambda = log(5), phi = log(1), sumto = 20)
#'
#' probs <- data.frame(y = 0:30)
#' within(probs, {
#'     py0 <- dpois(y, lambda = 15)
#'     py1 <- dcmp(y, loglambda = log(15), phi = log(1), sumto = 50)
#'     py2 <- dcmp(y, loglambda = log(915), phi = log(2.5), sumto = 50)
#'     py3 <- dcmp(y, loglambda = log(2.2), phi = log(0.3), sumto = 50)
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
    FUN = function(y, loglambda, phi, sumto = NULL, log = FALSE) {
        py <- sapply(y, function(yi) {
            -llcmp(c(phi, loglambda), y = yi, X = 1, sumto = sumto)
        })
        if (!log)
            py <- exp(py)
        return(py)
    }, vectorize.args = c("y", "loglambda", "phi"))

#' @title Calcula o Valor Esperado para a Distribuição
#'     Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @description Função para calcular a média do tipo \eqn{E(Y) = \mu =
#'     \sum y\cdot \Pr(y)} para uma variável aleatória COM-Poisson a
#'     partir dos parâmetros \eqn{\lambda > 0} e \eqn{\nu \geq 0}.
#' @param loglambda Valor de \eqn{\log(\lambda)}, sendo \eqn{\lambda} o
#'     parâmetro da função de distribuição de probabilidades
#'     COM-Poisson.
#' @param phi Parâmetro \eqn{\phi = \log{\nu}} da função de distribuição
#'     de probabilidades COM-Poisson.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @param tol Tolerância para interromper a procura pelo valor de
#'     \code{ymax}, valor cuja probabilidade correspondente é inferior a
#'     \code{tol}, para valores os valores de \code{lambda} e
#'     \code{nu} informados.
#' @return Um vetor de tamanho igual ao do maior vetor, \code{lambda} ou
#'     \code{nu} com os valores correspondentes de \eqn{\mu}.

calc_mean_cmp <- function(loglambda, phi, sumto = NULL, tol = 1e-5) {
    ## Faz com que os parâmetros sejam vetores de mesmo tamanho.
    names(loglambda) <- NULL
    names(phi) <- NULL
    pars <- data.frame(loglambda = loglambda, phi = phi)
    ## Calcula o ymax usando mu + 5 (sqrt(sigma)) aproximados
    lambda <- exp(loglambda)
    nu <- exp(phi)
    approxmu <- lambda^(1/nu) - (nu - 1)/(2 * nu)
    sigma <- (1/nu) * approxmu
    ymax <- ceiling(max(approxmu + 5 * sqrt(sigma)))
    ## Agora verifica se a prob(ymax) é de fato pequena, se não, soma 1.
    loglambdamax <- max(pars$loglambda)
    phimin <- min(pars$phi)
    pmax <- dcmp(y = ymax, loglambda = loglambdamax,
                 phi = phimin, sumto = sumto)
    while (pmax > tol) {
        ymax <- ymax + 1
        pmax <- dcmp(y = ymax, loglambdamax, phimin, sumto = sumto)
    }
    ## Define o vetor onde avaliar a densidade COM-Poisson.
    y <- 1:ymax
    estmean <- mapply(loglambda = pars$loglambda,
                      phi = pars$phi,
                      MoreArgs = list(y = y, sumto = sumto),
                      FUN = function(loglambda, phi, y, sumto) {
                          py <- dcmp(y, loglambda, phi, sumto)
                          sum(y * py)
                      },
                      SIMPLIFY = TRUE)
    names(estmean) <- NULL
    return(estmean)
}

#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @title Calcula o Valor da Variância para a Distribuição
#'     Conway-Maxwell-Poisson
#' @description Função para calcular a variância do tipo \eqn{V(Y) =
#'     E(Y^2) - E^2(Y) = \sum y^2\cdot \Pr(y) - \left ( \sum y\cdot
#'     \Pr(y) \right )^2} para uma variável aleatória COM-Poisson a
#'     partir dos parâmetros \eqn{\lambda > 0} e \eqn{\nu \geq 0}.
#' @param loglambda Valor de \eqn{\log(\lambda)}, sendo \eqn{\lambda} o
#'     parâmetro da função de distribuição de probabilidades
#'     COM-Poisson.
#' @param phi Parâmetro \eqn{\phi = \log{\nu}} da função de distribuição
#'     de probabilidades COM-Poisson.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo da constante normalizadora Z.
#' @param tol Tolerância para interromper a procura pelo valor de
#'     \code{ymax}, valor cuja probabilidade correspondente é inferior a
#'     \code{tol}, para valores os valores de \code{lambda} e \code{nu}
#'     informados.
#' @return Um vetor de tamanho igual ao do maior vetor, \code{lambda} ou
#'     \code{nu} com os valores correspondentes de \eqn{\mu}.

calc_var_cmp <- function(loglambda, phi, sumto = NULL, tol = 1e-5) {
    ## Faz com que os parâmetros sejam vetores de mesmo tamanho.
    names(loglambda) <- NULL
    names(phi) <- NULL
    pars <- data.frame(loglambda = loglambda, phi = phi)
    ## Calcula o ymax usando mu + 5 (sqrt(sigma)) aproximados
    lambda <- exp(loglambda)
    nu <- exp(phi)
    approxmu <- lambda^(1/nu) - (nu - 1)/(2 * nu)
    sigma <- (1/nu) * approxmu
    ymax <- ceiling(max(approxmu + 5 * sqrt(sigma)))
    ## Agora verifica se a prob(ymax) é de fato pequena, se não, soma 1.
    loglambdamax <- max(pars$loglambda)
    phimin <- min(pars$phi)
    pmax <- dcmp(y = ymax, loglambda = loglambdamax,
                 phi = phimin, sumto = sumto)
    while (pmax > tol) {
        ymax <- ymax + 1
        pmax <- dcmp(y = ymax, loglambdamax, phimin, sumto = sumto)
    }
    ## Define o vetor onde avaliar a densidade COM-Poisson.
    y <- 1:ymax
    estvar <- mapply(loglambda = pars$loglambda,
                      phi = pars$phi,
                      MoreArgs = list(y = y, sumto = sumto),
                      FUN = function(loglambda, phi, y, sumto) {
                          py <- dcmp(y, loglambda, phi, sumto)
                          diff(c(sum(y * py)^2, sum(y^2 * py)))
                      },
                      SIMPLIFY = TRUE)
    names(estvar) <- NULL
    return(estvar)
}

#' @title Ajuste do Modelo Conway-Maxwell-Poisson
#' @author Eduardo E. R. Junior, \email{edujrrib@gmail.com}
#' @export
#' @description Estima os parâmetros de um modelo COM-Poisson pela
#'     otimização da função de log-verossimilhança definida em
#'     \code{\link{llcmp}}. A sintaxe assemelha-se com a função
#'     \code{\link{glm}} (Generalized Linear Models).
#' @param formula Um objeto da classe \code{\link{formula}}.
#' @param data Um objeto de classe \code{data.frame}, que contém as
#'     variáveis descritas na \code{formula}.
#' @param start Um vetor nomeado com os valores iniciais para os
#'     parâmetros do modelo necessários para o início do procedimento de
#'     estimação. Se \code{NULL} as estimativas de um modelo log-linear
#'     Poisson, com \eqn{\phi = 0}, são utilizadas como valores
#'     iniciais, pois uma chamada da \code{\link[stats]{glm.fit}} é
#'     feita internamente para obtê-los. O parâmetro \eqn{\phi} deve ser
#'     o primeiro elemento do vetor. Os restantes devem estar na
#'     correspondente às colunas da matriz gerada pelo argumento
#'     \code{formula}.
#' @param sumto Número de incrementos a serem considerados para a
#'     cálculo das constantes normalizadoras. Como padrão, para cálculo
#'     dessa constante faz-se uso de um processo iterativo, porém esse
#'     processo é demasiadamente demorado quando se ajusta os modelos
#'     via otimização numérica. Portanto indicar esse valor tornará o
#'     procedimento de estimação dos parâmetros mais veloz.
#' @param ... Argumentos opcionais do framework de maximização numérica
#'     \code{\link[bbmle]{mle2}}.
#' @return Um objeto de classe \code{mle2}, retornado da função de
#'     \code{\link[bbmle]{mle2}}, usada para ajuste de modelos por
#'     máxima verossimilhança.
#' @importFrom bbmle parnames mle2

cmp <- function(formula, data, start = NULL, sumto = NULL, ...) {
    ##-------------------------------------------
    ## Constrói as matrizes do modelo
    frame <- model.frame(formula, data)
    terms <- attr(frame, "terms")
    y <- model.response(frame)
    X <- model.matrix(terms, frame)
    if(!is.null(model.offset(frame))) {
        stop("Este modelo ainda nao suporta offset")
    }
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
                         data = list(y = y, X = X, sumto = sumto,
                                     terms = terms),
                         vecpar = TRUE, ...)
    return(model)
}
