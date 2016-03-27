#' @name defoliation
#' @title Capulhos de Algodão em Função de Desfolha Artificial
#' @description Experimento conduzido sob delineamento interamente
#'     casualizado com 5 repetições em casa de vegetação com plantas de
#'     algodão \emph{Gossypium hirsutum} submetidas à diferentes níveis
#'     de desfolha artificial de remoção foliar, em combinação com o
#'     estágio fenológico no qual a desfolha foi aplicada. A unidade
#'     experimental foi um vaso com duas plantas onde avaliou-se o
#'     número de capulhos produzidos ao final da ciclo cultura.
#' @format Um code{data.frame} com 125 observações e 4 variáveis
#'
#' \describe{
#'
#' \item{\code{est}}{Um fator categórico com 5 níveis que representa o
#'     estágio fenológico da planto durante a aplicação da desfolha.}
#'
#' \item{\code{des}}{Um fator numérico com 5 níveis que representa o
#'     nível de desfolha artificial (percentual da área da folha}
#'     removida com tesoura) aplicada a todas as folhas na planta.
#'
#' \item{\code{rept}}{Inteiro que representa cada unidade experimental.}
#'
#' \item{\code{ncap}}{Inteiro que representa o número de capulhos de
#'     algodão produzidos ao final da ciclo cultura.}
#' }
#'
#' @docType data
#' @keywords subdispersão
#' @usage data(defoliation)
#'
#' @references Silva, A. M., Degrande, P. E., Suekane, R., Fernandes,
#'     M. G., & Zeviani, W. M. (2012). Impacto de diferentes níveis de
#'     desfolha artificial nos estádios fenológicos do
#'     algodoeiro. Revista de Ciências Agrárias, 35(1), 163–172.
#'
#' Zeviani, W. M., Ribeiro, P. J., Bonat, W. H., Shimakura, S. E., &
#'     Muniz, J. A. (2014). The Gamma-count distribution in the analysis
#'     of experimental underdispersed data. Journal of Applied
#'     Statistics, 41(12),
#'     1–11. http://doi.org/10.1080/02664763.2014.922168
#'
#' @examples
#' data(defoliation)
#'
#' library(lattice)
#'
#' xyplot(ncap ~ des | est,
#'        data = defoliation,
#'        layout = c(NA, 2),
#'        type = c("p", "smooth"),
#'        xlab = "Níveis de desfolha artificial",
#'        ylab = "Número de capulhos produzidos",
#'        xlim = extendrange(c(0:1), f = 0.15),
#'        jitter.x = TRUE,
#'        grid = TRUE)
#'
#' # Média e variância amostral para cada unidade experimental
#' (mv <- aggregate(ncap ~ est + des, data = defoliation,
#'                 FUN = function(x) c(mean = mean(x), var = var(x))))
#' xlim <- ylim <- extendrange(c(mv$ncap), f = 0.05)
#'
#' # Evidência de subdispersão
#' xyplot(ncap[, "var"] ~ ncap[, "mean"],
#'        data = mv,
#'        xlim = xlim,
#'        ylim = ylim,
#'        ylab = "Sample variance",
#'        xlab = "Sample mean",
#'        panel = function(x, y) {
#'            panel.xyplot(x, y, type = c("p", "r"), grid = TRUE)
#'            panel.abline(a = 0, b = 1, lty = 2)
#'        })
NULL

#' @name soyaBeans
#' @title Avaliação da Cultura de Soja
#' @description Experimento conduzido na Universidade Federal da Grande
#'     Dourados (UFGD) em casa de vegetação com a cultura de soja. Foram
#'     experimentados diferentes fatores de água do solo e de adubação
#'     potássica, em cada parcela que continha 2 plantas. Para controle
#'     de variação local, as parcelas foram arranjadas em blocos.
#' @format Um \code{data.frame} com 75 observações e 10 variáveis.
#'     \describe{
#'
#' \item{\code{potassio}}{Número inteiro que representa a adubação
#'     potássica experimentada em cada parcela. Foram 5 nivéis de
#'     adubação considerados.}
#'
#' \item{\code{agua}}{Numérico que representa o nível de água do solo
#'     que as parcelas foram cultivadas. Foram 3 níveis experimentados
#'     pouca água, água em quantidade ideal, água em abundância.}
#'
#' \item{\code{bloco}}{Fator com 5 níveis que representam os blocos
#'     utilizados para controle de variação local.}
#'
#' \item{\code{rend}}{Rendimento de grãos.}
#'
#' \item{\code{peso}}{Peso de grãos.}
#'
#' \item{\code{kgrao}}{Conteúdo de potássio no grão.}
#'
#' \item{\code{pgrao}}{Conteúdo de fósforo no grão.}
#'
#' \item{\code{ts}}{Total de sementes por planta.}
#'
#' \item{\code{nvi}}{Número de vagens inviáveis.}
#'
#' \item{\code{nv}}{Número de vagens total.}
#'
#' }
#'
#' @details As variáveis \code{ts} e \code{nv} possuem alta correlação
#'     (0.978), pois em média são esperados 3 sementes por cada
#'     vagem. Então \code{ts} pode ser tratada como função de \code{nv}.
#' @keywords equidispersão superdispersão efeito-aleatório
#' @examples
#' data(soyaBeans)
#'
#' library(lattice)
#'
#' splom(soyaBeans[, -c(1:3)],
#'       type = c("p", "smooth"),
#'       grid = TRUE)
#'
#' # Para variável número de vagens em função dos fatores experimentais
#' xyplot(nv ~ potassio | factor(agua),
#'        data = soyaBeans,
#'        type = c("p", "spline"),
#'        grid = TRUE,
#'        as.table = TRUE,
#'        layout = c(NA, 1))
#'
#' # Para variável número de vagens inviáveis em função dos fatores
#' # experimentais
#' xyplot(nvi ~ potassio | factor(agua),
#'        data = soyaBeans,
#'        type = c("p", "spline"),
#'        grid = TRUE,
#'        as.table = TRUE,
#'        layout = c(NA, 1))
NULL


#' @name whiteFly
#' @title Ninfas de Mosca-Branca em Lavoura de Soja
#' @description Experimento conduzido em casa de vegetação sob o
#'     delineamento de blocos casualizados. No experimento foram
#'     avaliadas plantas de diferentes cultivares de soja contabilizando
#'     o número de ninfas de mosca-branca nos folíolos dos terços
#'     superior, médio e inferior das plantas. As avaliações ocorreram
#'     em 6 datas dentre os 38 dias do estudo.
#' @format Um \code{data.frame} com 240 observações e 8 variáveis.
#'     \describe{
#'
#' \item{\code{data}}{Data em que foram avaliadas as plantas de soja.}
#'
#' \item{\code{dias}}{Inteiro que indica o número de dias após o
#'     experimento no ato da avaliação.}
#'
#' \item{\code{cult}}{Fator com a identificação da cultivar de
#'     soja. Foram 10 cultivares avaliadas neste experimento.}
#'
#' \item{\code{bloco}}{Fator com 4 níveis que representam os blocos
#'     utilizados para controle de variação local.}
#'
#' \item{\code{nsup}}{Número de ninfas de mosca-branca nos folíolos do
#'     terço superior.}
#'
#' \item{\code{nmed}}{Número de ninfas de mosca-branca nos folíolos do
#'     terço médio.}
#'
#' \item{\code{ninf}}{Número de ninfas de mosca-branca nos folíolos do
#'     terço inferior.}
#'
#' \item{\code{ntot}}{Número de ninfas de mosca-branca considerando
#'     todos os folíolos (soma de \code{nsup}, \code{nmed},
#'     \code{ntot}).}
#'
#' }
#' @keywords superdispersão
#' @references Suekane, R., Degrande, P. E., de Lima Junior, I. S., de
#'     Queiroz, M. V. B. M., & Rigoni, E. R. (2013). Danos da
#'     Mosca-Branca Bemisia Tabaci e distribuição vertical das ninfas em
#'     cultivares de soja em casa de vegetação. Arquivos do Instituto
#'     Biológico, 80(2), 151-158.
#' @examples
#' data(whiteFly)
#'
#' library(lattice)
#'
#' xyplot(ntot ~ dias | cult,
#'        data = whiteFly,
#'        type = c("p", "spline"),
#'        grid = TRUE,
#'        as.table = TRUE,
#'        layout = c(NA, 2))
#'
#' # Somente as cultivares que contém BRS na identificação
#' da <- droplevels(subset(whiteFly, grepl("BRS", x = cult)))
#'
#' xyplot(ntot ~ dias | cult,
#'        data = da,
#'        type = c("p", "spline"),
#'        grid = TRUE,
#'        as.table = TRUE,
#'        layout = c(NA, 2))
NULL

#' @name eggs
#' @title Produção de Ovos por Galinha
#' @description Dados provenientes de um experimento com 300 galinhas
#'     poedeiras da genética \emph{Isa Brown}, submetidas à regimes de
#'     iluminação contínua de 17 horas por dia (12h natural e 5h
#'     artificial) em regime 17L:7E (luz: escuridão) com diferentes
#'     fontes para luz artificial durante 70 dias. O objetivo do estudo
#'     foi avaliar o impacto de diferentes cores de luz na produção de
#'     ovos. Para esta avaliação o experimento foi conduzido sob o
#'     delineamento de quadrado latino com cinco tratamentos (cinco
#'     cores), divididos em cinco períodos e cinco parcelas com seis
#'     repetições de dez aves em cada parcela. As variáveis de interesse
#'     foram número de ovos produzidos e a massa dos ovos.
#' @format Um \code{data.frame} com 2100 observações e 7 variáveis.
#'     \describe{
#'
#' \item{\code{periodo}}{Fator com cinco níveis que indica o período, de
#'     14 dias, do qual a observação provém.}
#'
#' \item{\code{box}}{Fator com cinco níveis que indica o box, com seis
#'     gaiolas de 10 galinhas, do qual a observação provém. }
#'
#' \item{\code{luz}}{Fator com cinco níveis (amarelo, azul, branco,
#'     verde e vermelho) que indica a cor da luz artificial sob
#'     incidência.}
#'
#' \item{\code{gaiola}}{Fator com seis níveis que indica a gaiola,
#'     com 10 galinhas, da qual a observação provém.}
#'
#' \item{\code{dia}}{Inteiro que indica o dia de observação no
#'     respectivo período.}
#'
#' \item{\code{ovos}}{Contagem do número de ovos produzidos por 10
#'     galinhas na respectiva parcela experimental (gaiola).}
#'
#' \item{\code{massa}}{Massa dos ovos produzidos na respectiva parcela
#'     experimental (gaiola), mensurada em gramas.}
#'
#' }
#' @keywords efeito-aleatório
#' @references Borille, Rodrigo, Garcia, Rodrigo G., Nääs, Irenilza A.,
#'     Caldara, Fabiana R., & Santana, Mayara R.. (2015). Monochromatic
#'     light-emitting diode (LED) source in layers hens during the
#'     second production cycle. Revista Brasileira de Engenharia
#'     Agrícola e Ambiental, 19(9),
#'     877-881. https://dx.doi.org/10.1590/1807-1929/agriambi.v19n9p877-881
#'
#' Borille, Rodrigo (2013). Led de diferentes cores como alternativa
#'     sustentável para iluminação de poedeiras comerciais. Dissertação
#'     de Mestrado, Universidade Federal da Grande Dourados, Grande
#'     Dourados, MS, Brasil.
#'     Dissertação.
#' @examples
#'
#' data(eggs)
#'
#' # 6 repetições de cada parcela experimental (combinação de período, luz
#' # e box) totalizando 84 observações em cada parcela. Experimento
#' # balanceado
#' ftable(xtabs(~periodo + luz + box, data = eggs))
#'
#' library(latticeExtra)
#'
#' # Contagem de ovos
#' useOuterStrips(
#'     xyplot(ovos ~ dia | luz + periodo,
#'            data = eggs,
#'            jitter.x = TRUE,
#'            groups = gaiola,
#'            type = c("p", "smooth")))
#'
#' # Massa de ovos
#' useOuterStrips(
#'     xyplot(massa ~ dia | luz + periodo,
#'            data = eggs,
#'            jitter.x = TRUE,
#'            groups = gaiola,
#'            type = c("p", "smooth")))
#'
#' # Relação das variáveis de interesse
#' xyplot(massa ~ ovos,
#'        data = eggs,
#'        type = c("p", "g", "smooth"))
#'
#' # Agrupando as contagens a cada 14 dias
#' da <- aggregate(ovos ~ periodo + box + luz + gaiola,
#'                 data = eggs, FUN = sum)
#'
#' xyplot(ovos ~ periodo,
#'        groups = luz,
#'        data = da,
#'        jitter.x = TRUE,
#'        type = c("p", "g", "smooth"))
NULL

#' @name cottonBolls
#' @title Capulhos de Algodão em Função da Exposição à Mosca Branca
#' @description Experimento conduzido na Universidade Federal da Grande
#'     Dourados (UFGD) em 2007, cujo objetivo foi avaliar os impactos da
#'     exposição de plantas à alta infestação de Mosca-Branca
#'     \emph{Bemisia tabaci} em componentes de produção do algodão. No
#'     experimento, plantas de algodão foram expostas à alta infestação
#'     da praga por períodos diferentes e ao final do experimento
#'     avaliou-se o número de capulhos produzidos, o número de
#'     estruturas reprodutivas, o número de nós, a altura da planta e o
#'     peso dos capulhos por vaso. A condução do estudo deu-se via
#'     delineamento interamente casualizado com 5 vasos, contendo duas
#'     plantas, para cada período de exposição.
#' @format Um \code{data.frame} com 60 observações e 8 variáveis.
#'     \describe{
#'
#' \item{\code{dexp}}{Inteiro com 6 valores que representa os dias de
#'     exposição à alta infestação de Mosca-Branca.}
#'
#' \item{\code{vaso}}{Fator que indica o vaso no qual foram mensurados
#'     os componentes de produção do algodão.}
#'
#' \item{\code{planta}}{Fator que indica a planta na qual foram
#'     mensurados os componentes de produção do algodão.}
#'
#' \item{\code{alt}}{Altura da planta, mensurada em centímetros.}
#'
#' \item{\code{pesocap}}{Peso dos capulhos de algodão, mesurado para
#'     cada vaso (que contém duas plantas). No \code{data.frame} somente
#'     a primeira planta do vaso contém a observação do peso.}
#'
#' \item{\code{nerep}}{Contagem do número de estruturas reprodutivas da
#'     planta.}
#'
#' \item{\code{ncapu}}{Contagem do número de capulhos produzidos.}
#'
#' \item{\code{nnos}}{Contagem do número de nós da planta.}
#'
#' }
#' @keywords efeito-aleatório subdispersão
#' @examples
#'
#' data(cottonBolls)
#'
#' library(lattice)
#'
#' # Número de estruturas reprodutivas
#' xyplot(nerep ~ dexp | vaso,
#'        data = cottonBolls,
#'        jitter.x = TRUE,
#'        type = c("p", "g", "smooth"))
#'
#' # Número de capulhos produzidos
#' xyplot(ncapu ~ dexp | vaso,
#'        data = cottonBolls,
#'        jitter.x = TRUE,
#'        type = c("p", "g", "smooth"))
#'
#' # Número de nós
#' xyplot(nnos ~ dexp | vaso,
#'        data = cottonBolls,
#'        jitter.x = TRUE,
#'        type = c("p", "g", "smooth"))
#'
#' # Altura das plantas
#' xyplot(alt ~ dexp | vaso,
#'        data = cottonBolls,
#'        jitter.x = TRUE,
#'        type = c("p", "g", "smooth"))
#'
#' # Peso dos capulhos nos vasos
#' da <- cottonBolls[complete.cases(cottonBolls), ]
#' xyplot(pesocap ~ dexp,
#'        data = cottonBolls,
#'        jitter.x = TRUE,
#'        type = c("p", "g", "smooth"))
NULL
