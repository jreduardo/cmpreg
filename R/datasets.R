#-----------------------------------------------------------------------
#' @title Artificial defoliation in cotton phenology
#'
#' @description Cotton production can be drastically reduced by attack
#'   of defoliating insects. Depending on the growth stage, the plants
#'   can recover from the caused damage and keep production not affected
#'   or can have the production reduced by low intensity defoliation. In
#'   order to study the recovery of cotton plants (\emph{Gossypium
#'   hirsutum}) in terms of production, Silva (2012) conducted a
#'   greenhouse experiment under a completely randomized design with
#'   five replicates. The experimental unity was a pot with two plants
#'   and it was recorded the number of cotton bolls at five artificial
#'   defoliation levels (0\%, 25\%, 50\%, 75\%, and 100\%) and five
#'   growth stages: vegetative, flower-bud, blossom, fig and cotton
#'   boll.
#'
#' @format A \code{\link[tibble]{tibble}} with 125 observations and 4
#'   colums:
#'
#' \itemize{
#' \item \code{stage}: A ordered factor with the (phenological) growth
#'   stages;
#' \item \code{defol}: Numerical with the defoliation levels (percent in
#'   leaf area removed with scissors);
#' \item \code{rept}: Indexes of repetition;
#' \item \code{bolls}: Number of bolls produced at harvest of cotton.
#'
#' }
#'
#' @usage data(cotton, package = "cmpreg")
#' @references Silva, A. M., Degrande, P. E., Suekane, R., Fernandes,
#'   M. G., Zeviani, W. M. (2012). Impacto de diferentes niveis de
#'   desfolha artificial nos estagios fenologicos do
#'   algodoeiro. \strong{Revista de Ciencias Agrarias}, 35(1), 163–172.
#'
"cotton"

#-----------------------------------------------------------------------
#' @title Soil moisture and potassium doses on soybean culture
#'
#' @description A study of potassium doses and soil moisture levels on
#'   soybean (\emph{Glicine Max}) production. The tropical soils are
#'   usually poor in potassium (K) and demand potassium fertilization
#'   when cultivated with soybean to obtain satisfactory yields. Soybean
#'   production is affected by long exposition to water deficit. As
#'   potassium is a nutrient involved in the water balance in plant, by
#'   hyphotesis, a good supply of potassium avoids to reduce
#'   production. To evaluate the effects of potassium doses and soil
#'   humidity levels on soybean production, Serafim (2012) conducted a
#'   \eqn{5\times 3} factorial experiment in a randomized complete
#'   block design with 5 replicates. Five different potassium doses (0,
#'   30, 60, 120 and 180 \eqn{\times} mg dm\eqn{^{-3}}) were applied to
#'   the soil and soil moisture levels were controlled at (37.5, 50, and
#'   62.5\%). The experiment was carried out in a greenhouse and the
#'   experimental units were pots with two plants in each.
#'
#' @format A \code{\link[tibble]{tibble}} with 73 observations and 6
#'   colums:
#'
#' \itemize{
#'
#' \item \code{K}: Integer value indicated the potassium fertilization
#'   dose (in mg dm\eqn{^{-3}});
#' \item \code{water}: Numerical value of amount of water in the soil
#'   (soil moisture in percent);
#' \item \code{block}: Fatcor indicating block;
#' \item \code{seeds}: Number of bean seeds at pot;
#' \item \code{vpods}: Number of viable pods at pot;
#' \item \code{tpods}: Total of pods at pot.
#'
#' }
#'
#' @usage data(soybean, package = "cmpreg")
#' @references Serafim, M. E., F. B. Ono, W. M. Zeviani, J. O. Novelino,
#'   and J. V. Silva (2012). Umidade do solo e doses de potassio na
#'   cultura da soja. \strong{Revista Ciencia Agronomica}, 43(2),
#'   222-227.
#'
"soybean"

#-----------------------------------------------------------------------
#' @title Toxicity of nitrofen in aquatic systems
#'
#' @description Nitrofen is no longer in commercial use in the United
#'   States, having been the first pesticide to be withdrawn due to
#'   tetragenic effects (Bailer, 1994). This data set comes from an
#'   experiment to measure the reproductive toxicity of the herbicide,
#'   nitrofen, on a species of zooplankton (\emph{Ceriodaphnia
#'   dubia}). Fifty animals were randomized into batches of ten and each
#'   batch was put in a solution with a measured concentration of
#'   nitrofen (0, 80, 160, 235 and 310 \eqn{\mu}g/litre. Subsequently,
#'   the number of live offspring was recorded.
#'
#' @format A \code{\link[tibble]{tibble}} with 50 observations and 2
#'   colums:
#'
#' \itemize{
#'
#' \item \code{dose}: Numeric value of the nitrofen concentration level
#'   (in \eqn{\mu}g/litre);
#' \item \code{noffs}: Number of live offspring.
#'
#' }
#'
#' @usage data(nitrofen, package = "cmpreg")
#' @references Bailer, A. and J. Oris (1994). Assessing toxicity of
#'   pollutant in aquatic systems. \strong{In Case Studies in Biometry},
#'   25-40.
#'
"nitrofen"

#-----------------------------------------------------------------------
#' @title Annona mucosa for control of Sitophilus zeamaus
#'
#' @description New control methods are necessary for stored grain pest
#'   management programs due to both the widespread problems of
#'   insecticide-resistance populations and the increasing concerns of
#'   consumers regarding pesticide residues in food products. Ribeiro
#'   (2013) carried out an experiment to assess the bioactivity of
#'   extracts of \emph{Annona mucosa} (Annonaceae) for control
#'   \emph{Sitophilus zeamaus} (Coleoptera: Curculionidae), a major pest
#'   of stored maize in Brazil. Petri dishes containing 10g of corn were
#'   treated with extracts prepared with different parts of \emph{Annona
#'   mucosa} (seeds, leaves and branches) or just water (control) were
#'   completely randomized with 10 replicates. Then 20 animals adults
#'   were placed in each Petri dish and the numbers of emerged insects
#'   (progeny) after 60 days were recorded.
#'
#' @format A \code{\link[tibble]{tibble}} with 40 observations and 2
#'   colums:
#'
#' \itemize{
#'
#' \item \code{extract}: Factor indicating the extrated used in the
#'   solution;
#' \item \code{ninsect}: Number of emerged insects.
#'
#' }
#'
#' @usage data(sitophilus, package = "cmpreg")
#' @references Ribeiro, L. P., J. D. Vendramim, K. U. Bicalho,
#'   M. S. Andrade, J. B. Fernandes, R. A. Moral, and C. G. B. Demetrio
#'   (2013). Annona mucosa Jacq. (Annonaceae): A promising source of
#'   bioactive compounds against Sitophilus zeamais Mots. (Coleoptera:
#'   Curculionidae). \strong{Journal of Stored Products Research}, 55,
#'   6-14.
#'
"sitophilus"

#-----------------------------------------------------------------------
#' @title Alternative substrats for bromeliad production
#'
#' @description This dataset comes from a randomized experiment
#'   conducted in a greenhouse in four blocks design with objective of
#'   evaluate five different recipients of alternative substrates for
#'   bromeliads (Kanashiro, 2008). All treatments contained peat and
#'   perlite and differed in the third component: Pinus bark, Eucalyptus
#'   bark, Coxim, coconut fiber and Xaxim. The response variable was the
#'   number of leaves per experimental unit (pot with initially eight
#'   plants), which was registered at 4, 173, 229, 285, 341, and 435
#'   days after planting.
#'
#' @format A \code{\link[tibble]{tibble}} with 120 observations and 4
#'   colums:
#'
#' \itemize{
#'
#' \item \code{treat}: Factor indicating the alternative substrates;
#' \item \code{time}: Days after planting;
#' \item \code{block}: Factor indicating block;
#' \item \code{nleaves}: The median number of leaves.
#'
#' }
#'
#' @usage data(bromelia, package = "cmpreg")
#'
"bromelia"
