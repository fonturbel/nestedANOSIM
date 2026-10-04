#' Animal visitors to two mistletoe species at Las Chinchillas National Reserve
#'
#' Camera-trap records of animal visitors to the mistletoes *Tristerix
#' aphyllus* and *T. verticillatus* in Las Chinchillas National Reserve,
#' northern Chile, over two sampling years. Each row is one camera
#' (one mistletoe individual) in one sampling year. Counts are the number of
#' records of each taxon.
#'
#' The focal species (*Mimus thenca* and *Sephanoides sephaniodes*) and
#' unidentified records were removed, as in the original analysis. Seven
#' cameras recorded none of the remaining taxa; they are kept so that the
#' automatic removal of empty samples in [nested_anosim()] can be seen.
#'
#' @format A data frame with 87 rows and 43 columns (3 metadata columns and 40 visitor taxa):
#' \describe{
#'   \item{Year_study}{Sampling year: `"First"` or `"Second"`.}
#'   \item{Mistletoe}{Mistletoe species: `"T_aphyllus"` or `"T_verticillatus"`.}
#'   \item{Camera}{Camera (mistletoe individual) code. The same cameras were
#'     used in both years.}
#'   \item{Anairetes_parulus}{Records of *Anairetes parulus*.}
#'   \item{butterfly}{Records of butterfly.}
#'   \item{Callipepla_californica}{Records of *Callipepla californica*.}
#'   \item{Callopistes_maculatus}{Records of *Callopistes maculatus*.}
#'   \item{Cicadidae_sp}{Records of *Cicadidae* sp.}
#'   \item{Colorhamphus_parvirostris}{Records of *Colorhamphus parvirostris*.}
#'   \item{cricket}{Records of cricket.}
#'   \item{Curaeus_curaeus}{Records of *Curaeus curaeus*.}
#'   \item{Diptera}{Records of Diptera.}
#'   \item{Diuca_diuca}{Records of *Diuca diuca*.}
#'   \item{Elaenia_albiceps}{Records of *Elaenia albiceps*.}
#'   \item{Falco_sparverius}{Records of *Falco sparverius*.}
#'   \item{Geositta_cunicularia}{Records of *Geositta cunicularia*.}
#'   \item{Geranoaetus_polyosoma}{Records of *Geranoaetus polyosoma*.}
#'   \item{grasshopper}{Records of grasshopper.}
#'   \item{Hemiptera}{Records of Hemiptera.}
#'   \item{Hymenoptera}{Records of Hymenoptera.}
#'   \item{Leistes_loyca}{Records of *Leistes loyca*.}
#'   \item{Leptasthenura_aegithaloides}{Records of *Leptasthenura aegithaloides*.}
#'   \item{Liolaemus_sp}{Records of *Liolaemus* sp.}
#'   \item{Lizard}{Records of Lizard.}
#'   \item{Lycalopex_griseus}{Records of *Lycalopex griseus*.}
#'   \item{Molothrus_bonariensis}{Records of *Molothrus bonariensis*.}
#'   \item{moth}{Records of moth.}
#'   \item{Muscisaxicola_maclovianus}{Records of *Muscisaxicola maclovianus*.}
#'   \item{Nothoprocta_perdicaria}{Records of *Nothoprocta perdicaria*.}
#'   \item{Phrygilus_alaudinus}{Records of *Phrygilus alaudinus*.}
#'   \item{Phrygilus_atriceps}{Records of *Phrygilus atriceps*.}
#'   \item{Phrygilus_gayi}{Records of *Phrygilus gayi*.}
#'   \item{Phyllotis_darwini}{Records of *Phyllotis darwini*.}
#'   \item{Porphyrospiza_alaudina}{Records of *Porphyrospiza alaudina*.}
#'   \item{Pseudasthenes_humicola}{Records of *Pseudasthenes humicola*.}
#'   \item{Pteroptochos_megapodius}{Records of *Pteroptochos megapodius*.}
#'   \item{Pyrisitia_lisa}{Records of *Pyrisitia lisa*.}
#'   \item{Pyrope_pyrope}{Records of *Pyrope pyrope*.}
#'   \item{Thylamys_elegans}{Records of *Thylamys elegans*.}
#'   \item{Troglodytes_aedon}{Records of *Troglodytes aedon*.}
#'   \item{Veniliornis_lignarius}{Records of *Veniliornis lignarius*.}
#'   \item{Zenaida_auriculata}{Records of *Zenaida auriculata*.}
#'   \item{Zonotrichia_capensis}{Records of *Zonotrichia capensis*.}
#' }
#'
#' @source Las Chinchillas National Reserve camera-trap project. See the
#'   original study at \doi{10.1016/j.jaridenv.2025.105518}. A semicolon-separated
#'   copy is available at
#'   `system.file("extdata", "mistletoe_visitors.csv", package = "nestedANOSIM")`.
#'   Please cite the paper if you use these data.
#'
#' @examples
#' data(mistletoe_visitors)
#' table(mistletoe_visitors$Year_study, mistletoe_visitors$Mistletoe)
#' comm <- mistletoe_visitors[, -(1:3)]
#' sort(colSums(comm), decreasing = TRUE)[1:5]
"mistletoe_visitors"
