#' Trees from glottolog v4.3
#'
#' A multiPhylo object which provides a representation of 
#' the phylogenetic relationships of the languages in
#' glottolog.
#'
#' @format A multiPhylo object
#' @source \url{https://glottolog.org/meta/downloads}
"glottolog_trees_v4.3"


#' Geogrpahical data from glottolog v4.3
#'
#' A dataset of geographical information about the
#' languages in glottolog.
#'
#' @format A dataframe
#' \describe{
#'   \item{glottocode}{}
#'   \item{name}{name the lect}
#'   \item{isocodes}{the ISO-639-3 code of the lect}
#'   \item{level}{"language" or "dialect"}
#'   \item{macroarea}{glottolog's geographical macroarea}
#'   \item{latitude}{}
#'   \item{longitude}{}
#' }
#' @source \url{https://glottolog.org/meta/downloads}
"glottolog_geography_v4.3"