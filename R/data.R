#' Trees from glottolog v4.3
#'
#' A multiPhylo object which provides a representation of 
#' the phylogenetic relationships of the languages in
#' glottolog.
#'
#' @format A multiPhylo object
#' @source \url{https://glottolog.org/meta/downloads}
"glottolog_trees_v4.3"


#' Trees from glottolog v4.4
#'
#' A multiPhylo object which provides a representation of the phylogenetic
#' relationships of the languages in glottolog.
#'
#' @format A multiPhylo object
#' @source \url{https://glottolog.org/meta/downloads}
"glottolog_trees_v4.4"


#' Geographical data from glottolog v4.3
#'
#' A dataset of geographical information about the languages in glottolog.
#'
#' @format A dataframe:
#' \describe{
#'   \item{glottocode}{glottocode of the lect}
#'   \item{name}{name of the lect}
#'   \item{isocodes}{the ISO-639-3 code of the lect}
#'   \item{level}{"language" or "dialect"}
#'   \item{macroarea}{glottolog's geographical macroarea}
#'   \item{latitude}{}
#'   \item{longitude}{}
#' }
#' @source \url{https://glottolog.org/meta/downloads}
"glottolog_geography_v4.3"


#' Geographical data from glottolog v4.4
#'
#' A dataset of geographical information about the languages in glottolog.
#'
#' @format A dataframe:
#' \describe{
#'   \item{glottocode}{glottocode of the lect}
#'   \item{name}{name of the lect}
#'   \item{isocodes}{the ISO-639-3 code of the lect}
#'   \item{level}{"language" or "dialect"}
#'   \item{macroarea}{glottolog's geographical macroarea}
#'   \item{latitude}{}
#'   \item{longitude}{}
#' }
#' @source \url{https://glottolog.org/meta/downloads}
"glottolog_geography_v4.4"


#' Sonority sequencing violation data in 496 languages, 
#' from Yin (2020)
#'
#' A dataframe containing language identifiers and data on the presence of SSP
#' violations in word-initial onsets and word-final codas, as used in:
#'
#' Yin, Ruihua. 2020. “Violations of the Sonority Sequencing Principle: How, and
#' How Often?” 2020 Conference of the Australian Linguistic Society, Brisbane.
#' https://als.asn.au/Conference/Program
#'
#' @format A dataframe:
#' \describe{
#'   \item{name}{name of the lect}
#'   \item{glottocode}{glottocode of the lect}
#'   \item{has_onset_violation}{presence (1) or absence (0) of word-initial SSP violation}
#'   \item{has_coda_violation}{presence (1) or absence (0) of word-final SSP violation}}
#' @source Yin (2020)
#'   
"yin_2020_data"