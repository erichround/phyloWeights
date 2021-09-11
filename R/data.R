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