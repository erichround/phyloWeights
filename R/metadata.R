# Glottolog metadata


#' A simple glottolog language metadata table
#'
#' Returns a dataframe of glottolog metadata with columns: \code{glottocode},
#' \code{isocode}, \code{name}, \code{name_in_tree}, \code{position},
#' \code{tree} and \code{tree_name}.
#'
#' By default, returns metadata from the newest version of glottolog. Obtain
#' metadata for older versions by setting the parameter \code{glottolog_version}.
#'
#' @param glottolog_version A character string. Which glottolog version to use.
#'   Currently, available options are \code{'4.3'} and \code{'4.4'}.
get_glottolog_languages = function(
  glottolog_version = get_newest_version()
) {
  get_glottolog_phylo_geo(glottolog_version) %>%
    select(
      glottocode,
      isocode = isocodes,
      name,
      name_in_tree = vertex_name,
      position = vertex_type,
      tree = tree,
      tree_name = family_name
    )
}


#' A simple glottolog language family metadata table
#'
#' Returns a dataframe of glottolog metadata with columns: \code{tree},
#' \code{tree_name}, \code{n_tips}, \code{n_nodes} and \code{main_macroarea}. 
#'
#' By default, returns metadata from the newest version of glottolog. Obtain
#' metadata for older versions by setting the parameter
#' \code{glottolog_version}.
#'
#' @inheritParams get_glottolog_languages
get_glottolog_families = function(
  glottolog_version = get_newest_version()
) {
  main_macro <-
    get_glottolog_phylo_geo(glottolog_version) %>%
    filter(!is.na(tree)) %>%
    group_by(tree, macroarea) %>%
    summarise(n = n()) %>%
    group_by(tree) %>%
    arrange(-n) %>%
    slice(1) %>%
    select(tree, main_macroarea = macroarea)
  
  get_glottolog_phylo_geo(glottolog_version) %>%
    select(
      tree = tree,
      tree_name = family_name,
      vertex_type
    ) %>%
    filter(!is.na(tree)) %>%
    group_by(tree, tree_name) %>%
    summarise(
      n_tips = sum(vertex_type == "tip"),
      n_nodes = sum(vertex_type == "node")
      ) %>%
    arrange(tree) %>%
    left_join(main_macro, by = "tree") %>%
    as.data.frame()
}


#' An extended glottolog metadata table
#'
#' Returns a dataframe of glottolog geographical and phylogenetic metadata with
#' columns: \code{glottocode}, \code{name}, \code{vertex_name},
#' \code{family_name}, \code{isocodes}, \code{level}, \code{macroarea},
#' \code{latitude}, \code{longitude}, \code{vertex_type}, \code{vertex_label},
#' \code{tree} and \code{family_glottocode}.
#'
#' By default, returns metadata from the newest version of glottolog. Obtain
#' metadata for older versions by setting the parameter
#' \code{glottolog_version}.
#' 
#' @inheritParams get_glottolog_languages
get_glottolog_phylo_geo = function(
  glottolog_version = get_newest_version()
) {
  
  # Check glottolog_version
  error_msg <- .check_glottolog_version(glottolog_version)
  if (!is.na(error_msg)) { stop(error_msg) }
  
  # Choose appropriate dataset
  glottolog_version <- as.character(glottolog_version)
  if (glottolog_version == "4.3") {
    phylo_geo <- glottolog_phylo_geo_v4.3
  } else if (glottolog_version == "4.4") {
    phylo_geo <- glottolog_phylo_geo_v4.4
  }
  
  phylo_geo
}


#' The current glottolog version
#'
#' Returns the newest version of glottolog for which data is included in this
#' package, which is \code{'4.4'}.
#'
#' @return A character string.
get_newest_version = function() {
  "4.4"
}


#' Tree numbers of glottolog families
#' 
#' Returns the tree number of one or more glottolog families.
#' 
#' By default, returns metadata from the newest version of glottolog. Obtain
#' metadata for older versions by setting the parameter
#' \code{glottolog_version}.
#' 
#' @inheritParams get_glottolog_trees
#' @return A named vector of integers, giving the tree numbers and the
#'   family names as the vector names.
which_tree = function(
  family = NULL,
  glottolog_version = get_newest_version()
) {
  
  if (is.null(family)) {
    stop(str_c("A value for `family` needs to be supplied.\n",
               "You didn't supply one."))
  }
  
  # Check glottolog_version
  error_msg <- .check_glottolog_version(glottolog_version)
  if (!is.na(error_msg)) { stop(error_msg) }
  glottolog_version <- as.character(glottolog_version)
  
  # Check family
  check_result <- 
    .check_family(family, glottolog_version)
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  f <- get_glottolog_families(glottolog_version)
  tree_nums <- f$tree[match(family, f$tree_name)]
  names(tree_nums) <- family
  tree_nums
}


#' Vector of glottolog families
#' 
#' @inheritParams get_glottolog_languages
#' @return A vector.
#' @noRd
.get_family_vector = function(
  glottolog_version
) {
  get_glottolog_phylo_geo(glottolog_version) %>%
    .$family_name %>%
    unique()
}
