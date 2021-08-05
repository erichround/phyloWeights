#' Current glottolog version
#'
#' @return A character string. The
#'   newest version of glottolog for
#'   which data is included in this
#'   package.
get_glottolog_version = function() {
  "4.4"
}


#' Simple metadata table
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
#' @return A dataframe with columns
#'   glottocode, isocode, name, 
#'   name_in_tree, position,
#'   tree, tree_name
get_glottolog_languages = function(
  glottolog_version = get_glottolog_version()
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


#' Simple metadata table
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
#' @return A dataframe with columns
#'   tree, tree_name
get_glottolog_families = function(
  glottolog_version = get_glottolog_version()
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


#' Glottolog geographical data by version
#' 
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
#' @return A multiPhylo object, the
#'   glottolog trees of the corresponding
#'   version.
get_glottolog_phylo_geo = function(
  glottolog_version = get_glottolog_version()
) {
  
  if (is.numeric(glottolog_version[1])) {
    glottolog_version <- as.character(glottolog_version)
  }
  
  if (glottolog_version == "4.3") {
    phylo_geo <- glottolog_phylo_geo_v4.3
  } else if (glottolog_version == "4.4") {
    phylo_geo <- glottolog_phylo_geo_v4.4
  } else {
    stop("Available values for glottolog_version are '4.3' and '4.4'.")
  }
  
  phylo_geo
}


#' Get the tree number of one or more glottolog families
#' 
#' @param family A character string or
#'   vector of strings.
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
#' @return An integer or vector of integers.
#'   The tree number(s).
which_tree = function(
  family = NULL,
  glottolog_version = get_glottolog_version()
) {
  if (is.null(family)) {
    stop("family needs to be supplied.")
  }
  f <- get_glottolog_families(glottolog_version)
  f$tree[match(family, f$tree_name)]
}