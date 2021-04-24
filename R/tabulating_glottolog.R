#' Simple metadata table
#' @return A dataframe with columns
#'   glottocode, isocode, name, 
#'   name_in_tree, position,
#'   tree, tree_name
glottolog_metadata = function() {
  glottolog_phylo_geo_v4.3 %>%
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
#' @return A dataframe with columns
#'   tree, tree_name
glottolog_families = function() {
  glottolog_phylo_geo_v4.3 %>%
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
    as.data.frame()
}


#' Extract parts from a glottolog tree label
#' @param label A string
#' @return A string
extract_name = function(labels) {
  regex <- "^[^\\[]+"
  str_extract(labels, regex)
}


#' Extract parts from a glottolog tree label
#' @param label A string
#' @return A string
extract_glottocode = function(labels) {
  regex <- "(?<=\\[)[a-z]{4}[0-9]{4}(?=\\])"
  str_extract(labels, regex)
}