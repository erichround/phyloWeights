#' Simple metadata table
#' @return A dataframe with columns
#'   glottocode, isocode, name, 
#'   name_in_tree, position,
#'   tree, tree_name
glottolog_languages = function() {
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
  main_macro <-
    glottolog_phylo_geo_v4.3 %>%
    filter(!is.na(tree)) %>%
    group_by(tree, macroarea) %>%
    summarise(n = n()) %>%
    group_by(tree) %>%
    arrange(-n) %>%
    slice(1) %>%
    select(tree, main_macroarea = macroarea)
  
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
    left_join(main_macro, by = "tree") %>%
    as.data.frame()
}