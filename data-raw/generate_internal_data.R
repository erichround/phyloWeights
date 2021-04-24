# Tabulated glottolog data

#' Tabulate gottolog tree vertices and geo data
#' To be saved as sysdata under the name
#' glottolog_phylo_geo_v4.3
#' @return A dataframe of gottolog
#'   nodes, tips, and their geo information
tabulate_language_geo = function(
  phy = glottolog_trees_v4.3,
  geo = languages_and_dialects_geo
) {
  
  # Compile a table of vertices in trees
  vert <- 
    lapply(
      1:length(phy),
      function(i) {
        bind_rows(
          data.frame(vertex_type = "node",
                     vertex_label = phy[[i]]$node.label,
                     stringsAsFactors = FALSE),
          data.frame(vertex_type = "tip",
                     vertex_label = phy[[i]]$tip.label,
                     stringsAsFactors = FALSE)
        ) %>% 
          mutate(tree = i) }
    ) %>% 
    bind_rows() %>%
    mutate(
      glottocode = extract_glottocode(vertex_label),
      vertex_name = extract_name(vertex_label)
    )
  
  # Add a column vertex_name, the equivalent of name
  # that we'd expect to see in a vertex label
  geo <- geo %>%
    mutate(vertex_name = name %>% 
             str_remove_all("[' ]") %>%
             str_replace_all("\\(", "{") %>%
             str_replace_all("\\)", "}"))
  
  # Combine the geo and vertex info
  full_join(
    geo, vert, 
    by = c("vertex_name", "glottocode")
  ) %>% 
    # Add family names
    left_join(
      tabulate_family_labels(phy), by = "tree"
    ) %>%
    select(glottocode, name, vertex_name,
           family_name, everything()) %>%
    arrange(glottocode)
}


#' Tabulate the labels of families' trees
#' #' To be saved as sysdata under the name
#' glottolog_family_labels_v4.3
#' @phy A multiPhylo object holding multiple
#'   trees in glottolog format
#' @return A dataframe with columns
#'   tree (an integer), family_name,
#'   family_glottocode
tabulate_family_labels = function(
  phy = glottolog_trees_v4.3
) {
  
  root_labels <- 
    lapply(phy, function(p) p$node.label[1]) %>%
    unlist()
  
  data.frame(
    tree = 1:length(phy),
    family_name = extract_name(root_labels),
    family_glottocode = extract_glottocode(root_labels),
    stringsAsFactors = FALSE
  )
}



#' Tabulate gottolog families and macroareas
#' #' To be saved as sysdata under the name
#' glottolog_family_geo_v4.3
#' @return A dataframe of gottolog
#'   nodes, tips, and their geo information
tabulate_family_geo = function(
  phy = glottolog_trees_v4.3,
  geo = glottolog_geography_v4.3
) {
  
  tabulate_language_geo(phy, geo) %>%
    filter(!is.na(macroarea), !is.na(tree)) %>% 
    group_by(tree, family_name, family_glottocode, macroarea) %>% 
    summarise(n = n()) %>% 
    group_by(tree, family_name, family_glottocode) %>%
    arrange(-n) %>%
    mutate(macroarea_n = str_c(macroarea, ":", n)) %>% 
    summarise(
      main_macroarea = macroarea[1],
      all_macroareas = str_c(macroarea_n, collapse = ", ")
    ) %>% 
    arrange(tree)
}