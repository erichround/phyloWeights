# Ruihua's tree

#' Yin (2020) tree
yin_2020_tree = function() {

  # altered macroareas
  phylo_geo_modified <- 
    glottolog_phylo_geo() %>%
    mutate(
      main_macroarea = 
        main_macroarea %>%
        str_replace(".*America", "Americas") %>%
        str_replace("Africa", "Eurasia")
    )
  
  glottocode_tree <-
    glottolog_super_tree(
      glottocodes = yin_2021_langs$glottocode,
      phylo_geo = phylo_geo_modified
    ) %>%
    collapse_node("mada1298") %>%
    set_branch_lengths_exp(ultra = FALSE) %>%
    set_first_branch_lengths(0.025)

  name_tree <-
    glottocode_tree %>%
    label_by_name(
      language_table = yin_2021_langs
    )
  
  # Return the two trees as a multiPhylo object
  c(glottocode_tree, name_tree)
}