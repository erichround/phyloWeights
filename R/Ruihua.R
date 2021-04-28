# Ruihua's tree

#' Yin (2020) tree
yin_2020_tree = function() {

  macro <- 
    list(c("South America", "North America"),
         c("Africa", "Eurasia"),
         "Papunesia",
         "Australia")

  y <-
    glottolog_supertree(macro_groups = macro) %>%
    keep_tips(yin_2021_langs$glottocode) %>%
    collapse_node("mada1298") %>%
    set_branch_lengths_exp() %>%
    set_first_branch_lengths(0.025) %>%
    label_by_glottocode(g)
    
  # Return the two trees as a multiPhylo object
  y
}


#' Yin (2020) tree
# yin_2020_tree_old = function() {
#   
#   # altered macroareas
#   phylo_geo_modified <- 
#     glottolog_phylo_geo() %>%
#     mutate(
#       main_macroarea = 
#         main_macroarea %>%
#         str_replace(".*America", "Americas") %>%
#         str_replace("Africa", "Eurasia")
#     )
#   
#   glottocode_tree <-
#     glottolog_supertree(
#       glottocodes = yin_2021_langs$glottocode,
#       phylo_geo = phylo_geo_modified
#     ) %>%
#     collapse_node("mada1298") %>%
#     set_branch_lengths_exp(ultra = FALSE) %>%
#     set_first_branch_lengths(0.025)
#   
#   name_tree <-
#     glottocode_tree %>%
#     label_by_name(
#       language_table = yin_2021_langs
#     )
#   
#   # Return the two trees as a multiPhylo object
#   c(glottocode_tree, name_tree)
# }


