# Creating a glottolog super-tree

#' Main function
#' @param phy A multiPhylo object holding
#'   trees in glottolog format
#' @param geo A data frame of glottolog
#'   langauge and dialect data
#' @param glottocodes A vector of strings
#' @param phy A multiphylo object in 
#'   glottolog format
#' @param has_macro_groups A logical
#' @param phylo_geo A dataframe whose
#'   columns include tree and main_macroarea
#' @param label_type A string
glottolog_super_tree = function(
  phy = glottolog_trees_v4.3,
  geo = glottolog_geography_v4.3,
  glottocodes = NULL,
  has_macro_groups = TRUE,
  phylo_geo = glottolog_phylo_geo(phy),
  label_type = c("glottocode", "name")
) {
  
  if (is.null(glottocodes)) {
    
    # If no glottocodes are specified,
    # use all available in the trees
    
    glottocodes <-
      glottolog_phylo_geo_v4.3 %>%
      filter(!is.na(vertex_type)) %>%
      .$glottocode
  } else {
    
    # If they are specified, provide some info
    
    vertex_types <- 
      glottolog_phylo_geo_v4.3 %>%
      filter(glottocode %in% glottocodes)$vertex_type
    is_dup_tip <- duplicated(glottocodes) & vertex_types == "tip"
    is_dup_node <- duplicated(glottocodes) & vertex_types == "node"
    dup_tip <- unique(glottocodes[is_dup_tip])
    dup_node <- unique(glottocodes[is_dup_node])
    n_dup_tip <- length(dup_tip)
    n_dup_node <- length(dup_node)
    dup_tip_msg <- 
      str_c(", of which ",
        str_c(head(dup_tip, 4), collapse = ", "),
        ifelse(n_dup_tip > 5, str_c("+", (n_dup_tip - 4), " more"),
               ifelse(n_dup_tip == 0, "none", "")),
        ifelse(n_dup_tip < 2, " occurs", " occur"),
        " more than once.")
    dup_node_msg <- 
      str_c(", of which ",
            str_c(head(dup_node, 4), collapse = ", "),
            ifelse(n_dup_tip > 5, str_c("+", (n_dup_node - 4), " more"),
                   ifelse(n_dup_node == 0, "none", "")),
            ifelse(n_dup_node < 2, " occurs", " occur"),
            " more than once.")   
      
    message(str_c(
      length(glottocodes), " glottocodes supplied, of which ",
      length(vertex_types), " found as tip or node.\n",
      sum(vertex_types == "tip"), " are glottolog tips", 
      dup_tip_msg, "\n",
      sum(vertex_types == "node"), " are glottolog nodes", 
      dup_node_msg
      ))
  }
  
  ## Repair nodes with one child
  # Add a dummy child to be deleted later
  
  for (i in 1:length(phy)) {
    phy[[i]] <- repair_bad_nodes(phy[[i]])
  }
  
  ## Establish macroarea grouping
  
  if (!has_macro_groups) {
    # When not grouping by macroarea,
    # put all families in one group, "World"
    phylo_geo$main_macroarea <- "World"
  }
  
  # Group family trees by macroarea
  macroareas <- unique(phylo_geo$main_macroarea)
  macroarea_multiPhylos <- list()
  for (m in macroareas) {
    tree_set <- 
      filter(phylo_geo, main_macroarea == m)$tree
    macroarea_multiPhylos[[m]] <- phy[tree_set]
  }
  
  ## Combine into high-level rake structure
  
  # First within each macroarea
  macro_rakes <- list()
  class(macro_rakes) <- "multiPhylo"
  for (i in 1:length(macroarea_multiPhylos)) {
    rake <- bind_as_rake(macroarea_multiPhylos[[i]])
    rake_name <- macroareas[[i]]
    rake$node.label[1] <- rake_name
    macro_rakes[[rake_name]] <- rake
  }
  # Then globally
  super_phy <- bind_as_rake(macro_rakes)
  
  ## Change node and tip labels to glottocodes
  
  is_rake_node <- super_phy$node.label %in% macroareas
  rake_node_labels <- super_phy$node.label[is_rake_node]
  super_phy <- label_by_glottocode(super_phy)
  super_phy$node.label[is_rake_node] <- rake_node_labels
  
  ## Filter by needed vertices, and duplicateÏ
  ## vertices to self-sisters as needed
  
  super_phy <- super_phy %>% keep_vertices(glottocodes)
  
  # Label root and set branch lengths to 1
  super_phy$node.label[1] <- "World"
  super_phy$edge.length <- 
    rep(1, length(super_phy$edge.length))
  
  super_phy
}
