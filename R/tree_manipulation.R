# Tree manipulation


#' Shorten tree labels to the glottocode
#' 
#' Retain only the glottocode substring
#' within the tip and node labels
#' 
#' @param phy A phylo or multiPhylo object
#' @return A phylo object
label_by_glottocode = function(phy) {

  if (class(phy) == "phylo") {
    
    # Method for a single tree
    
    phy$tip.label <- extract_glottocode(phy$tip.label)
    phy$node.label <- extract_glottocode(phy$node.label) 
  
  } else if (class(phy) == "multiPhylo") {
    
    # Method for multiple trees: apply single
    # tree method to each
    
    for (i in 1:length(phy)) {
      phy[[i]] <- label_by_glottocode(phy[[i]])
    }
  }
  
  phy
}



#' Change glottocode tip and node labels to names
#' 
#' @param phy A phylo object
#' @param languages A dataframe of with columns
#'   glottocode and name
#' @return A phylo object
label_by_name = function(
  phy, 
  language_table = NULL
) {
  
  if (is.null(language_table)) {
    language_table <- 
      glottolog_phylo_geo_v4.3 %>%
      select(glottocode, name = vertex_name) %>%
      distinct()
  }
  
  # Add _1, _2 etc to repeated glottocodes
  language_table <-
    language_table %>%
    group_by(glottocode) %>%
    mutate(
      n = n(), 
      i = row_number(),
      suffix = ifelse(n > 1, str_c("_", i), "")
    ) %>% 
    ungroup() %>%
    mutate(glottocode = str_c(glottocode, suffix))
  
  # Tabulate glottocodes in the tree
  # and match names to them
  if (class(phy) == "multiPhylo") {
    
    labels_df <-
      lapply(phy, function(p) {
        bind_rows(
          data.frame(
            glottocode = extract_glottocode(p$tip.label),
            stringsAsFactors = FALSE),
          data.frame(
            glottocode = extract_glottocode(p$node.label),
            stringsAsFactors = FALSE)
        )
      }) %>%
      bind_rows() %>%
      left_join(language_table, by = "glottocode")
    
  } else {
    
    labels_df <-
      bind_rows(
        data.frame(
          glottocode = extract_glottocode(phy$tip.label),
          stringsAsFactors = FALSE),
        data.frame(
          glottocode = extract_glottocode(phy$node.label),
          stringsAsFactors = FALSE)
      ) %>%
      left_join(language_table, by = "glottocode")
  }
  
  n_missing <- sum(is.na(labels_df$name))
  
  # If a name is missing, give warning and retain 
  # the glottocode
  if (n_missing > 0) {
    warning(
      str_c(
        n_missing, " glottocode(s) in the tree ",
        "not recognised; retaining the glottocode ",
        "as the label in that case."
        )
    )
    labels_df <- labels_df %>%
      mutate(name = ifelse(is.na(name), 
                           glottocode, name))
  }
    
  if (class(phy) == "multiPhylo") {
    
    # Method for multiple trees: apply single
    # tree method to each
    
    phy <- label_by_glottocode(phy)
    
    for (i in 1:length(phy)) {
      phy_i <- phy[[i]]
      tip_matches <- match(phy_i$tip.label, labels_df$glottocode)
      node_matches <- match(phy_i$node.label, labels_df$glottocode)
      phy[[i]]$tip.label <- labels_df$name[tip_matches]
      phy[[i]]$node.label <- labels_df$name[node_matches]
    }
    
  } else if (class(phy) == "phylo") {
      
    # Method for a single tree
    
    phy <- label_by_glottocode(phy)
    tip_matches <- match(phy$tip.label, labels_df$glottocode)
    node_matches <- match(phy$node.label, labels_df$glottocode)
    phy$tip.label <- labels_df$name[tip_matches]
    phy$node.label <- labels_df$name[node_matches]
  
  }
  
  phy
}



#' Add a root to a tree
#' @param phy A phylo object
#' @reutrn A phylo object with
add_root = function(phy) { 
  if (length(phy) == 1) {
    phy$root.edge <- phy$edge[,1]
  } else {
    phy$root.edge <- setdiff(phy$edge[,1], phy$edge[,2])
  }
  phy 
}



#' Set lengths of first branches
#' @param phy A phlyo object
#' @return A phlyo object
set_first_branch_lengths = function(phy, branch_length = 1) {
  root_edge <- add_root(phy)$root.edge
  first_edges <- which(phy$edge[,1] == root_edge)
  phy$edge.length[first_edges] <- branch_length
  phy
}


#' Set branch length to 1
#' @param phy A phylo object
#' @param ultrametric A logical. Whether or not
#'   to lengthen branches immediately above
#'   tips to make the tree ultrametric
#' @return A phylo object
set_branch_lengths_1 = function(phy) {
  
  phy$edge.length <- rep(1, Nedge(phy))
  phy
  
}


#' Ultrametricise tree by stretching final edges
#' @param phy A phylo object
#' @return A phylo object
ultrametricise = function(phy) {
  
  if (any(is.na(phy$edge.length))) {
    warning("Converting NA branch lengths to 1")
    phy$edge.length[is.na(phy$edge.length)] <- 1
  }
  
  n_tips <- Ntip(phy)
  tip_depths <- node.depth.edgelength(phy)[1:n_tips]
  added_depth <- max(tip_depths) - tip_depths
  tip_edges <- match(1:n_tips, phy$edge[,2])
  phy$edge.length[tip_edges] <- 
    phy$edge.length[tip_edges] + added_depth
  
  phy
}


#' Exponentialise branch lengths
#' @param phy A phylo object
#' @return A phylo object
set_branch_lengths_exp = function(phy, ultrametric = FALSE) {
  
  nonroot <- phy$edge[,2]
  
  # Get node depths for non-roots
  phy$edge.length <- rep(1, Nedge(phy))
  depth <- node.depth.edgelength(phy)[nonroot]
  
  # Assign new branch length according to depth
  phy$edge.length <- (1/2) ^ depth
  phy
  
}


#' Exponentialise branch lengths
#' @param phy A phylo object
#' @return A phylo object
set_branch_lengths_exp_old = function(phy, ultrametric = FALSE) {

  n_tips <- Ntip(phy)
  root_node <- n_tips + 1
  
  # Depths of the nodes
  phy$edge.length <- rep(1, Nedge(phy))
  nonroot_depths <- node.depth.edgelength(phy)[-root_node]
  n_levels <- max(nonroot_depths)
  
  # Lengths of branches at each depth
  exp_lengths <- 1 / 2 ^ (1:n_levels)
  
  # Length for tip branches at each depth,
  # if ultrametric = TRUE
  ultra_lengths <- exp_lengths * 2 - exp_lengths[n_levels]
  
  # Assign lengths to edges
  phy$edge.length <-
    as.data.frame(phy$edge) %>%
    rename(tipward_node = V2) %>%
    mutate(original_order = row_number()) %>%
    arrange(tipward_node) %>%
    mutate(depth = nonroot_depths) %>%
    arrange(original_order) %>%
    mutate(
      is_tip_branch = (tipward_node < root_node),
      length = 
        ifelse(is_tip_branch & ultrametric,
               ultra_lengths[depth],
               exp_lengths[depth])
    ) %>%
    .$length
  
  phy
}



#' Bind trees with high-level rake
#' @param phy A multiphylo object
#' @param rake_length A numeric, the length
#'    of the high-level rake branches
#' @return A phylo object, a single tree
bind_as_rake_old = function(phy, rake_length = 1) {
  
  if (class(phy) != "multiPhylo") {
    stop("phy must be of class multiPhylo.")
  }
  
  n_trees <- length(phy)

  # The rake to which the trees will be added
  # Include two extra branches which will be dropped
  # at the end, to avoid errors from bind.tree()
  rake_size <- n_trees + 2
  drop_labels <- c("drop1", "drop2")
  rake_labels <- c(str_c("temp_", 1:n_trees), drop_labels)
  rake <- starTree(rake_labels, 
                   branch.lengths = rep(1, rake_size))
  
  # Order should be 2,1,3,4,5...
  tree_order <- 
    if (n_trees == 1) { 1 }
    else if (n_trees == 2) { c(2, 1) }
    else { c(2, 1, 3:n_trees) }
  
  # Bind trees to the rake, 
  for (i in tree_order) {
    target_position <- which(rake$tip.label == str_c("temp_", i))
    rake <- bind.tree(rake, add_root(phy[[i]]), target_position)
  }
  
  # Drop the extra rake branches
  drop_tips <- which(rake$tip.label %in% drop_labels)
  rake <- drop.tip(rake, drop_tips)
  
  # set rake edge lengths and return
  rake %>% set_first_branch_lengths(rake_length)
}


#' Bind trees with high-level rake
#' @param phy A multiphylo object
#' @return A phylo object, a single tree
bind_as_rake = function(phy) {
  
  if (class(phy) != "multiPhylo") {
    stop("phy must be of class multiPhylo.")
  }
  
  n_trees <- length(phy)
  n_tip_vec <- Ntip(phy)
  n_tips <- sum(n_tip_vec)
  root <- n_tips + 1
  tip_offset_vec <- cumsum(c(0, n_tip_vec[-n_trees]))
  
  n_node_vec <- Nnode(phy) - (n_tip_vec == 1) # Remove extraneous node above isolate
  n_node_cum <- cumsum(c(1, n_node_vec[-n_trees]))
  node_offset_vec <- n_node_cum  + n_tips - n_tip_vec
  
  edge <-
    lapply(1:n_trees, function(i) {
      # Adjust node numbers
      e <- eo <- phy[[i]]$edge
      nt <- n_tip_vec[i]
      e[eo <= nt] <- eo[eo <= nt] + tip_offset_vec[i]
      e[eo > nt] <- eo[eo > nt] + node_offset_vec[i]
      # Add edge connecting to new root
      edge_to_root <- c(root, root + n_node_cum[i])
      e <- rbind(edge_to_root, e, deparse.level = 0)
      # Remove extraneous node above isolates
      if (n_tip_vec[i] == 1) {
        e <- matrix(c(root, e[2,2]), nrow = 1)
      }
      e
    }) %>%
    do.call(rbind, .)
    
  edge_lengths <- 
    lapply(phy, function(x) {
      if (!("edge.length" %in% names(x))) {
        x$edge.length <- rep(1, Nedge(x))
      }
      e <- c(1, x$edge.length)
      # Remove extraneous node above isolates
      if (Ntip(x) == 1) { e <- 1 }
      e
    } ) %>% 
    unlist()
  node_labels <- 
    lapply(phy, function(x) {
      if (Ntip(x) == 1) { 
        # Remove extraneous node above isolates
        NULL 
        } else { x$node.label }
    }) %>% 
    unlist() %>% c("", .)
  tip_labels <- unlist(lapply(phy, function(x) x$tip.label))
  
  bound_tree <- list(
    edge = edge,
    Nnode = length(node_labels),
    node.label = node_labels,
    tip.label = tip_labels,
    edge.length = edge_lengths
  )
  class(bound_tree) <- "phylo"
  bound_tree
}


#' Collapse a node rootwards
#' 
#' Note this will also multichotomise
#' all nodes with a 0-length branch above 
#' them
#' @param phy A phylo object
#' @param node_label A string, the label
#'   of the node to collapse
#' @reutrn A phylo object
collapse_node = function(phy, node_label) {
  
  n_tips <- length(phy$tip.label)
  target_node <- which(phy$node.label == node_label) + n_tips
  target_branch <- which(phy$edge[,2] == target_node)
  
  # Reduce the branch above the target node to 0
  phy$edge.length[target_branch] <- 0
  
  # Remove 0-length branches
  di2multi(phy)
}



#' Clone named nodes to self-daughter tips
#' 
#' A node is cloned once, even if listed
#' more than once in labels
#' 
#' @param phy A phlyo object
#' @param labels A vector of strings,
#'   some of which might be node labels
#'   of phy
#' @return phy
clone_nodes_by_label = function(phy, labels) {
  
  node_labels <- labels[labels %in% phy$node.label]
  
  if (length(node_labels) == 0) { return(phy) }

  for (l in unique(node_labels)) {
    # clone one node 
    n_tips <- length(phy$tip.label)
    location <- which(phy$node.label == l) + n_tips
    phy <-
      phy %>%
      bind.tip(edge.length = 1, 
               tip.label = l,
               where = location)
  }
  
  phy
}



#' Clone named tips
#' @param phy A phlyo object
#' @param labels A vector of strings,
#'   some of which might be tip labels
#'   of phy
#' @return phy
clone_tips_by_label = function(phy, labels) {
  
  tip_labels <- labels[labels %in% phy$tip.label]
  dup_labels <- tip_labels[duplicated(tip_labels)]
  
  if (length(dup_labels) == 0) { return(phy) }
  
  # To get the new tips in the right order, we
  # handle tip_2 first, then (funnily enough) 1,3,4...
  for (l in unique(dup_labels)) {
    n_dup <- sum(dup_labels == l) + 1
    label_2 = str_c(l, "_2")
    phy$tip.label[which(phy$tip.label == l)] <- label_2
    
    for(i in (1:n_dup)[-2]) {
      # clone a tip
      clone_label <- str_c(l, "_", i) 
      n_tips <- length(phy$tip.label)
      location <- 
        ifelse(i == 1,
               # The first time, bind to tip
               which(phy$tip.label == label_2),
               # Afterwards, bind to newly created node
               which(phy$node.label == "NA") + n_tips)
      phy <-
        phy %>%
        bind.tip(edge.length = 1,
                 tip.label = clone_label,
                 where = location)
    }
    
    # label the newly created node
    phy$node.label[phy$node.label == "NA"] <- l
    
  }
  
  phy
}



#' Keep a set of nodes & tips, as tip set
#' @param phy A phylo object
#' @param labels A vector of strings,
#'   which are tip and node labels
keep_vertices = function(phy, labels) {
  
  tip_labels <- labels[labels %in% phy$tip.label]
  
  # Clone all named nodes as self-sister tips
  phy <- clone_nodes_by_label(phy, labels)
  
  # Keep only the named tips, include the
  # nodes newly cloned as tips
  phy <- keep.tip(phy, unique(labels))
  
  # Duplicate tips that appear twice in
  # the labels parameter
  phy <- clone_tips_by_label(phy, labels)
  
  phy
}



#' Add dummy tips below nodes with one child
#'
#' The glottolog tree contains nodes with just
#' one child, which are illegal structures, and
#' will be removed by methods such as ape::bind.tree. 
#' Create dummy tips below them to prevent this.
repair_bad_nodes = function(phy) {
  
  # Identify any bad nodes
  edge <- phy$edge
  bad_nodes <-
    as.data.frame(edge) %>%
    group_by(V1) %>%
    summarise(n = n()) %>%
    filter(n == 1) %>%
    .$V1
  n_bad_nodes <- length(bad_nodes)

  if (n_bad_nodes == 0) { return(phy) }

  # New information for dummy tip(s)
  old_n_tips <- length(phy$tip.label)
  new_tips <- old_n_tips + (1:n_bad_nodes)
  new_edges <- matrix(c(bad_nodes, new_tips), ncol = 2)
  new_edgelengths <- rep(1, n_bad_nodes)
  new_tiplabels <- rep("dummy_tip", n_bad_nodes)
  
  # Raise the node numbers to accommodate
  # the new tip(s)
  edge[edge > old_n_tips] <- edge[edge > old_n_tips] + n_bad_nodes
  bad_nodes <- bad_nodes + n_bad_nodes
  
  # Add edges
  for (i in 1:n_bad_nodes) {
    bad_node <- bad_nodes[i]
    break_row <- which(edge[,1] == bad_node)
    new_row <- c(bad_node, new_tips[i])
    edge <- 
      rbind(deparse.level = 0,
        edge[1:break_row,],
        new_row,
        edge[-(1:break_row),]
      )
  }
  phy$edge <- edge

  # Add other information
  phy$edge.length <- phy$edge.length %>% c(new_edgelengths)
  phy$tip.label <- phy$tip.label %>% c(new_tiplabels)
  
  phy
}