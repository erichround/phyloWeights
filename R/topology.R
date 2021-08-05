# Manipulating tree topoplogy


#' Glottolog trees by version
#' 
#' @param family A character string.
#'   Which families' trees to return. If left
#'   unspecified, all trees are returned.
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
#' @return A multiPhylo object, the
#'   glottolog trees of the corresponding
#'   version; unless just one tree was
#'   requested, in which case the return
#'   is a phylo object of just one tree.
get_glottolog_trees = function(
  family = NULL,
  glottolog_version = get_glottolog_version()
) {
  
  if (is.numeric(glottolog_version[1])) {
    glottolog_version <- as.character(glottolog_version)
  }
  
  if (glottolog_version == "4.3") {
    phy <- glottolog_trees_v4.3
  } else if (glottolog_version == "4.4") {
    phy <- glottolog_trees_v4.4
  } else {
    stop("Available values for glottolog_version are '4.3' and '4.4'.")
  }
  
  if (is.null(family)) { return(phy) }
  
  # Select the requested subset of the families
  phy <- phy[which_tree(family)]
  
  if (length(phy) == 1) { 
    phy[[1]]
  } else {
    phy
  }
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


#' Create a glottolog super-tree
#' @param macro_groups A list
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
assemble_supertree = function(
  macro_groups = 
    as.list(unique(get_glottolog_families(glottolog_version)$main_macroarea)),
  glottolog_version = get_glottolog_version()
) {
  
  phy <- get_glottolog_trees(glottolog_version = glottolog_version)
  
  # Get the predominant macro_area for each tree
  main_macro <- 
    get_glottolog_families(glottolog_version)$main_macroarea
  
  # Group family trees by macroarea
  if (is.null(macro_groups)) {
    macro_phys <- phy
  } else {
    macro_phys <-
      lapply(macro_groups, function(m) {
        tree_set <- which(main_macro %in% m)
        macro_tree <- bind_as_rake(phy[tree_set])
        macro_label <- str_c(.name_to_label(m), collapse = "-")
        macro_tree$node.label[1] <- macro_label
        macro_tree
      }) %>%
      do.call(c, .)
  }
  
  # Group globally
  super_phy <- bind_as_rake(macro_phys)
  super_phy$node.label[1] <- "World"
  super_phy$edge.length <- rep(1, Nedge(super_phy))
  
  super_phy
}


#' Keep a set of nodes & tips, as tip set
#' @param phy A phylo object
#' @param labels A vector of strings,
#'   which are tip and node labels
select_tips = function(phy, label) {
  
  tip_labels <- phy$tip.label
  node_labels <- phy$node.label
  
  # Clone all named nodes as self-sister tips
  phy <- clone_node(phy, label)
  
  # Keep only the named tips, including the
  # nodes newly cloned as tips
  drop_tips <- setdiff(tip_labels, label)
  phy <- drop.tip(phy, drop_tips, collapse.singles = FALSE)
  
  # Duplicate tips that appear twice or more 
  # in the labels parameter, including the
  # nodes newly cloned as tips
  dup_tips <- label[duplicated(label)]
  phy <- clone_tip(phy, dup_tips)
  
  phy
}


#' Clone named nodes to self-daughter tips
#' 
#' A node is cloned once, even if listed
#' more than once in labels
#' 
#' @param phy A phlyo object
#' @param label A vector of strings,
#'   the labels of the nodes to clone
#' @return phy
clone_node = function(phy, label) {
  
  label <- 
    unique(label[label %in% phy$node.label])
  n_new <- length(label)
  
  if (n_new == 0) { return(phy) }
  
  n_edge <- Nedge(phy)
  n_tip <- Ntip(phy)
  n_node <- Nnode(phy)
  edges <- phy$edge
  new_edge <- matrix(NA, ncol = 2, nrow = n_edge + n_new)
  new_tip.label <- rep("", n_tip + n_new)
  new_edge.length <- rep(NA, n_edge + n_new)
  insert_e <- rep(NA, n_new)
  insert_t <- rep(NA, n_new)
  offset_e <- rep(0, n_edge)
  offset_t <- rep(0, n_tip)
  nodes <- rep(NA, n_new)
  
  # Work out how many positions to offset the current
  # edges and tips
  for (i in 1:n_new) {
    
    nodes[i] <- which(phy$node.label == label[i])
    u <- nodes[i] + Ntip(phy)
    
    # Find the first edge from u. This is the row, before 
    # which to insert a new row in phy$edge etc.
    insert_e[i] <- min(which(edges[,1] == u)) - 1
    # Update the offsets
    inc <- c(rep(0, insert_e[i]), rep(1, n_edge - insert_e[i]))
    offset_e <- offset_e + inc
    
    # Find previous tip in $edge. This is the tip number,
    # before which to insert a new tip in $tip.label etc.
    earlier_v <- edges[1:insert_e[i], 2]
    earlier_tip <- earlier_v[earlier_v <= n_tip]
    insert_t[i] <- max(c(0, earlier_tip))
    # Update the offsets
    inc <- c(rep(0, insert_t[i]), rep(1, n_tip - insert_t[i]))
    offset_t <- offset_t + inc
  }
  
  # What the old values in phy$edge map to:
  remap <- c((1:n_tip) + offset_t, (1:n_node) + n_tip + n_new)
  # Remap the old values to the right values and places:
  new_edge[(1:n_edge) + offset_e, ] <- remap[edges]
  new_edge.length[(1:n_edge) + offset_e] <- phy$edge.length
  new_tip.label[(1:n_tip) + offset_t] <- phy$tip.label
  
  # The insertion positions are now empty:
  empty_e <- sort(insert_e) + (1:n_new)
  empty_t <- sort(insert_t) + (1:n_new)
  # For empty edge rows, copy u from next row
  new_edge[empty_e, ] <- c(new_edge[empty_e + 1, 1], empty_t)
  new_edge.length[empty_e] <- 1
  # For empty tip labels, use node labels in order
  new_tip.label[empty_t] <- label[order(nodes)]
  
  new_tree <- list(
    edge = new_edge,
    Nnode = phy$Nnode,
    node.label = phy$node.label,
    tip.label = new_tip.label,
    edge.length = new_edge.length
  )
  class(new_tree) <- "phylo"
  new_tree
}


#' Clone named tips
#' @param phy A phlyo object
#' @param label A vector of strings,
#'   some of which might be tip labels
#'   of phy
#' @return phy
clone_tip = function(phy, label) {
  
  for (l in unique(label)) {
    n_dup <- sum(label == l)
    
    for(i in (1:n_dup)) {
      # clone a tip
      n_tips <- length(phy$tip.label)
      location <- 
        ifelse(i == 1,
               # The first time, bind to tip
               which(phy$tip.label == l),
               # Afterwards, bind to newly created node
               which(phy$node.label == "NA") + n_tips)
      phy <-
        phy %>%
        bind.tip(edge.length = 1,
                 # This is not the final label, but it
                 # keeps the tip labels distinct for now
                 tip.label = str_c(l, "-", i),
                 where = location)
    }
    
    # The original tip's branch length gets changed to 0,
    # so ensure it's also 1:
    dup_tips <- which(.remove_copy_suffix(phy$tip.label) == l)
    dup_tip_edges <- which(phy$edge[,2] %in% dup_tips)
    phy$edge.length[dup_tip_edges] <- 1
    
    # label the newly created node
    phy$node.label[phy$node.label == "NA"] <- l
  }
  
  phy$tip.label <- .add_copy_suffix(phy$tip.label)
  phy$node.label <- .add_copy_suffix(phy$node.label)
  
  phy
}


#' Collapse a node rootwards
#' 
#' Note this will also multichotomise
#' all nodes with a 0-length branch above 
#' them
#' @param phy A phylo object
#' @param label A string, the labels
#'   of the nodes to collapse
#' @return A phylo object
collapse_node = function(phy, label) {
  
  n_tips <- length(phy$tip.label)
  target_node <- match(label, phy$node.label) + n_tips
  target_branch <- match(target_node, phy$edge[,2])
  
  # Reduce the branch above the target node to 0
  phy$edge.length[target_branch] <- 0
  
  # Remove 0-length branches
  di2multi(phy)
}


.which_nonsister_copies = function(phy) {
  # Sister nodes in edge[,2] will share a
  # parent in edge[,1]. Check if this is true, 
  # and that the branch lengths are equal.
  df <-
    data.frame(
      u = phy$edge[,1],
      v = phy$edge[,2],
      len = phy$edge.length,
      v_label_orig = phy$tip.label[phy$edge[,2]],
      v_label = phy$tip.label[phy$edge[,2]] %>%
        .remove_copy_suffix(),
      stringsAsFactors = FALSE
      ) %>%
    filter(!is.na(v_label)) %>%
    group_by(v_label) %>%
    mutate(
      n_u = length(unique(u)),
      n_len = length(unique(len))
      ) %>%
    filter(n_u > 1 | n_len > 1)
  df$v_label_orig
}