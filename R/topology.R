# Manipulating tree topoplogy


#' Glottolog trees by version
#'
#' Returns a multiPhylo object containing all, or a requested subset, of the
#' glottolog trees.
#'
#' By default, trees are returned from the most recent version of glottolog.
#' Alternatively, an older version of glottolog can be specified.
#'
#' @param family A character vector. Elements are names of glottolog families
#'   whose trees are to be returned. If \code{family} is left unspecified, all
#'   trees are returned.
#' @inheritParams get_glottolog_languages
#' @return A \code{phylo} object containing one glottolog tree, or a
#'   \code{multiPhylo} object containing multiple glottolog trees.
get_glottolog_trees = function(
  family,
  glottolog_version
) {
  
  # Check glottolog_version
  if (missing(glottolog_version)) {
    glottolog_version <- .get_newest_version()
  } else {
    error_msg <- .check_glottolog_version(glottolog_version)
    if (!is.na(error_msg)) { stop(error_msg) }
    glottolog_version <- as.character(glottolog_version)
  }
  
  # Choose appropriate dataset
  if (glottolog_version == "4.3") {
    phy <- glottolog_trees_v4.3
  } else if (glottolog_version == "4.4") {
    phy <- glottolog_trees_v4.4
  }
  
  # If family is missing, return the whole dataset
  if (missing(family)) { return(phy) }
  
  # Else, select the requested subset of the families
  
  # Check family
  check_result <- .check_family(family, glottolog_version)
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }

  phy <- phy[which_tree(family, glottolog_version)]
  
  if (length(phy) == 1) { 
    # If only one family, return as phylo object.
    phy[[1]]
  } else {
    # Otherwise, return a multiPhylo.
    phy
  }
}


#' Bind trees as a high-level rake
#'
#' Takes a multiPhylo object containing multiple trees and combines them into a
#' single tree with a rake structure at its root, below which each tree appears
#' on its own branch.
#'
#' @param phy A multiphylo object containing the trees to be combined.
#' @return A phylo object, a single tree.
bind_as_rake = function(phy) {
  
  if (class(phy) != "multiPhylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class multiPhylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  n_trees <- length(phy)
  if (n_trees == 1) {
    warning("`phy` contained only one tree. No changes were made to it.")
    return(phy[[1]])
  }
  
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
#'
#' Combining glottolog family trees into one large tree. Families can be
#' assembled directly below a rake structure at the root, or can be grouped, so
#' that the root first branches into groups, and the families then branch out
#' below the group nodes.
#'
#' Grouping is controlled by the \code{macro_groups} parameter. Groups can
#' comprise a single glottolog macroarea, or multiple macroareas. Current
#' macroareas are \code{Africa}, \code{Australia}, \code{Eurasia}, \code{North
#' America}, \code{Papunesia} and \code{South America}. Setting
#' \code{macro_groups} to \code{NULL} causes the tree to be assembled without
#' groups.
#'
#' @param macro_groups A list of character vectors, in whic each vector contains
#'   the names of one or more macroareas which define a group. Alternatively,
#'   setting \code{macro_groups} to \code{NULL} causes the tree to be assembled
#'   without groups.
#' @inheritParams get_glottolog_languages
assemble_supertree = function(
  macro_groups,
  glottolog_version
) {
  
  # Check glottolog_version
  if (missing(glottolog_version)) {
    glottolog_version <- .get_newest_version()
  } else {
    error_msg <- .check_glottolog_version(glottolog_version)
    if (!is.na(error_msg)) { stop(error_msg) }
    glottolog_version <- as.character(glottolog_version)
  }
  
  phy <- get_glottolog_trees(glottolog_version = glottolog_version)
  
  # Check macro_groups
  available_macros <-
    unique(get_glottolog_families(glottolog_version)$main_macroarea)
  if (missing(macro_groups)) { 
    # No argument given; use the available macro groups without
    # further grouping
    macro_groups <- list(available_macros) 
  }
  if (!is.null(macro_groups)){
    if (!is.list(macro_groups)) {
      if (is.character(macro_groups) & 
          all(macro_groups %in% available_macros)) {
        stop(str_c("`macro_groups` must be a list.\n",
                     "You have supplied character vector.\n",
                     "Did you mean `list(", 
                     str_c(str_c("'", macro_groups, "'"), 
                           collapse = ","),
                     ")`?"))
      }
      mc <- class(macro_groups)
      stop(str_c("`macro_groups` must be a list.\n",
                   "You have supplied an object of class ", mc, "."))
    }
  }

  # Get the predominant macro_area for each tree
  main_macro <- 
    get_glottolog_families(glottolog_version)$main_macroarea
  
  # Group family trees by macroarea
  if (is.null(macro_groups)) {
    # if phy = NULL, use no grouping at all
    macro_phys <- phy
  } else {
    # Do grouping
    macro_phys <-
      lapply(macro_groups, function(m) {
        tree_set <- which(main_macro %in% m)
        macro_tree <- 
          suppressWarnings(bind_as_rake(phy[tree_set]))
        macro_label <- str_c(.name_to_label(m), collapse = "-")
        macro_tree$node.label[1] <- macro_label
        macro_tree
      }) %>%
      do.call(c, .)
  }
  
  # Group globally
  super_phy <- 
    suppressWarnings(bind_as_rake(macro_phys))
  super_phy$node.label[1] <- "World"
  super_phy$edge.length <- rep(1, Nedge(super_phy))
  
  super_phy
}


#' Select tips
#'
#' From a tree, select which tips are to be kept.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip labels.
#' @return A phylo object containing the modified tree.
select_tip = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  # Keep only the named tips, including the
  # nodes newly cloned as tips
  drop_tips <- setdiff(phy$tip.label, label)
  phy <- drop.tip(phy, drop_tips, collapse.singles = FALSE)
  
  phy
}


#' Clone named nodes to self-daughter tips
#'
#' Clones internal nodes in a tree as self-daughter tips.
#'
#' The length of any new branch, between node n and its new clone, is set equal to the
#' longest of the original branches directly below node n.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing node labels.
#' @return A phylo object containing the modified tree.
clone_node = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  phy$node.label[phy$node.label == ""] <- "##NOLABEL##"
  phy$tip.label[phy$tip.label == ""] <- "##NOLABEL##"
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  # Check for duplicate nodes that match label
  l_nodes <- which(phy$node.label %in% label)
  dup_l_nodes <- duplicated(phy$node.label[l_nodes])
  dup_node_labels <- unique(phy$node.label[l_nodes[dup_l_nodes]])
  if (any(dup_l_nodes)) {
    stop(str_c("Cannot clone already-duplicated nodes: ",
               "`label` must not contain names of node labels ",
               "that occur more than once in `phy`.\n",
               "You supplied one or more labels that match ",
               "more than one node: ",
               str_c(head(dup_node_labels, 4), collapse = ","),
               ifelse(length(dup_node_labels) > 4, "..", ""), "."
    ))
  }
  
  label <- unique(label[label %in% phy$node.label])
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
  insert_e_length <- rep(NA, n_new)
  offset_e <- rep(0, n_edge)
  offset_t <- rep(0, n_tip)
  nodes <- rep(NA, n_new)
  
  # Work out how many positions to offset the current
  # edges and tips
  for (i in 1:n_new) {
    
    nodes[i] <- which(phy$node.label == label[i])
    u <- nodes[i] + Ntip(phy)
    
    u_edges <- which(edges[,1] == u)
    insert_e_length[i] <- max(phy$edge.length[u_edges])
    
    # Find the first edge from u. This is the row, before 
    # which to insert a new row in phy$edge etc.
    insert_e[i] <- min(u_edges) - 1
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
  new_edge.length[empty_e] <- insert_e_length
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
  
  new_tree$node.label[phy$node.label == "##NOLABEL##"] <- ""
  new_tree$tip.label[phy$tip.label == "##NOLABEL##"] <- ""
  
  new_tree
}


#' Clone named tips
#'
#' Clones tips as sisters of the original. Optionally, places the new clones and
#' the original in ther own subgroup, in which case the node for the new
#' subgroup is assigned the same label as the original tip.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing tip labels.
#' @param n A numeric vector. The number of clones to make.
#' @param subgroup A logical. Whether to create a subgroup containing the new
#'   clones and their original.
#' @return A phylo object containing the modified tree.
clone_tip = function(
  phy, 
  label, 
  n = 1, 
  subgroup = FALSE
) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  phy$node.label[phy$node.label == ""] <- "##NOLABEL##"
  phy$tip.label[phy$tip.label == ""] <- "##NOLABEL##"
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  # Check for duplicate tips that match label
  l_tips <- which(phy$tip.label %in% label)
  dup_l_tips <- duplicated(phy$tip.label[l_tips])
  dup_tip_labels <- unique(phy$tip.label[l_tips[dup_l_tips]])
  if (any(dup_l_tips)) {
    stop(str_c("Cannot clone already-duplicated tips: ",
               "`label` must not contain names of tip labels ",
               "that occur more than once in `phy`.\n",
               "You supplied one or more labels that match ",
               "more than one tip: ",
               str_c(head(dup_tip_labels, 4), collapse = ","),
               ifelse(length(dup_tip_labels) > 4, "..", ""), "."
               ))
  }
  
  # Check n
  if (!is.numeric(n)) {
    cn <- class(n)
    stop(str_c("`n` must be numeric.\n",
               "You supplied an object of class ", cn, "."))
  }
  n_n <- length(n)
  n_label <- length(label)
  if (n_n != 1 & n_n != n_label) {
    stop(str_c("`n` must be length 1 or the same length as `label`.\n",
               "`label` is length ", n_label, " but `n` is length ", n_n, "."))
  } else if (n_n == 1) {
    n <- rep(n, n_label)
  }
  
  n_label <- length(label)
  is_dupl <- duplicated(label)
  
  for (i in (1:n_label)[!is_dupl]) {
    l <- label[i]
    l_tip <- which(phy$tip.label == l)[1]
    orig_edge_length <- phy$edge.length[which(phy$edge[,2] == l_tip)]
    
    for(j in 1:n[i]) {
      
      n_tips <- length(phy$tip.label)
      l_tip <- which(phy$tip.label == l)[1]
      l_parent <- 
        ifelse(subgroup, 
               which(phy$node.label == "NA") + n_tips,
               phy$edge[phy$edge[,2] == l_tip, 1])
      
      # Add a tip at the destination
      phy <-
        phy %>%
        .bind_tip(
          edge.length = 1,
          # Temporary label avoids unexpected behaviour from bind.tip():
          tip.label = str_c("#TEMP#"),
          where = ifelse(subgroup & (j == 1), l_tip, l_parent)
          )
      phy$tip.label[phy$tip.label == "#TEMP#"] <- l
    }
    
    # Set total edge lengths the same as the original:
    n_tips <- length(phy$tip.label)
    is_target_daughter <- c(phy$tip.label == l, phy$node.label == "NA")
    target_edges <- which(phy$edge[,2] %in% which(is_target_daughter))
    if (subgroup) {
      phy$edge.length[target_edges] <- orig_edge_length / 2
    } else {
      phy$edge.length[target_edges] <- orig_edge_length
    }
    
    # label the newly created node
    phy$node.label[phy$node.label == "NA"] <- l
  }
  
  phy$node.label[phy$node.label == "##NOLABEL##"] <- ""
  phy$tip.label[phy$tip.label == "##NOLABEL##"] <- ""
  
  phy
}


#' Collapse a node rootwards
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param label A character vector containing node labels.
#' @return A phylo object containing the modified tree.
collapse_node = function(phy, label) {
  
  # Check phy
  check <- .check_phy(phy)
  if (!is.na(check$error_msg)) { stop(check$error_msg) }
  # Note, if needed, this will add missing $node.label, $tip.label, $edge.length
  phy <- check$phy
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
  n_tips <- length(phy$tip.label)
  target_node <- match(label, phy$node.label) + n_tips
  target_branch <- match(target_node, phy$edge[,2])
  
  # Reduce the branch above the target node to 0
  phy$edge.length[phy$edge.length == 0] <- 0.11123
  phy$edge.length[target_branch] <- 0
  # Remove 0-length branches
  phy <- di2multi(phy)
  phy$edge.length[phy$edge.length == 0.11123] <- 0
  
  phy
}


#' Safely bind tips
#'
#' \code{ape::bind.tip()} deletes tip and node label substrings enclosed in
#' square brackets. This is a safe version which doesn't.
#'
#' @noRd
.bind_tip = function(phy, tip.label, ...) {
  
  tip.label <- 
    tip.label %>%
    str_replace_all("\\[", "〔") %>%
    str_replace_all("\\]", "〕")
  
  phy$tip.label <- 
    phy$tip.label %>%
    str_replace_all("\\[", "〔") %>%
    str_replace_all("\\]", "〕")
  
  phy$node.label <- 
    phy$node.label %>%
    str_replace_all("\\[", "〔") %>%
    str_replace_all("\\]", "〕")
  
  phy <- bind.tip(phy, tip.label, ...)
  
  phy$tip.label <- 
    phy$tip.label %>%
    str_replace_all("〔", "\\[") %>%
    str_replace_all("〕", "\\]")
  
  phy$node.label <- 
    phy$node.label %>%
    str_replace_all("〔", "\\[") %>%
    str_replace_all("〕", "\\]")
  
  phy
}