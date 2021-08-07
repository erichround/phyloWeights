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
#' @return A multiPhylo object containing the glottolog trees.
get_glottolog_trees = function(
  family,
  glottolog_version = get_newest_version()
) {
  
  # Check glottolog_version
  error_msg <- .check_glottolog_version(glottolog_version)
  if (!is.na(error_msg)) { stop(error_msg) }
  
  # Choose appropriate dataset
  glottolog_version <- as.character(glottolog_version)
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

  phy <- phy[which_tree(family)]
  
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
    warning("`phy` contained only one tree.")
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
#' Combing glottolog family trees into one large tree. Families can be assembled
#' directly below a rake structure at the root, or can be grouped, so that the
#' root first branches into groups, and the families then branch out below the
#' group nodes.
#'
#' Grouping is controlled by the \code{macro_groups} parameter. Groups can
#' comprise a single glottolog macroarea, or multiple macroareas. Current
#' macroareas are \code{Africa}, \code{Australia}, \code{Eurasia}, \code{North
#' America}, \code{Papunesia} and \code{South America}. Setting
#' \code{macro_groups} to \code{NULL} causes the tree to be assembled without
#' groups.
#'
#' @param macro_groups A list of character vectors, in whic each vector contains the
#'   names of one or more macroareas which define a group. Alternatively, setting
#'   \code{macro_groups} to \code{NULL} causes the tree to be assembled without
#'   groups.
#' @inheritParams get_glottolog_languages
assemble_supertree = function(
  macro_groups,
  glottolog_version = get_newest_version()
) {
  
  # Check glottolog_version
  error_msg <- .check_glottolog_version(glottolog_version)
  if (!is.na(error_msg)) { stop(error_msg) }
  glottolog_version <- as.character(glottolog_version)
  
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


#' Select tips from current tips and nodes
#'
#' From a tree, select existing tips to be kept, or to be cloned, and existing
#' internal nodes to be cloned as tips. To keep a tip, include its tip label
#' once in \code{label}. To clone a tip, resulting in \eqn{n} copies, includes
#' its tip label \ean{n} times in \code{label}. To clone an internal node as a
#' self-daughter tip, include its node label in \coded{label}.
#'
#' Cloned tips will appear under a new node of the same name. Duplicate tips and
#' nodes will be assigned labels with a suffix \code{-1}, \code{-2}, \code{-3},
#' ...
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param labels A character vector containing tip and node labels.
#' @return A phylo object containing the modified tree.
select_tips = function(phy, label) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "both")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }

  # Clone all named nodes as self-sister tips
  node_labels <- label[label %in% phy$node.label]
  if (length(node_labels) > 0) {
    phy <- clone_node(phy, node_labels)
  }
  
  # Keep only the named tips, including the
  # nodes newly cloned as tips
  drop_tips <- setdiff(phy$tip.label, label)
  phy <- drop.tip(phy, drop_tips, collapse.singles = FALSE)
  
  # Duplicate tips that appear twice or more 
  # in the labels parameter, including the
  # nodes newly cloned as tips
  dup_labels <- label[duplicated(label)]
  dup_tips <- dup_labels[dup_labels %in% phy$tip.label]
  if (length(dup_tips) > 0) {
    phy <- clone_tip(phy, dup_tips)
  }
  
  phy
}


#' Clone named nodes to self-daughter tips
#'
#' Clones internal nodes in a tree as self-daughter tips.
#'
#' A node is cloned only once, even if its label appears more than once in
#' \code{label}.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param labels A character vector containing node labels.
#' @return A phylo object containing the modified tree.
clone_node = function(phy, label) {
  
  # Check phy
  if (class(phy) != "phylo") {
    stop("phy must be of class phylo.")
  }
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "node")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
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
#'
#' #' Clones tips in a tree as self-sisters, under a new parent node.
#'
#' To clone a tip, resulting in \eqn{n} copies, include its tip label \ean{n-1}
#' times in \code{label}.
#'
#' Cloned tips will appear under a new node of the same name. Duplicate tips
#' will be assigned labels with a suffix \code{-1}, \code{-2}, \code{-3}, ...
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param labels A character vector containing tip labels.
#' @return A phylo object containing the modified tree.
#' @return phy
clone_tip = function(phy, label) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check labels
  check_result <- .check_labels(phy, label, type = "tip")
  if (!is.na(check_result$error_msg)) {
    stop(check_result$error_msg)
  } else if (!is.na(check_result$warning_msg)) {
    warning(check_result$warning_msg) 
  }
  
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
        .bind_tip(edge.length = 1,
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
#' Note this will also multichotomise all nodes with a 0-length branch above
#' them.
#'
#' @param phy A phylo object. The tree to manipulate.
#' @param labels A character vector containing node labels.
#' @return A phylo object containing the modified tree.
collapse_node = function(phy, label) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
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
  phy$edge.length[target_branch] <- 0
  
  # Remove 0-length branches
  di2multi(phy)
}


#' Find topologically distinct duplicate tips
#'
#' Identifies duplicate tips that aren't sisters on equal-length branches. Is
#' used by \code{phylo_average()} to check assumptions.
#'
#' @param phy A phylo object, the tree to examine.
#' @return A character vector of tip labels.
#' @noRd
.which_nonsister_copies = function(phy) {
  # Sister nodes in edge[,2] will share a parent in edge[,1]. Check if this is
  # true, and that the branch lengths are equal.
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
    # Keep only edges that end in tips
    filter(!is.na(v_label)) %>%
    # Group by tip label
    group_by(v_label) %>%
    mutate(
      # Count the number of different parents
      n_u = length(unique(u)),
      # Count the number of different lengths
      n_len = length(unique(len))
      ) %>%
    # Keep only those with >1 of either
    filter(n_u > 1 | n_len > 1)
  # Return their labels
  df$v_label_orig
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