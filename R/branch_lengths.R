# Setting branch lengths

#' Set length of deepest branches
#'
#' Sets lengths of branches immediately below the root to the same,
#' user-specified length.
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @param branch_length A numeric stating the branch length.
#' @return A phlyo object, the manipulated tree.
set_deepest_branch_lengths = function(
  phy, 
  branch_length = 1
) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check branch_length
  if (!is.numeric(branch_length)) {
    cb <- class(branch_length)
    stop(str_c("`branch_length` must be numeric.\n",
               "You supplied an object of class ", cb, "."))
  }
  if (length(branch_length) != 1) {
    stop(str_c("`branch_length` must be length 1.\n",
               "You supplied a vector length ", length(branch_length), "."))
  }
  if (branch_length < 0) {
    stop(str_c("`branch_length` must be non-negative.\n",
               "You supplied the negative value ", branch_length, "."))
  }
  
  root <- Ntip(phy) + 1
  first_edges <- which(phy$edge[,1] == root)
  phy$edge.length[first_edges] <- branch_length
  phy
}


#' Set all branch length to 1
#' 
#' Sets all branch lengths in a tree to 1.
#' 
#' @param phy A phlyo object, the tree to manipulate.
#' @return A phlyo object, the manipulated tree.
set_branch_lengths_1 = function(phy) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  phy$edge.length <- rep(1, Nedge(phy))
  phy
}


#' Exponentialize branch lengths
#'
#' Sets the deepest branches to length 1/2, the next deepest to 1/4, the next to
#' 1/8, etc.
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @return A phlyo object, the manipulated tree.
set_branch_lengths_exp = function(phy) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  nonroot <- phy$edge[,2]
  
  # Get node depths for non-roots
  phy$edge.length <- rep(1, Nedge(phy))
  depth <- node.depth.edgelength(phy)[nonroot]
  
  # Assign new branch length according to depth
  phy$edge.length <- (1/2) ^ depth
  phy
}


#' Ultrametricize tree by stretching final edges
#'
#' Alters branches ending in a tip in such a way that all tips are equidistant
#' from the root. Does this by lengthening branches above all but the existing,
#' most-distance tip(s).
#' 
#' Identical to \code{\link{ultrametricise}}.
#'
#' @param phy A phlyo object, the tree to manipulate.
#' @return A phlyo object, the manipulated tree.
ultrametricize = function(phy) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
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


#' Ultrametricize tree by stretching final edges
#'
#' Alters branches ending in a tip in such a way that all tips are equidistant
#' from the root. Does this by lengthening branches above all but the existing,
#' most-distance tip(s).
#' 
#' Identical to \code{\link{ultrametricize}}.
#' 
#' @param phy A phlyo object, the tree to manipulate.
#' @return A phlyo object, the manipulated tree.
#' @param phy A phylo object
#' @return A phylo object
ultrametricise = ultrametricize