# Setting branch lengths

#' Set lengths of first branches
#' @param phy A phlyo object
#' @return A phlyo object
set_first_branch_lengths = function(phy, branch_length = 1) {
  
  root <- Ntip(phy) + 1
  first_edges <- which(phy$edge[,1] == root)
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