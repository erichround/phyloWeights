# Phylogenetic weights

#' Phylogenetically-weighted averages
#' @param phy A phylo or multiPhylo object
#' @param data A dataframe
#' @return A list, containing phy, data,
#'   and dataframes: ACL_weights, BM_weights,
#'   ACL_averages, BM_averages
phylo_average = function(
  phy = NULL,
  data = NULL
) {

  # Check classes
  if (!(class(phy) %in% c("phylo", "multiPhylo"))) {
    stop("phy must be of class phylo or multiPhylo.")
  }
  if (!is.data.frame(data_Fig2)) {
    stop("data must be a data.frame.")
  }
  
  
  # Convert phy to multiPhlyo
  phy_original <- phy
  if (class(phy) == "phylo") { phy <- c(phy) }
  n_tree <- length(phy)
  
  
  # Check phy contents
  for (i in 1:n_tree) {
    if (!("edge.length" %in% names(phy[[i]]))) {
      stop(paste("phy[[", i , "]] lacks branch lengths.", 
                 sep = ""))
    } else if(any(is.na(phy[[i]]$edge.length))) {
      stop(paste("phy[[", i , "]] has at least one ",
                 "undefined branch length.", sep = ""))       
    }
  }
  
  
  # Check data contents
  if (!("language" %in% colnames(data))) {
    stop("data must contain a column named 'language'.")
  }
  data$language <- as.character(data$language)
  if (any(duplicated(data$language))) {
    stop("data contains duplicate labels in column 'language'.")
  }
  data <- data %>% arrange(language)
  if (ncol(data) < 2) {
    stop("data must contain at least one column other than 'language'.")
  }
  for (i in 1:ncol(data)) {
    if (colnames(data)[i] == "language") { next }
    if (!is.numeric(data[, i])) {
      stop(paste("Column", i , "of data is not ",
                 "numeric but should be.", sep = ""))
    }
  }
  
  
  # Check that languages match if phy and data
  for (i in 1:n_tree) {
    tips <- phy[[i]]$tip.label
    missing_tips  <- setdiff(data$language, tips)
    missing_langs <- setdiff(tips, data$language)
    if (length(missing_tips) > 0) {
      stop(paste("The following languages are in ",
                 "data but not in phy[[", i , "]]: ",
                 paste(missing_tips, collapse = ", "), 
                 sep = ""))
    }
    if (length(missing_langs) > 0) {
      stop(paste("The following languages are in phy[[", 
                 i , "]] but not in data: ",
                 paste(missing_langs, collapse = ", "), 
                 sep = ""))
    }
  }
  
  
  # Get weights
  ACL_weights <- BM_weights <-
    data.frame(language = tips, stringsAsFactors = FALSE)
  for (i in 1:length(phy)) {
    ACL_i <- data.frame(x = ACL(phy[[i]]))
    BM_i  <- data.frame(x = BM(phy[[i]]))
    colnames(ACL_i) <- colnames(BM_i) <- 
      paste("tree", i, sep = "")
    ACL_weights <- bind_cols(ACL_weights, ACL_i)
    BM_weights  <- bind_cols(BM_weights,  BM_i)
  }
  
  
  # Calculate averages using matrix multiplication
  data_mat <- t(as.matrix(data %>% select(-language)))
  ACL_mat <- as.matrix(ACL_weights %>% select(-language))
  BM_mat  <- as.matrix(BM_weights  %>% select(-language))
  ACL_ave_mat <- t(data_mat %*% ACL_mat)
  BM_ave_mat  <- t(data_mat %*% BM_mat)
  
  # Prepare results as dataframes
  ACL_averages <- 
    as.data.frame(ACL_ave_mat) %>%
    mutate(tree = paste("tree", 1:n_tree, sep = "")) %>%
    select(tree, everything())
  BM_averages <- 
    as.data.frame(BM_ave_mat) %>%
    mutate(tree = paste("tree", 1:n_tree, sep = "")) %>%
    select(tree, everything())
  rownames(ACL_averages) <- 
    rownames(BM_averages) <- 
    rownames(ACL_weights) <-
    rownames(BM_weights) <-
    NULL
  
  # Return results
  list(
    phy  = phy_original,
    data = data,
    ACL_weights  = ACL_weights,
    BM_weights   = BM_weights,
    ACL_averages = ACL_averages,
    BM_averages  = BM_averages
  )
}


#' ACL weights
#' @param phy A phylo object.
#' @return A vector of weights
ACL = function(phy) {
  if (class(phy) != "phylo") 
    stop("object \"phy\" is not of class \"phylo\"")
  
  S <- vcv.phylo(phy)
  rowSums(solve(S))/sum(solve(S))
}


#' BranchMaster weights
#' @param phy A phylo object
#' @return A vector of weights
BM <- function(phy) {
  if (class(phy) != "phylo") 
    stop("object \"phy\" is not of class \"phylo\"")
  
  # Dichotomise the tree
  phy <- multi2di(phy)
  
  # Remove interior node labels
  phy$node.label <- rep("", length(phy$node.label))
  
  ###
  # From https://static-content.springer.com/esm/art%3A10.1186%2F1471-2105-8-222/MediaObjects/12859_2007_1594_MOESM3_ESM.txt
  # However, the function edgeroot() has been re-implemented to repair
  # a dependency on an old version of the ape package.
  numedges <- dim(phy$edge)[1]
  numtaxa <- length(phy$tip.label)
  W <- rep(0,numtaxa)
  totallength <- 0
  for (i in 1:numedges) {
    curedgelen <- phy$edge.length[i]
    W <- W + curedgelen*ACL(edgeroot(i, phy))
    totallength <- totallength + curedgelen
  }
  W <- W/totallength
  ###
  
  W[order(names(W))]
}


#' Edgeroot
#' The function edgeroot takes a rooted 
#' tree (input phy) and reroots that tree 
#' at the midpoint of the specified edge 
#' @param x An integer, specifying an egde.
#' @param phy A phylo object.
edgeroot = function(x, phy) {
  
  # Add root.edge to phy if it is missing
  if(is.null(phy$root.edge)) {
    phy$root.edge <- setdiff(phy$edge[,1], phy$edge[,2])
  }
  
  # Insert a new node with a child
  # tip "temp_tip", appended from the mid-
  # point of edge x
  tipward_node <- phy$edge[x,2]
  current_edge_length <- phy$edge.length[x]
  newtip_phy <- 
    bind.tip(phy,
             tip.label = "temp_tip",
             edge.length = 1, # Arbitrary, since the edge will soon be removed
             where = tipward_node,
             position = current_edge_length/2)
  newtip <- which(newtip_phy$tip.label == "temp_tip")
  
  # Reroot, to place the new tip as an outlier
  rerooted_phy <-
    phytools::reroot(newtip_phy, newtip)
  
  # Drop that outlier and return the phylogeny
  newtip <- which(rerooted_phy$tip.label == "temp_tip")
  drop.tip(rerooted_phy, newtip)
}