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
  if (!is.data.frame(data)) {
    stop("data must be a data.frame.")
  }
  
  # Convert phy to multiPhlyo
  phy_original <- phy
  if (class(phy) == "phylo") { phy <- c(phy) }
  n_tree <- length(phy)
  
  # Check data contents
  col_name <- colnames(data)
  if (!("language" %in% col_name)) {
    stop("data must contain a column named 'language'.")
  }
  
  is_num_col <- 
    sapply(1:ncol(data), 
           function(i) is.numeric(data[,i]))
  if (sum(is_num_col) == 0) {
    stop("data must contain at least numeric column.")
  }
  
  is_extra_col <- !is_num_col & (col_name != "language")
  if (any(is_extra_col)) {
    warning(
      str_c("data contains non-numeric columns other than",
            " 'language', which have been ignored: ",
            str_c(head(col_name[is_extra_col]), collapse = ", "),
            ifelse(sum(is_extra_col) > 6, "...", "")
            )
      )
  }
  data$.tip.label = .add_copy_suffix(data$language)
  
  # Check phy contents
  for (i in 1:n_tree) {
    if (!("edge.length" %in% names(phy[[i]]))) {
      stop(paste("phy[[", i , "]] lacks branch lengths.", 
                 sep = ""))
    } else if(any(is.na(phy[[i]]$edge.length))) {
      stop(paste("phy[[", i , "]] has at least one ",
                 "undefined branch length.", sep = ""))       
    }
    if (Ntip(phy[[i]]) != nrow(data)) {
      stop(
        str_c("Number of tips in phy[[", i, "]] ",
              "differs from the number of rows in data")
      )
    }
    phy[[i]]$tip.label <- .add_copy_suffix(phy[[i]]$tip.label)
    if (any(sort(data$.tip.label) != sort(phy[[i]]$tip.label))) {
      stop(
        str_c("Tip labels do not match data$language ",
              "in phy[[", i, "]]")
      )
    }
    nonsister_copies <- which_nonsister_copies(phy[[i]])
    if (length(nonsister_copies) > 0) {
      stop(
        str_c("Trees can contain duplicate tip labels, ",
              "but only if the tips are sisters with ",
              "identical branch lengths above them. This ",
              "is not true in phy[[", i, "]] for sisters: ",
              str_c(sort(nonsister_copies), collapse = ", "))
      )
    }
  }
  
  # Get weights
  ACL_weights <- BM_weights <- data %>% select(.tip.label)
  for (i in 1:length(phy)) {
    ACL_i <- data.frame(x = ACL(phy[[i]]), 
                        .tip.label = phy[[i]]$tip.label)
    BM_i  <- data.frame(x = BM(phy[[i]]),
                        .tip.label = phy[[i]]$tip.label)
    colnames(ACL_i)[1] <- colnames(BM_i)[1] <-  str_c("tree", i)
    ACL_weights <- left_join(ACL_weights, ACL_i, by = ".tip.label")
    BM_weights  <- left_join(BM_weights,  BM_i, by = ".tip.label")
  }
  
  # Calculate averages using matrix multiplication
  data_mat <- t(as.matrix(data[, is_num_col]))
  ACL_mat <- as.matrix(ACL_weights %>% select(-.tip.label))
  BM_mat  <- as.matrix(BM_weights  %>% select(-.tip.label))
  ACL_ave_mat <- t(data_mat %*% ACL_mat)
  BM_ave_mat  <- t(data_mat %*% BM_mat)
  
  # Prepare results as dataframes
  ACL_averages <- 
    as.data.frame(ACL_ave_mat) %>%
    mutate(tree = str_c("tree", 1:n_tree)) %>%
    select(tree, everything())
  BM_averages <- 
    as.data.frame(BM_ave_mat) %>%
    mutate(tree = str_c("tree", 1:n_tree)) %>%
    select(tree, everything())
  rownames(ACL_averages) <- 
    rownames(BM_averages) <- 
    rownames(ACL_weights) <-
    rownames(BM_weights) <-
    NULL
  ACL_weights <- 
    left_join(data[, !is_num_col], 
              ACL_weights, by = ".tip.label") %>%
    select(-.tip.label)
  BM_weights <- 
    left_join(data[, !is_num_col], 
              BM_weights, by = ".tip.label") %>%
    select(-.tip.label)
  
  # Return results
  list(
    phy  = phy_original,
    data = data %>% select(-.tip.label),
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