# Phylogenetic weights

#' Phylogenetically-weighted averages
#'
#' Calculates phylogenetic weights according to the ACL method (see \link{ACL})
#' and BranchManager method (see \link{BM}), for one or more reference trees,
#' and applies the resulting sets of weights to one ore more numerical variables
#' which characterize the languages in the trees, resulting in
#' phylogenetically-sensitive (i.e., phylogenetically-weighted) averages. If the
#' data is in the form of binary {0,1} values, this is equivalent to a
#' phylogenetically-sensitive proportion.
#'
#' The function returns a list, containing: \code{phy}, the input tree(s);
#' \code{data}, the input dataframe; and dataframes \code{ACL_weights},
#' \code{BM_weights}, \code{ACL_averages}, \code{BM_averages}. Dataframe
#' \code{ACL_weights} contains one column for each tree in \code{phy}. In those
#' columns are the phylogenetic weights obtained using the ACL method.
#' Additionally, \code{ACL_weights} contains all non-numeric columns of the
#' input \code{data} dataframe. Dataframe \code{BM_weights} is similar but with
#' phylogenetic weights obtained using the BM method. Dataframe
#' \code{ACL_averages} has one row per tree and one column for each numerical
#' column in \code{data}, and contains phylogenetically-sensitive averages
#' obtained using the ACL method. Dataframe \code{BM_averages} is similar but
#' with phylogenetically-sensitive averages obtained using the BM method.
#'
#' @param phy A phylo or multiPhylo object, containing one or more trees,
#'   representing one or more hypotheses on the phylogenetic relatedness of a
#'   set of languages.
#' @param data A dataframe, containing data on the languages. Must contain a
#'   column \code{language}, whose contents match the tip labels in \code{phy},
#'   and at least one column of numerical data.
#' @return A list, containing phy, data, and dataframes: ACL_weights,
#'   BM_weights, ACL_averages, BM_averages
phylo_average = function(
  phy = NULL,
  data = NULL
) {

  # Check phy class
  if (class(phy) != "phylo" & 
      class(phy) != "multiPhylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo or multiPhylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Check data
  if (!is.data.frame(data)) {
    cd <- class(data)
    stop(str_c("`data` must be a dataframe.\n",
               "You supplied an object of class ", cd, "."))
  }
  col_name <- colnames(data)
  if (!("language" %in% col_name)) {
    stop("Can't find column `language` in `data`.")
  } else {
    if (is.factor(data$language)) {
      data$language <- as.character(data$language)
    } else {
      cl <- class(data$language) 
      if (!is.character(data$language)) {
        stop(str_c("Column `language` in `data` must contain character strings.\n",
                   "You supplied a column of class ", cl, "."))
      }
    }
  }
  is_num_col <- 
    sapply(1:ncol(data), 
           function(i) is.numeric(as.data.frame(data)[,i]))
  if (sum(is_num_col) == 0) {
    stop(str_c("`data` must contain at least one numeric column.\n",
               "You supplied a data frame with no numeric columns."))
  }
  is_extra_col <- !is_num_col & (col_name != "language")
  if (any(is_extra_col)) {
    warning(
      str_c("`data` contains non-numeric columns other than",
            " `language`, which have been ignored: ",
            str_c(head(col_name[is_extra_col]), collapse = ", "),
            ifelse(sum(is_extra_col) > 6, "...", "")
      )
    )
  }
  
  # Convert phy to multiPhlyo
  phy_original <- phy
  if (class(phy) == "phylo") { phy <- c(phy) }
  n_tree <- length(phy)
  
  # Add column of tip copy suffixes
  data$.tip.label = .add_copy_suffix(data$language)
  
  # Check phy's contents and match with data
  for (i in 1:n_tree) {
    if (!("edge.length" %in% names(phy[[i]]))) {
      stop(str_c("All trees in `phy` must have branch lengths.\n",
                 "`phy[[", i , "]]` lacks branch lengths."))
    } else if(any(is.na(phy[[i]]$edge.length))) {
      stop(str_c("Trees in `phy` must have lengths for all branches.\n", 
                 "`phy[[", i , "]]` has at least one ",
                 "undefined branch length."))       
    }
    if (Ntip(phy[[i]]) != nrow(data)) {
      stop(str_c("The number of tips for all trees in `phy` must ",
                 "match the number of rows in `data`.",
                 "Number of tips in `phy[[", i, "]]` ",
                 "is ", Ntip(phy[[i]]), " but the number of rows ",
                 "in `data` is ", nrow(data), ".")
      )
    }
    phy[[i]]$tip.label <- .add_copy_suffix(phy[[i]]$tip.label)
    if (any(sort(data$.tip.label) != sort(phy[[i]]$tip.label))) {
      stop(str_c("Tip labels in all trees in `phy` must match ",
                 "the contents of the `language` column of `data`.\n",
                 "Tip labels do not match `data$language` ",
                 "in `phy[[", i, "]]`.")
      )
    }
    nonsister_copies <- .which_nonsister_copies(phy[[i]])
    if (length(nonsister_copies) > 0) {
      stop(str_c("Trees in `phy` can contain duplicate tip labels, ",
                 "but only if the tips are sisters with ",
                 "identical branch lengths above them.\n",
                 "This is not true in `phy[[", i, "]]` for sisters: ",
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
  colnames(ACL_averages) <-
    colnames(BM_averages) <-
    c("tree", colnames(data)[is_num_col])
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
#' 
#' Calculates phylogenetic weights (a.k.a. phylogenetic means) according to the classic method of Altschul, Carroll & Lipman DJ (1987).
#' 
#' @param phy A phylo object.
#' @return A vector of weights.
ACL = function(phy) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  S <- vcv.phylo(phy)
  rowSums(solve(S))/sum(solve(S))
}


#' BranchManager weights
#' 
#' Calculates phylogenetic weights according the to BranchManager method of Stone & Sidow (2007).
#' 
#' Adapted from original R code by Eric A. Stone available at \url{https://static-content.springer.com/esm/art%3A10.1186%2F1471-2105-8-222/MediaObjects/12859_2007_1594_MOESM3_ESM.txt}.
#' The function \code{edgeroot()} has been re-implemented to repair a dependency on an old version of the ape package.
#' 
#' @param phy A phylo object
#' @return A vector of weights
BM <- function(phy) {
  
  # Check phy
  if (class(phy) != "phylo") {
    cp <- class(phy)
    stop(str_c("`phy` must be of class phylo.\n",
               "You supplied an object of class ", cp, "."))
  }
  
  # Dichotomise the tree
  phy <- multi2di(phy)
  
  # Remove interior node labels
  phy$node.label <- rep("", length(phy$node.label))
  
  ###

  numedges <- dim(phy$edge)[1]
  numtaxa <- length(phy$tip.label)
  W <- rep(0,numtaxa)
  totallength <- 0
  for (i in 1:numedges) {
    curedgelen <- phy$edge.length[i]
    W <- W + curedgelen*ACL(.edgeroot(i, phy))
    totallength <- totallength + curedgelen
  }
  W <- W/totallength
  W
}


#' Edgeroot
#' 
#' The function edgeroot takes a rooted 
#' tree (input phy) and reroots that tree 
#' at the midpoint of the specified edge 
#' 
#' @param x An integer, specifying an egde.
#' @param phy A phylo object, the tree to manipulate.
#' @return A phylo object, the manipulated tree. 
#' @noRd
.edgeroot = function(x, phy) {
  
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