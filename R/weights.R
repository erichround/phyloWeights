# Phylogenetic weights

#' Phylogenetically-weighted averages
#'
#' Calculates phylogenetic weights according to the ACL method (see \link{ACL})
#' and BranchManager method (see \link{BM}), for one or more reference trees,
#' and applies the resulting sets of weights to one or more numerical variables
#' which characterize the languages in the trees, resulting in
#' phylogenetically-sensitive (i.e., phylogenetically-weighted) averages. If the
#' data is in the form of binary \{0,1\} values, this is equivalent to a
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
#'   column \code{tip}, whose contents match the tip labels in \code{phy}, and
#'   at least one column of numerical data.
#' @param na.rm A logical. Whether to ignore languages with NA values in the
#'   data. If TRUE then for each column of data, tips with NA values are
#'   ignored and the weights of the remaining tips are re-scaled to sum to 1.
#' @return A list, containing phy, data, and dataframes: ACL_weights,
#'   BM_weights, ACL_averages, BM_averages
#' @examples
#'
#' library(ape)
#'
#' # A dataframe of language data, with columns: tip, is_SOV, n_consonants
#' dat <- data.frame(
#'   tip = c("A", "B", "C", "D"),
#'   is_SOV = c(1, 1, 1, 0),
#'   n_consonants = c(18, 20, 22, 40))
#' # Three trees, representing three genealogical hypotheses, that relate
#' # the same languages A,B,C,D:
#' trees <- c(
#'   read.tree(text = "(((A:1,B:1,C:1):1,D:2):0.3);"),
#'   read.tree(text = "(((A:0.2,B:0.2,C:0.2):1.8,D:2):0.3);"),
#'   read.tree(text = "(((A:1.8,B:1.8,C:1.8):0.2,D:2):0.3);"))
#' p_ave <- phylo_average(phy = trees, data = dat)
#' # Phylogenetic weights, calculated from the trees:
#' p_ave$ACL_weights
#' p_ave$BM_weights
#' # Phylogenetically-weighted averages, or, for binary data, phylogenetically-
#' # weighted proportions
#' p_ave$ACL_averages
#' p_ave$BM_averages
phylo_average = function(phy, data, na.rm = FALSE) {

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
  # Check the tip column
  if (!("tip" %in% col_name)) {
    stop("Can't find column `tip` in `data`.")
  } else {
    if (is.factor(data$tip)) {
      data$tip <- as.character(data$tip)
    } else {
      cl <- class(data$tip) 
      if (!is.character(data$tip)) {
        stop(str_c("Column `tip` in `data` must contain character ",
                   "strings.\n You supplied a column of class ", cl, "."))
      }
    }
    is_dupl <- duplicated(data$tip) 
    if (any(is_dupl)) {
      stop(
        str_c("`Column `tip` in `data` must not contain duplicates.\n",
              "You supplied duplicates of: ",
              str_c(head(unique(data$tip[is_dupl]), 4), collapse = ", "),
              ifelse(sum(is_dupl) > 4, "..", ""), ".")
      )
    }
  }
  # Check the other columns
  is_num_col <- 
    sapply(1:ncol(data), 
           function(i) is.numeric(as.data.frame(data)[,i]))
  if (sum(is_num_col) == 0) {
    stop(str_c("`data` must contain at least one numeric column.\n",
               "You supplied a data frame with no numeric columns."))
  }
  is_extra_col <- !is_num_col & (col_name != "tip")
  if (any(is_extra_col)) {
    warning(
      str_c("`data` contains non-numeric columns other than",
            " `tip`, which have been ignored: ",
            str_c(head(col_name[is_extra_col], 4), collapse = ", "),
            ifelse(sum(is_extra_col) > 4, "..", ""), ".")
    )
  }
  data_numeric <- data %>% select(which(is_num_col))
  data_nonnumeric <- data %>% select(which(!is_num_col))
  
  # Convert phy to multiPhlyo
  is_phylo <- class(phy) == "phylo"
  phy_original <- phy
  if (is_phylo) { phy <- c(phy) }
  n_tree <- length(phy)
  
  
  # Check phy's contents and their match with data
  for (i in 1:n_tree) {
    phy_str <- ifelse(is_phylo, "`phy`", str_c("`phy[[", i, "]]`"))
    Each_phy <- ifelse(is_phylo, "`phy`", "Each tree in `phy`")
    each_phy <- ifelse(is_phylo, "`phy`", "each tree in `phy`")
    
    # Check presence of branch lengths
    if (!("edge.length" %in% names(phy[[i]]))) {
      stop(str_c(Each_phy, " must have branch lengths.\n",
                 phy_str, " lacks branch lengths."))
    } else if(any(is.na(phy[[i]]$edge.length))) {
      stop(str_c(Each_phy, " must have lengths for all branches.\n", 
                 phy_str, " has at least one undefined branch length."))       
    }
    # Check tip labels and their matches
    is_dupl <- duplicated(phy$tip.label) 
    dupl <- unique(phy$tip.label[is_dupl])
    if (any(is_dupl)) {
      stop(
        str_c("In ", each_phy, ", all tips must have distinct labels.\n",
              "In ", phy_str, " you supplied duplicates of: ",
              str_c(head(dupl, 4), collapse = ", "),
              ifelse(length(dupl) > 4, "..", ""), ".")
      )
    }
    extra <- setdiff(phy$tip.label, data$tip)
    if (length(extra) > 0) {
      stop(
        str_c("In ", each_phy, ", all tips must be represented in ",
              "the `tip` column of `data`.\n",
              "There is no match in `data` for: ",
              str_c(head(extra, 4), collapse = ", "),
              ifelse(length(extra) > 4, "...", ""), " in ", phy_str, ".")
      )
    }
    missing <- setdiff(phy$tip.label, data$tip)
    if (length(missing) > 0) {
      stop(
        str_c(Each_phy, " must have a tip for every value in ",
              "the `tip` column of `data`.\n",
              "In ", phy_str, "there is no tip that matches: ",
              str_c(head(missing, 4), collapse = ", "),
              ifelse(length(missing) > 4, "..", ""), ".")
      )
    }
  }
  
  # Get weights
  ACL_weights <- BM_weights <- data %>% select(tip)
  for (i in 1:length(phy)) {
    ACL_i <- data.frame(x = ACL(phy[[i]]), tip = phy[[i]]$tip.label)
    BM_i  <- data.frame(x = BM(phy[[i]]),  tip = phy[[i]]$tip.label)
    colnames(ACL_i)[1] <- colnames(BM_i)[1] <-  str_c("tree", i)
    ACL_weights <- left_join(ACL_weights, ACL_i, by = "tip")
    BM_weights  <- left_join(BM_weights,  BM_i,  by = "tip")
  }
  
  # Calculate averages 
  data_mat <- t(as.matrix(data_numeric))
  ACL_mat <- as.matrix(ACL_weights %>% select(-tip))
  BM_mat  <- as.matrix(BM_weights  %>% select(-tip))
  # Use matrix multiplication
  ACL_ave_mat <- t(data_mat %*% ACL_mat)
  BM_ave_mat  <- t(data_mat %*% BM_mat)
  has_na <- apply(data_mat, 1, function(i) any(is.na(i)))
  if (na.rm & any(has_na)) {
    # Use loops for data columns with NAs
    for (d in which(has_na)) {
      dat_d <- data_mat[d, ]
      is_good <- !is.na(dat_d)
      if (any(is_good)) {
        dat_good <- dat_d[is_good]
        for (p in 1:n_tree) {
          wt_ACL_good <- ACL_mat[, p][is_good]
          wt_BM_good <- BM_mat[, p][is_good]
          ACL_ave_mat[p, d] <- sum(dat_good * wt_ACL_good) / sum(wt_ACL_good)
          BM_ave_mat[p, d] <- sum(dat_good * wt_BM_good) / sum(wt_BM_good)
        }
      }
    }
  }
  
  # Prepare results as dataframes
  ACL_averages <- 
    as.data.frame(ACL_ave_mat) %>%
    mutate(tree = str_c("tree", 1:n_tree)) %>%
    select(tree, everything())
  BM_averages <- 
    as.data.frame(BM_ave_mat) %>%
    mutate(tree = str_c("tree", 1:n_tree)) %>%
    select(tree, everything())
  colnames(ACL_averages) <- colnames(BM_averages) <-
    c("tree", colnames(data_numeric))
  rownames(ACL_averages) <- rownames(BM_averages) <- 
    rownames(ACL_weights) <- rownames(BM_weights) <- NULL
  ACL_weights <- left_join(data_nonnumeric, ACL_weights, by = "tip") 
  BM_weights  <- left_join(data_nonnumeric, BM_weights,  by = "tip")
  
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
#'
#' Calculates phylogenetic weights (a.k.a. phylogenetic means) according to the
#' classic method of Altschul, Carroll & Lipman (1987).
#'
#' @param phy A phylo object.
#' @return A vector of weights.
#' @examples 
#' 
#' library(ape)
#' 
#' # Three trees, representing three genealogical hypotheses, that relate
#' # languages A,B,C,D:
#' tree1 <- read.tree(text = "(((A:1,B:1,C:1):1,D:2):0.3);")
#' plot(tree1)
#' tree2 <- read.tree(text = "(((A:0.2,B:0.2,C:0.2):1.8,D:2):0.3);")
#' plot(tree2)
#' tree3 <- read.tree(text = "(((A:1.8,B:1.8,C:1.8):0.2,D:2):0.3);")
#' plot(tree2)
#' # The ACL weights of each language in each tree:
#' ACL(tree1)
#' ACL(tree2)
#' ACL(tree3)
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
#' Calculates phylogenetic weights according the to BranchManager method of
#' Stone & Sidow (2007).
#'
#' Adapted from original R code by Eric A. Stone, published as supplementary
#' material to Stone & Sidow (2007).
#'
#' @param phy A phylo object
#' @return A vector of weights
#' @examples 
#' 
#' library(ape)
#' 
#' # Three trees, representing three genealogical hypotheses, that relate
#' # languages A,B,C,D:
#' tree1 <- read.tree(text = "(((A:1,B:1,C:1):1,D:2):0.3);")
#' plot(tree1)
#' tree2 <- read.tree(text = "(((A:0.2,B:0.2,C:0.2):1.8,D:2):0.3);")
#' plot(tree2)
#' tree3 <- read.tree(text = "(((A:1.8,B:1.8,C:1.8):0.2,D:2):0.3);")
#' plot(tree2)
#' # The BM weights of each language in each tree:
#' BM(tree1)
#' BM(tree2)
#' BM(tree3)
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
    # Here, the function .edgeroot() has been re-implemented to repair a 
    # dependency on an old version of the ape package.
    W <- W + curedgelen*ACL(.edgeroot(i, phy))
    totallength <- totallength + curedgelen
  }
  W <- W/totallength
  W_order <- match(phy$tip.label, names(W))
  W[W_order]
}


#' Edgeroot
#'
#' The function edgeroot takes a rooted tree (input phy) and reroots that tree
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
  
  # Insert a new node with a child tip "temp_tip", appended from the mid- point
  # of edge x
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