# Labelling

#' Shorten tree labels to the glottocode
#' 
#' Retain only the glottocode substring
#' within the tip and node labels
#' 
#' @param phy A phylo or multiPhylo object
#' @return A phylo object
label_by_glottocode = function(phy) {
  
  tip_orig <- .get_tip_labels(phy)
  node_orig <- .get_node_labels(phy)
  
  is_multi <- (class(phy) == "multiPhylo")
  if (!is_multi) { phy <- c(phy) }
  
  tip_fail <- node_fail <- logical(0)
  
  for (i in 1:length(phy)) {
    p <- phy[[i]]
    tip_glotto  <- extract_glottocode(p$tip.label)
    node_glotto <- extract_glottocode(p$node.label)
    tip_matched  <- !is.na(tip_glotto)
    node_matched <- !is.na(node_glotto)
    phy[[i]]$tip.label[tip_matched]  <- tip_glotto[tip_matched]
    phy[[i]]$node.label[node_matched] <- node_glotto[node_matched]
    tip_fail <- c(tip_fail, !tip_matched)
    node_fail <- c(node_fail, !node_matched)
  }
  
  if (!is_multi) { phy <- phy[[1]] }
  
  if (any(c(tip_fail, node_fail))) {
    warning(
      str_c("Labels without glottocodes detected and left unchanged: ",
            sum(tip_fail), " tip(s)",
            ifelse(!any(tip_fail), "", 
                   str_c(": ",
                         str_c(str_replace_na(head(tip_orig[tip_fail])),
                               collapse = ", "))),
            ifelse(sum(tip_fail) > 6, "...", ""),
            "; ", sum(node_fail), " node(s)",
            ifelse(!any(node_fail), "", 
                   str_c(": ",
                         str_c(str_replace_na(head(node_orig[node_fail])),
                               collapse = ", "))),
            ifelse(sum(node_fail) > 6, "...", "")
      ))    
  }
  
  phy
}


#' Change glottocode tip and node labels to names
#' 
#' @param phy A phylo object
#' @param glottolog_version A character
#'   string. Which glottolog version to
#'   use.
#' @return A phylo object
label_by_name = function(
  phy,
  glottolog_version = "4.4"
) {
  
  tip_labels <- .get_tip_labels(phy)
  node_labels <- .get_node_labels(phy)
  labels <- c(tip_labels, node_labels)
  labels_df <-
    data.frame(
      label = labels,
      vertex = c(rep("tip", length(tip_labels)),
                 rep("node", length(node_labels))),
      stringsAsFactors = FALSE
    ) %>%
    mutate(
      glottocode = extract_glottocode(labels),
      copy_suffix = .extract_copy_suffix(glottocode),
      glotto_base = .remove_copy_suffix(glottocode)
    ) %>%
    left_join(
      get_glottolog_phylo_geo(glottolog_version) %>%
        select(glotto_base = glottocode, 
               name_base = vertex_name),
      by = "glotto_base"
    ) %>%
    mutate(
      name = str_c(name_base, copy_suffix)
    )
  
  missing_tips <- 
    filter(labels_df, is.na(name) & vertex == "tip")$label
  missing_nodes <- 
    filter(labels_df, is.na(name) & vertex == "node")$label
  n_missing_tips <- length(missing_tips)
  n_missing_nodes <- length(missing_nodes)
  
  # If a name is missing, give warning and retain 
  # the glottocode
  if (n_missing_tips + n_missing_nodes > 0) {
    warning(
      str_c("Labels without glottocodes detected and left unchanged: ",
            n_missing_tips, " tip(s)",
            ifelse(n_missing_tips == 0, "", 
                   str_c(": ",
                         str_c(str_replace_na(head(missing_tips)),
                               collapse = ", "))),
            ifelse(n_missing_tips > 6, "...", ""),
            "; ", n_missing_nodes, " node(s)",
            ifelse(n_missing_nodes == 0, "", 
                   str_c(": ",
                         str_c(str_replace_na(head(missing_nodes)),
                               collapse = ", "))),
            ifelse(n_missing_nodes > 6, "...", "")
      ))   
    labels_df <- labels_df %>%
      mutate(name = ifelse(is.na(name), label, name))
  }
  
  is_multi <- (class(phy) == "multiPhylo")
  if (!is_multi) { phy <- c(phy) }
  for (i in 1:length(phy)) {
    tip_matches <- match(phy[[i]]$tip.label, labels_df$label)
    node_matches <- match(phy[[i]]$node.label, labels_df$label)
    phy[[i]]$tip.label <- labels_df$name[tip_matches]
    phy[[i]]$node.label <- labels_df$name[node_matches]
  }
  if (!is_multi) { phy <- phy[[1]] }
  
  phy
}


#' Extract parts from a glottolog tree label
#' @param label A string
#' @return A string
extract_glottocode = function(labels) {
  # note: two glottocodes contain numbers in the intial
  #       four characters: b10b and 3adt
  regex <- "(?<=(^|\\[))([a-z]{4}|b10b|3adt)[0-9]{4}(?=(\\]|-[0-9]{1,3}$|$))"
  g <- str_extract(labels, regex)
  copy_suffix <- .extract_copy_suffix(labels)
  str_c(g, copy_suffix)
}


#' Add copy suffix to duplicated items
#' in a string of vectors
#' @param x A vector of strings
#' @param strip_first A logical, whether
#'   to remove copy suffixes first.
#' @return A vector of strings
.add_copy_suffix = function(
  x,
  strip_first = TRUE
) {
  if (strip_first) {
    x <- .remove_copy_suffix(x)
  }
  data.frame(x = x) %>%
    group_by(x) %>%
    mutate(
      n = n(), 
      suffix = ifelse(n > 1, str_c("-", row_number()), ""),
      x = str_c(x, suffix)
      ) %>%
    .$x
}


#' Remove copy suffixes
#' @param x A vector of strings
#' @return A vector of strings
.remove_copy_suffix = function(x) {
  x <- str_remove(x, "-[0-9]+$")
}


#' Extract copy suffix from a tree label
#' @param label A string
#' @return A string
.extract_copy_suffix = function(label) {
  regex <- "-[0-9]+$"
  e <- str_extract(label, regex)
  e[is.na(e)] <- ""
  e
}


#' Get tip labels from phylo or multiPhylo
#' @param phy A phlyo or multiPhlyo
.get_tip_labels = function(phy) {
  if (class(phy) == "phylo") {
    phy$tip.label
  } else if (class(phy) == "multiPhylo") {
    lapply(phy, function(x) x$tip.label) %>%
      unlist()
  }
}


.get_node_labels = function(phy) {
  if (class(phy) == "phylo") {
    phy$node.label
  } else if (class(phy) == "multiPhylo") {
    lapply(phy, function(x) x$node.label) %>%
      unlist()
  }
}


.name_to_label = function(name) {
  name %>% 
    str_replace_all(", ", "/") %>%
    str_replace_all("\\(", "{") %>%
    str_replace_all("\\)", "}") %>%
    str_remove_all("[' ]")
}