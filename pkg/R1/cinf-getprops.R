# Returns the list of property names from molecular database
get_prop_names <- function(moldbase) {
  prop_count <- list()
  nmol <- length(moldbase)
  for (imol in 1:nmol) {
    pnamelist <- names(moldbase[[imol]]$props)
	for (pname in pnamelist) {
	  if (is.null(prop_count[[pname]])) {
	    prop_count[[pname]] <- 1
	  } else {
	    prop_count[[pname]] <- prop_count[[pname]] + 1
	  }
	}
  }
  names(prop_count)
}

# Returns vector of property values
get_prop_vec <- function(moldbase, propname) {
  prop_vec <- numeric()
  nmol <- length(moldbase)
  for (imol in 1:nmol) {
    if (!is.null(moldbase[[imol]]$props[[propname]])) {
	  prop_vec[imol] <- moldbase[[imol]]$props[[propname]]
	} else {
	  prop_vec[imol] <- NA
	}
  }
  prop_vec
}

# Returns a data frame for the specified property names
get_props <- function(moldbase, propnames=get_prop_names(moldbase)) {
  pval_list <- list()
  for (pname in propnames) {
    prop_vec <- get_prop_vec(moldbase, pname)
	pval_list[[pname]] <- prop_vec
  }
  pval_df <- as.data.frame(pval_list)
  pval_df
}

# Returns a data frame for a specified property name
get_prop <- function(moldbase, propname) get_props(moldbase, list(propname))

