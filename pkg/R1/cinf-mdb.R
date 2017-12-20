# Manipulations with mdb (molecular data base)

# Get the number of compounds in mdb
mdb_get_num_comp <- function(mdb) {
  length(mdb)
}

# Returns complete list of properties in mdn
mdb_get_prop_names <- function(mdb) {
  glob_prop_names <- list()
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    loc_prop_names <- names(mdb[[imol]]$props)
    glob_prop_names <- unique(c(glob_prop_names, loc_prop_names))
  }
  glob_prop_names
}

# Keep in mdb only compounds containing values of certain property
mdb_keep_with_prop <- function(mdb, prop_to_keep) {
  new_mdb <- list()
  imol1 <- 0
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    prop_names <- names(mdb[[imol]]$props)
    if (prop_to_keep %in% prop_names) {
      imol1 <- imol1 + 1
      new_mdb[imol1] <- mdb[imol]
    }
  }
  new_mdb
}

# Keep in matrix only rows and columns corresponding to
# compounds containing values of certain property
mdb_keep_matr_with_prop <- function(matr, mdb, prop_to_keep) {
  nmols <- length(mdb)
  to_keep <- rep(FALSE, nmols)
  for (imol in 1:nmols) {
    prop_names <- names(mdb[[imol]]$props)
    if (prop_to_keep %in% prop_names) {
      to_keep[imol] <- TRUE
    }
  }
  matr[to_keep, to_keep]
}

# Extracts property vector by name
mdb_get_prop_vect <- function(mdb, propname) {
  ncomp <- length(mdb)
  propvect <- array(0, ncomp)
  for (imol in 1:ncomp) {
    propvect[imol] <- mdb[[imol]]$props[[propname]]
  }
  propvect
}

# Extracts property one-column matrix by name
mdb_get_prop_matr1 <- function(mdb, propname) {
  ncomp <- length(mdb)
  propmatr1 <- matrix(0, ncomp, 1)
  for (imol in 1:ncomp) {
    propmatr1[imol,1] <- mdb[[imol]]$props[[propname]]
  }
  propmatr1
}

# Transfers properties from one mdb to another
transfer_props <- function(mdb_target, mdb_source) {
  nmols <- length(mdb_target)
  for (imol in 1:nmols) {
    mdb_target[[imol]]$props <- mdb_source[[imol]]$props
  }
  mdb_target
}


