# Grid operations in CMF

# Returns dimensions of molecule
cmf_moldim <- function(mol) {
  moldim <- list()
  natoms <- length(mol$atoms)
  if (natoms > 0) {
    atom <- mol$atoms[[1]]
    moldim$xmin <- atom$x
    moldim$xmax <- atom$x
    moldim$ymin <- atom$y
    moldim$ymax <- atom$y
    moldim$zmin <- atom$z
    moldim$zmax <- atom$z
    if (natoms > 1) {
      for (iatom in 2:natoms) {
        atom <- mol$atoms[[iatom]]
        if (atom$x < moldim$xmin) moldim$xmin <- atom$x
        if (atom$x > moldim$xmax) moldim$xmax <- atom$x
        if (atom$y < moldim$ymin) moldim$ymin <- atom$y
        if (atom$y > moldim$ymax) moldim$ymax <- atom$y
        if (atom$z < moldim$zmin) moldim$zmin <- atom$z
        if (atom$z > moldim$zmax) moldim$zmax <- atom$z
      }
    }
  }
  moldim
}

# Initializes grid for a given step and margin around molecules from database
cmf_init_grid <- function(mdb, step=1.0, margin=2.0) {
  grid <- list()
  nmols <- length(mdb)
  if (nmols > 0) {
    mdbdim <- cmf_moldim(mdb[[1]])
    if (nmols > 1) {
      for (imol in 2:nmols) {
        mol <- mdb[[imol]]
        moldim <- cmf_moldim(mol)
        if (moldim$xmin < mdbdim$xmin) mdbdim$xmin <- moldim$xmin
        if (moldim$xmax > mdbdim$xmax) mdbdim$xmax <- moldim$xmax
        if (moldim$ymin < mdbdim$ymin) mdbdim$ymin <- moldim$ymin
        if (moldim$ymax > mdbdim$ymax) mdbdim$ymax <- moldim$ymax
        if (moldim$zmin < mdbdim$zmin) mdbdim$zmin <- moldim$zmin
        if (moldim$zmax > mdbdim$zmax) mdbdim$zmax <- moldim$zmax
      }
    }
    mdbdim$xmin <- mdbdim$xmin - margin
    mdbdim$xmax <- mdbdim$xmax + margin
    mdbdim$ymin <- mdbdim$ymin - margin
    mdbdim$ymax <- mdbdim$ymax + margin
    mdbdim$zmin <- mdbdim$zmin - margin
    mdbdim$zmax <- mdbdim$zmax + margin

    sizex <- mdbdim$xmax - mdbdim$xmin
    grid$ngridx <- ceiling(sizex / step) + 1
    grid$gridx <- numeric(grid$ngridx)
    grid$gridx[1] <- mdbdim$xmin
    grid$gridx[grid$ngridx] <- mdbdim$xmax
    if (grid$ngridx > 2) {
      for (igrid in 1:(grid$ngridx-2)) {
        grid$gridx[igrid+1] <- mdbdim$xmin + igrid * step
      }
    }

    sizey <- mdbdim$ymax - mdbdim$ymin
    grid$ngridy <- ceiling(sizey / step) + 1
    grid$gridy <- numeric(grid$ngridy)
    grid$gridy[1] <- mdbdim$ymin
    grid$gridy[grid$ngridy] <- mdbdim$ymax
    if (grid$ngridy > 2) {
      for (igrid in 1:(grid$ngridy-2)) {
        grid$gridy[igrid+1] <- mdbdim$ymin + igrid * step
      }
    }

    sizez <- mdbdim$zmax - mdbdim$zmin
    grid$ngridz <- ceiling(sizez / step) + 1
    grid$gridz <- numeric(grid$ngridz)
    grid$gridz[1] <- mdbdim$zmin
    grid$gridz[grid$ngridz] <- mdbdim$zmax
    if (grid$ngridz > 2) {
      for (igrid in 1:(grid$ngridz-2)) {
        grid$gridz[igrid+1] <- mdbdim$zmin + igrid * step
      }
    }

    grid$val <- array(0.0, c(grid$ngridx, grid$ngridy, grid$ngridz))
  }
  grid
}


