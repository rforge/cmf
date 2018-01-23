# Displays molecule

require(rgl)

#source("cinf-ptable.R")
load("data/PT.RData")

mol_view_cpk <- function(mol, reset=TRUE, alpha=1, rfactor=0) {
  natoms <- length(mol$atoms)
  x <- double(natoms)
  y <- double(natoms)
  z <- double(natoms)
  colors <- double(natoms)
  radius <- double(natoms)
  for (i in 1:natoms) {
    atom <- mol$atoms[[i]]
    x[i] <- atom$x
    y[i] <- atom$y
    z[i] <- atom$z
    colors[i] <- PT.Color[[atom$el]]
    radius[i] <- PT.AtRad[[atom$el]] * rfactor
  }
  if (reset) {
    open3d()
    bg3d("white")
  }
  spheres3d(x, y, z, color=colors, radius=radius, alpha=alpha)
}

mol_view_lines <- function(mol) {
  x <- double(2)
  y <- double(2)
  z <- double(2)
  nbonds <- length(mol$bonds)
  for (ibond in 1:nbonds) {
    bond <- mol$bonds[[ibond]]
    at1 <- bond$at1
    at2 <- bond$at2
    atom1 <- mol$atoms[[at1]]
    atom2 <- mol$atoms[[at2]]
    mx <- (atom1$x + atom2$x) / 2
    my <- (atom1$y + atom2$y) / 2
    mz <- (atom1$z + atom2$z) / 2
    x[1] <- atom1$x
    x[2] <- mx
    y[1] <- atom1$y
    y[2] <- my
    z[1] <- atom1$z
    z[2] <- mz
    lines3d(x, y, z, color=PT.Color[[atom1$el]])
    x[1] <- atom2$x
    y[1] <- atom2$y
    z[1] <- atom2$z
    lines3d(x, y, z, color=PT.Color[[atom2$el]])
  }
}

mol_view_cylindres <- function(mol) {
  x <- double(2)
  y <- double(2)
  z <- double(2)
  nbonds <- length(mol$bonds)
  for (ibond in 1:nbonds) {
    bond <- mol$bonds[[ibond]]
    at1 <- bond$at1
    at2 <- bond$at2
    atom1 <- mol$atoms[[at1]]
    atom2 <- mol$atoms[[at2]]
    mx <- (atom1$x + atom2$x) / 2
    my <- (atom1$y + atom2$y) / 2
    mz <- (atom1$z + atom2$z) / 2
    x[1] <- atom1$x
    x[2] <- mx
    y[1] <- atom1$y
    y[2] <- my
    z[1] <- atom1$z
    z[2] <- mz
    shade3d( cylinder3d(cbind(x,y,z), radius=0.1), color=PT.Color[[atom1$el]])
    x[1] <- atom2$x
    y[1] <- atom2$y
    z[1] <- atom2$z
    shade3d( cylinder3d(cbind(x,y,z), radius=0.1), color=PT.Color[[atom2$el]])
  }
}
