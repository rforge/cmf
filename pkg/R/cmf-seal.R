# Implementation of the SEAL algorithm

#source("cinf-mol2.R")
#source("cmf-triposff.R")
#source("cmf-align.R")

calc_wij <- function(mol1, mol2, we=1, ws=1) {
  atoms1 <- mol1$atoms
  natoms1 <- length(atoms1)
  atoms2 <- mol2$atoms
  natoms2 <- length(atoms2)
  wij <- matrix(0, natoms1, natoms2)
  for (iatom in 1:natoms1) {
    atom1 <- atoms1[[iatom]]
    q1 <- atom1$pch
    v1 <- tripos_Rvdw[[atom1$syb]]
    for (jatom in 1:natoms2) {
      atom2 <- atoms2[[jatom]]
      q2 <- atom2$pch
      v2 <- tripos_Rvdw[[atom2$syb]]
      wij[iatom,jatom] <- we*q1*q2 + ws*v1*v2 
    }
  }
  wij
}

calc_af <- function(xyz1, xyz2, wij, alpha=0.29) {
  size1 <- dim(xyz1)[2]
  size2 <- dim(xyz2)[2]
  af <- 0
  for (i in 1:size1) {
    for (j in 1:size2) {
      rij2 <- 0
      for (k in 1:3) {
        rij2 <- rij2 + (xyz1[k,i]-xyz2[k,j])^2
      }
      af <- af - wij[i,j] * exp(-alpha * rij2)
    }
  }
  af
}

# Converts quaternions to rotation mstrix
q2r <- function(qw=1, qx=0, qy=0, qz=0) {
  r <- matrix(0, 3, 3)
  nq <- sqrt(qw*qw + qx*qx + qy*qy + qz*qz)
  qw <- qw / nq
  qx <- qx / nq
  qy <- qy / nq
  qz <- qz / nq
  r[1,1] <- 1 - 2*qy*qy - 2*qz*qz
  r[1,2] <- 2*qx*qy - 2*qz*qw
  r[1,3] <- 2*qx*qz + 2*qy*qw
  r[2,1] <- 2*qx*qy + 2*qz*qw
  r[2,2] <- 1 - 2*qx*qx - 2*qz*qz
  r[2,3] <- 2*qy*qz - 2*qx*qw
  r[3,1] <- 2*qx*qz - 2*qy*qw
  r[3,2] <- 2*qy*qz + 2*qx*qw
  r[3,3] <- 1 - 2*qx*qx - 2*qy*qy
  r
}

superpose_mol_seal <- function(mol_m, mol_t, verbose=TRUE, maxit=100) {
  
  calc_mol2mol_transvec <- function(mol_m, mol_r) {
    xyz_m <- mol2xyz(mol_m)
    xyz_r <- mol2xyz(mol_r)
    mean_m <- rowMeans(xyz_m)
    mean_r <- rowMeans(xyz_r)
    mean_r - mean_m
  }
  
  apply_par_list <- function(xyz_m, par_list) {
    qw <- par_list[1]
    qx <- par_list[2]
    qy <- par_list[3]
    qz <- par_list[4]
    tx <- par_list[5]
    ty <- par_list[6]
    tz <- par_list[7]
    R <- q2r(qw, qx, qy, qz)
    Tr <- c(tx, ty, tz)
    xyz_a <- R %*% xyz_m + Tr
    xyz_a
  }
  
  fr <- function(par_list, ...) {
    xyz_a <- apply_par_list(xyz_m, par_list)
    af <- calc_af(xyz_a, xyz_t, wij, ...)
    if (verbose) cat(sprintf("af=%g qw=%g qx=%g qy=%g qz=%g tx=%g ty=%g tz=%g\n", 
                             af, par_list[1], par_list[2], par_list[3], par_list[4], par_list[5], par_list[6], 
                             par_list[7])); flush.console()
    af
  }
  
  wij <- calc_wij(mol_m, mol_t)
  xyz_m <- mol2xyz(mol_m)
  xyz_t <- mol2xyz(mol_t)
  t1 <- calc_mol2mol_transvec(mol_m, mol_t)
  
  res1 <- optim(c(1,0,0,0,t1[1],t1[2],t1[3]), fr, control=list(maxit=maxit),
                method = "BFGS", alpha=0.1)
  res2 <- optim(res1$par, fr, control=list(maxit=maxit),
                method = "BFGS", alpha=0.29)
  
  res_best <- res2
  res_value_best <- res2$value
  
  xyz_a <- apply_par_list(xyz_m, res_best$par)
  mol_a <- xyz2mol(mol_m, xyz_a)
  list(mol=mol_a, af=res_best$value)
}

# Aligns molecular database mdb using template mol_t and algorithm SEAL
align_mdb_seal <- function(mdb, mol_t, verbose=TRUE) {
  ncomp <- length(mdb)
  mdb_a <- list()
  for (imol in 1:ncomp) {
    mol_m <- mdb[[imol]]
    cat(sprintf("imol=%d", imol))
    res <- superpose_mol_seal(mol_m, mol_t, verbose=FALSE)
    cat(sprintf(" af=%g\n", res$af))
    mdb_a[[imol]] <- res$mol
  }
  mdb_a
}

#test_seal_1 <- function() {
#  mdb <- read_mol2("ligands.mol2")
#  mol_m <- mdb[[1]]
#  mol_m <- pert_mol(mol_m)
#  mol_t <- mdb[[1]]
#  superpose_mol_seal(mol_m, mol_t)
#}

#test_seal_2 <- function() {
#  mdb <- read_mol2("ligands.mol2")
#  mdb_p <- pert_mdb(mdb, tcoef=1)
#  write_mol2(mdb_p, "ligands-perturbed.mol2")
#  templ <- mdb[[1]]
#  mdb_a <- align_mdb_seal(mdb_p, templ)
#  write_mol2(mdb_a, "ligands-aligned-seal.mol2")
#}



