# Computes molecular field coefficients

#source("cmf-atomlists.R")

# Computes coefficient  at point (x,y,z) using support vectors
cmf_coef_xyz_sv <- function(mdb, a, ai, alpha, x, y, z, field) {
  coef <- 0.0
  nsv <- length(ai)
  for (isv in 1:nsv) {
    imol <- ai[isv]
	mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
   for (iatom in 1:natoms) {
      atom <- mol$atoms[[iatom]]
      dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
	  if (field == "q") {
        coef <- coef + a[isv] * atom$pch * exp(- alpha * dist2 / 2.0)
	  } else if (field == "vdw") {
        dist2rel <- dist2 / (tripos_Rvdw[[atom$syb]])^2
        coef <- coef + a[isv] * tripos_Evdw[[atom$syb]] * exp(- alpha * dist2rel / 2.0)
	  } else if (field == "logp") {
        coef <- coef + a[isv] * atom$hydroph * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abra") {
        coef <- coef + a[isv] * atom$abraham_a * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abrb") {
        coef <- coef + a[isv] * atom$abraham_b * exp(- alpha * dist2 / 2.0)
	  } else if (field == "vdwr") {
        coef <- coef + a[isv] * tripos_Rvdw[[atom$syb]] * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abrs") {
        coef <- coef + a[isv] * atom$abraham_s * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abre") {
        coef <- coef + a[isv] * atom$abraham_e * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_q") {
        coef <- coef + a[isv] * atom$mop_q * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_dn") {
        coef <- coef + a[isv] * atom$mop_dn * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_de") {
        coef <- coef + a[isv] * atom$mop_de * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_pis") {
        coef <- coef + a[isv] * atom$mop_pis * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_homo") {
        coef <- coef + a[isv] * atom$mop_homo * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_lumo") {
        coef <- coef + a[isv] * atom$mop_lumo * exp(- alpha * dist2 / 2.0)
	  } else if (field == "ind") {
        coef <- coef + a[isv] * exp(- alpha * dist2 / 2.0)
	  } else if (field %in% tripos_atom_types) {
	    if (atom$syb == field) {
          coef <- coef + a[isv] * exp(- alpha * dist2 / 2.0)
		}
	  }
    }
  }
  coef
}

# Computes coefficient  at point (x,y,z)
cmf_coef_xyz <- function(mdb, a, alpha, x, y, z, field) {
  nmols <- length(a)
  ai <- 1:nmols
  cmf_coef_xyz_cv(mdb, a, ai, alpha, x, y, z, field)
}

# Computes coefficients for grid using support vectors
cmf_coef_grid_sv <- function(mdb, a, ai, alpha, grid, field) {
  for (igridx in 1:grid$ngridx) {
    cat(sprintf("x: %d of %d\n", igridx, grid$ngridx)); flush.console()
    x <- grid$gridx[igridx]
    for (igridy in 1:grid$ngridy) {
      y <- grid$gridy[igridy]
      for (igridz in 1:grid$ngridz) {
        z <- grid$gridz[igridz]
        grid$val[igridx,igridy,igridz] <- cmf_coef_xyz_sv(mdb,a,ai,alpha,x,y,z,field)
      }
    }
  }
  grid
}

# Computes coefficients for grid
cmf_coef_grid <- function(mdb, a, alpha, grid, field) {
  nmols <- length(a)
  ai <- 1:nmols
  cmf_coef_grid_sv(mdb, a, ai, alpha, grid, field)
}

### Extended versions with field families ###

# Computes coefficient  at point (x,y,z) using support vectors
# for physico-chemical molecular fields
cmf_coef_xyz_sv_ex_phch <- function(mdb, a, ai, alpha, x, y, z, field, atomlists) {
  coef <- 0.0
  nsv <- length(ai)
  for (isv in 1:nsv) {
    imol <- ai[isv]
	mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
    for (iatom in atomlists[[imol]]) {
      atom <- mol$atoms[[iatom]]
      dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
	  if (field == "q") {
        coef <- coef + a[isv] * atom$pch * exp(- alpha * dist2 / 2.0)
	  } else if (field == "vdw") {
        dist2rel <- dist2 / (tripos_Rvdw[[atom$syb]])^2
        coef <- coef + a[isv] * tripos_Evdw[[atom$syb]] * exp(- alpha * dist2rel / 2.0)
	  } else if (field == "logp") {
        coef <- coef + a[isv] * atom$hydroph * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abra") {
        coef <- coef + a[isv] * atom$abraham_a * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abrb") {
        coef <- coef + a[isv] * atom$abraham_b * exp(- alpha * dist2 / 2.0)
	  } else if (field == "vdwr") {
        coef <- coef + a[isv] * tripos_Rvdw[[atom$syb]] * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abrs") {
        coef <- coef + a[isv] * atom$abraham_s * exp(- alpha * dist2 / 2.0)
	  } else if (field == "abre") {
        coef <- coef + a[isv] * atom$abraham_e * exp(- alpha * dist2 / 2.0)
	  }
    }
  }
  coef
}

# Computes coefficient  at point (x,y,z) using support vectors
# for continuous indicator field
cmf_coef_xyz_sv_ex_ind <- function(mdb, a, ai, alpha, x, y, z, field, atomlists) {
  coef <- 0.0
  nsv <- length(ai)
  for (isv in 1:nsv) {
    imol <- ai[isv]
	mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
    for (iatom in atomlists[[imol]]) {
      atom <- mol$atoms[[iatom]]
      dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
      coef <- coef + a[isv] * exp(- alpha * dist2 / 2.0)
    }
  }
  coef
}

# Computes coefficient  at point (x,y,z) using support vectors
# for mopac molecular fields
cmf_coef_xyz_sv_ex_mopac <- function(mdb, a, ai, alpha, x, y, z, field, atomlists) {
  coef <- 0.0
  nsv <- length(ai)
  for (isv in 1:nsv) {
    imol <- ai[isv]
	mol <- mdb[[imol]]
    natoms <- length(mol$atoms)
    for (iatom in atomlists[[imol]]) {
      atom <- mol$atoms[[iatom]]
      dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
	  if (field == "mop_q") {
        coef <- coef + a[isv] * atom$mop_q * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_dn") {
        coef <- coef + a[isv] * atom$mop_dn * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_de") {
        coef <- coef + a[isv] * atom$mop_de * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_pis") {
        coef <- coef + a[isv] * atom$mop_pis * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_homo") {
        coef <- coef + a[isv] * atom$mop_homo * exp(- alpha * dist2 / 2.0)
	  } else if (field == "mop_lumo") {
        coef <- coef + a[isv] * atom$mop_lumo * exp(- alpha * dist2 / 2.0)
	  }
    }
  }
  coef
}

# Computes coefficient  at point (x,y,z)
cmf_coef_xyz_ex <- function(mdb, a, alpha, x, y, z, field, field_family="PHCH") {
  atomlists <- gen_atomlists(mdb, field, field_family)
  nmols <- length(a)
  ai <- 1:nmols
  if (field_family == "PHCH") {
    cmf_coef_xyz_sv_ex_phch(mdb, a, ai, alpha, x, y, z, field, atomlists)
  } else if (field_family == "IND") {
    cmf_coef_xyz_sv_ex_ind(mdb, a, ai, alpha, x, y, z, field, atomlists)
  } else if (field_family == "MOPAC") {
    cmf_coef_xyz_sv_ex_mopac(mdb, a, ai, alpha, x, y, z, field, atomlists)
  }
}

# Computes coefficients for grid using support vectors
cmf_coef_grid_sv_ex <- function(mdb, a, ai, alpha, grid, field, field_family="PHCH") {
  atomlists <- gen_atomlists(mdb, field, field_family)
  for (igridx in 1:grid$ngridx) {
    cat(sprintf("x: %d of %d\n", igridx, grid$ngridx)); flush.console()
    x <- grid$gridx[igridx]
    for (igridy in 1:grid$ngridy) {
      y <- grid$gridy[igridy]
      for (igridz in 1:grid$ngridz) {
        z <- grid$gridz[igridz]
		if (field_family == "PHCH") {
          grid$val[igridx,igridy,igridz] <- cmf_coef_xyz_sv_ex_phch(mdb, a, ai, alpha, x, y, z, field, atomlists)
		} else if (field_family == "IND") {
          grid$val[igridx,igridy,igridz] <- cmf_coef_xyz_sv_ex_ind(mdb, a, ai, alpha, x, y, z, field, atomlists)
		} else if (field_family == "MOPAC") {
          grid$val[igridx,igridy,igridz] <- cmf_coef_xyz_sv_ex_mopac(mdb, a, ai, alpha, x, y, z, field, atomlists)
		}
      }
    }
  }
  grid
}

# Computes coefficients for grid
cmf_coef_grid_ex <- function(mdb, a, alpha, grid, field, field_family="PHCH") {
  nmols <- length(a)
  ai <- 1:nmols
  cmf_coef_grid_sv_ex(mdb, a, ai, alpha, grid, field, field_family=field_family)
}


