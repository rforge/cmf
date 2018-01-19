# Computes molecular fields 

#source("cmf-triposff.R")

# Computes field value at point (x,y,z)
cmf_fval_xyz <- function(ft, mol, alpha, x, y, z) {
  fval <- 0.0
  natoms <- length(mol$atoms)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
    dist2 <- (x - atom$x)^2 + (y - atom$y)^2 + (z - atom$z)^2
	if (ft == "q") {
      fval <- fval + atom$pch * exp(- alpha * dist2 / 2.0)
	} else if (ft == "vdw") {
      dist2rel <- dist2 / (tripos_Rvdw[[atom$syb]])^2
      fval <- fval + tripos_Evdw[[atom$syb]] * exp(- alpha * dist2rel / 2.0)
	} else if (ft == "logp") {
      fval <- fval + atom$hydroph * exp(- alpha * dist2 / 2.0)
	} else if (ft == "abra") {
      fval <- fval + atom$abraham_a * exp(- alpha * dist2 / 2.0)
	} else if (ft == "abrb") {
      fval <- fval + atom$abraham_b * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_q") {
      fval <- fval + atom$mop_q * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_dn") {
      fval <- fval + atom$mop_dn * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_de") {
      fval <- fval + atom$mop_de * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_pis") {
      fval <- fval + atom$mop_pis * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_homo") {
      fval <- fval + atom$mop_homo * exp(- alpha * dist2 / 2.0)
    } else if (ft == "mop_lumo") {
      fval <- fval + atom$mop_lumo * exp(- alpha * dist2 / 2.0)
	} else if (ft == "ind") {
      fval <- fval + exp(- alpha * dist2 / 2.0)
	} else if (ft %in% tripos_atom_types) {
	  if (ft == atom$syb) {
        fval <- fval + exp(- alpha * dist2 / 2.0)
	  }
	}
  }
  fval
}

# Computes field values for grid
cmf_fval_grid <- function(ft, mol, alpha, grid, verbose=1) {
  for (igridx in 1:grid$ngridx) {
    if (verbose) {cat("."); flush.console()}
    x <- grid$gridx[igridx]
    for (igridy in 1:grid$ngridy) {
      y <- grid$gridy[igridy]
      for (igridz in 1:grid$ngridz) {
        z <- grid$gridz[igridz]
        grid$val[igridx,igridy,igridz] <- cmf_fval_xyz(ft, mol,alpha,x,y,z)
      }
    }
  }
  if (verbose) {cat("\n"); flush.console()}
  grid
}

# Get product of two fields
cmf_multiply_fields <- function(grid1, grid2) {
  grid <- list()
  grid$ngridx <- grid1$ngridx
  grid$ngridy <- grid1$ngridy
  grid$ngridz <- grid1$ngridz
  grid$gridx <- grid1$gridx
  grid$gridy <- grid1$gridy
  grid$gridz <- grid1$gridz
  grid$val <- grid1$val * grid2$val
  grid
}


