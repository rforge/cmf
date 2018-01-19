# CMF kernels

heavy_atom_fields <- c("logp", "abra", "abrb", "abrs", "abre",
  "mop_dn", "mop_de", "mop_pis", "mop_homo", "mop_lumo")

# Calculation of the square of the Euclidean distance between two atoms
eucldist2 <- function(atom1, atom2) {
  dist2 <- (atom1$x - atom2$x)^2 + (atom1$y - atom2$y)^2 + (atom1$z - atom2$z)^2
}

# Normalization of the Gram matrix
normalize_gram <- function(gram) {
  ndim <- dim(gram)[1]
  for (irow in 1:(ndim-1)) {
    for (icol in (irow+1):ndim) {
      gram[irow,icol] <- gram[irow,icol] / (sqrt(gram[irow,irow]) * sqrt(gram[icol,icol]))
      gram[icol,irow] <- gram[irow,icol]
    }
  }
  for (i in 1:ndim) {
    gram[i,i] <- 1.0
  }
  gram
}

# Computation of the kernel that compares fields of two atoms
cmf_aa_kernel <- function(ft, atom1, atom2, alpha) {
  dist2 <- eucldist2(atom1, atom2)
  val <- 0.0
  if (ft == "q") {
    val <- atom1$pch * atom2$pch * exp(- alpha * dist2 / 4.0)
  } else if (ft == "vdwr") {
    val <- tripos_Rvdw[[atom1$syb]] * tripos_Rvdw[[atom2$syb]] * exp(- alpha * dist2 / 4.0)
  } else if (ft == "vdw") {
    w1 <- tripos_Evdw[[atom1$syb]]
    w2 <- tripos_Evdw[[atom2$syb]]
    a1 <- alpha / tripos_Rvdw[[atom1$syb]]^2
    a2 <- alpha / tripos_Rvdw[[atom2$syb]]^2
    coef <- (4 * sqrt(pi^3)) / ((a1+a2) * sqrt(2*a1+2*a2))
    expart <- exp( - (a1*a2*dist2) / (2 * (a1+a2)) )
    val <- w1 * w2 * coef * expart
  } else if (ft == "logp") {
    val <- atom1$hydroph * atom2$hydroph * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abra") {
    val <- atom1$abraham_a * atom2$abraham_a * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abrb") {
    val <- atom1$abraham_b * atom2$abraham_b * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abrs") {
    val <- atom1$abraham_s * atom2$abraham_s * exp(- alpha * dist2 / 4.0)
  } else if (ft == "abre") {
    val <- atom1$abraham_e * atom2$abraham_e * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_q") {
    val <- atom1$mop_q * atom2$mop_q * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_dn") {
    val <- atom1$mop_dn * atom2$mop_dn * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_de") {
    val <- atom1$mop_de * atom2$mop_de * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_pis") {
    val <- atom1$mop_pis * atom2$mop_pis * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_homo") {
    val <- atom1$mop_homo * atom2$mop_homo * exp(- alpha * dist2 / 4.0)
  } else if (ft == "mop_lumo") {
    val <- atom1$mop_lumo * atom2$mop_lumo * exp(- alpha * dist2 / 4.0)
  } else if (ft == "ind") {
    if (atom1$syb == atom2$syb) {
      val <- exp(- alpha * dist2 / 4.0)
	}
  }
  if (ft != "vdw") {
    val <- val * sqrt(pi^3 / alpha^3)
  }
  val
}

# Computation of the kernel that compares fields of two molecules
# with specified atom lists
cmf_kernel_al <- function(ft, mol1, mol2, alpha, atomlist1, atomlist2) {
  res <- 0.0
  for (iatom1 in atomlist1) {
    atom1 <- mol1$atoms[[iatom1]]
    for (iatom2 in atomlist2) {
      atom2 <- mol2$atoms[[iatom2]]
	  res <- res + cmf_aa_kernel(ft, atom1, atom2, alpha)
    }
  }
  res
}

# Computation of the kernel that compares fields of two molecules
cmf_kernel <- function(ft, mol1, mol2, alpha) {
  natoms1 <- length(mol1$atoms)
  natoms2 <- length(mol2$atoms)
  atomlist1 <- 1:natoms1
  atomlist2 <- 1:natoms2
  cmf_kernel_al(ft, mol1, mol2, alpha, atomlist1, atomlist2)
}

# Creation of atom lists
make_atom_lists <- function(ft, mdb) {
  atomlists <- list()
  ncomp <- length(mdb)
  for (imol in 1:ncomp) {
    mol <- mdb[[imol]]
	natoms <- length(mol$atoms)
	atomlist <- 1:natoms
	if (ft %in% heavy_atom_fields) {
	  isheavy <- logical(natoms)
	  for (i in 1:natoms) isheavy[i] <- mol$atoms[[i]]$el != "H"
      atomlist <- atomlist[isheavy]	  
	}
	atomlists[[imol]] <- atomlist
  }
  atomlists
}

# Computation of the kernel (Gram) matrix for the training set
cmf_kernel_matrix_tt <- function(ft, mdb, alpha, verbose=1) {
  nmol <- length(mdb)
  gram <- matrix(0, nmol, nmol)
  atomlists <- make_atom_lists(ft, mdb)
  for (imol in 1:nmol) {
    if (verbose) {cat("."); flush.console()}
    mol <- mdb[[imol]]
    gram[imol,imol] <- cmf_kernel_al(ft, mol, mol, alpha, atomlists[[imol]], atomlists[[imol]])
  }
  if (verbose) {cat("\n"); flush.console()}
  for (imol1 in 1:(nmol-1)) {
    mol1 <- mdb[[imol1]]
	if (verbose) {cat("."); flush.console()}
    for (imol2 in (imol1+1):nmol) {
      mol2 <- mdb[[imol2]]
      gram[imol1,imol2] <- cmf_kernel_al(ft, mol1, mol2, alpha, atomlists[[imol1]], atomlists[[imol2]])
      gram[imol2,imol1] <- gram[imol1,imol2]
    }
  }
  if (verbose) {cat("\n"); flush.console()}
  gram
}

# Computation of the kernel matrix between the training and prediction sets
cmf_kernel_matrix_tp <- function(ft, mdb1, mdb2, alpha, verbose=1) {
  nmol1 <- length(mdb1)
  nmol2 <- length(mdb2)
  gram <- matrix(0, nmol1, nmol2)
  atomlists1 <- make_atom_lists(ft, mdb1)
  atomlists2 <- make_atom_lists(ft, mdb2)
  for (imol1 in 1:nmol1) {
    mol1 <- mdb1[[imol1]]
    for (imol2 in 1:nmol2) {
	  mol2 <- mdb2[[imol2]]
	  gram[imol1,imol2] <- cmf_kernel_al(ft, mol1, mol2, alpha, atomlists1[[imol1]], atomlists2[[imol2]])
	}
    if (verbose) {cat("."); flush.console()}
  }
  if (verbose) {cat("\n"); flush.console()}
  gram
}

## Indicator fields 

cmf_aa_ind_kernel <- function(field, atom1, atom2, alpha) {
  val <- 0.0
  if (field == atom1$syb && field == atom2$syb) {
    dist2 <- eucldist2(atom1, atom2)
    val <- exp(- alpha * dist2 / 4.0)
    val <- val * sqrt(pi^3 / alpha^3)
  }
  val
}

# Code of Gleb Sitnikov - begin
cmf_indicator_kernel <- function(mol1, mol2, alpha, syb_type) {
  res <- 0.0
  natoms1 <- length(mol1$atoms)
  natoms2 <- length(mol2$atoms)
  for (iatom1 in 1:natoms1) {
    atom1 <- mol1$atoms[[iatom1]]
    if(atom1$syb != syb_type) next
    for (iatom2 in 1:natoms2) {
      atom2 <- mol2$atoms[[iatom2]]
      if(atom2$syb != syb_type) next
      dist2 <- eucldist2(atom1, atom2)
      res <- res + exp(- alpha * dist2 / 4.0)
    }
  }
  coef <- sqrt(pi^3 / alpha^3)
  res <- coef * res
  res
}

cmf_indicator_kernel_matrix <- function(mdb, alpha, syb_type, verbose=1) {
  nmol <- length(mdb)
  gram <- matrix(0, nmol, nmol)
  for (imol in 1:nmol) {
    if (verbose) {cat("."); flush.console()}
    mol <- mdb[[imol]]
    gram[imol,imol] <- cmf_indicator_kernel(mol, mol, alpha, syb_type)
  }
  if (verbose) {cat("\n"); flush.console()}
  for (imol1 in 1:(nmol-1)) {
    mol1 <- mdb[[imol1]]
	if (verbose) {cat("."); flush.console()}
    for (imol2 in (imol1+1):nmol) {
      mol2 <- mdb[[imol2]]
      gram[imol1,imol2] <- cmf_indicator_kernel(mol1, mol2, alpha, syb_type)
      gram[imol2,imol1] <- gram[imol1,imol2]
    }
  }
  if (verbose) {cat("\n"); flush.console()}
  gram
}

cmf_indicator_kernel_matrix_pred <- function(mdb1, mdb2, alpha, syb_type, verbose=1) {
  nmol1 <- length(mdb1)
  nmol2 <- length(mdb2)
  gram <- matrix(0, nmol1, nmol2)
  for (imol1 in 1:nmol1) {
    mol1 <- mdb1[[imol1]]
    for (imol2 in 1:nmol2) {
	  mol2 <- mdb2[[imol2]]
	  gram[imol1,imol2] <- cmf_indicator_kernel(mol1, mol2, alpha, syb_type)
	}
    if (verbose) {cat("."); flush.console()}
  }
  if (verbose) {cat("\n"); flush.console()}
  gram
}
# Code of Gleb Sitnikov - end

## Kernel combination and interpolation 

# Linear interpolation of kernel values
cmf_kernels_interpolate <- function(kernels_a, alpha, alphas)
{
  nalphas <- length(alphas)
  if (alpha >= alphas[1] && alpha <= alphas[nalphas]) {
    iamin <- 1
    iamax <- 2
    for (ia in 1:(nalphas-1)) {
      amin <- alphas[ia]
      amax <- alphas[ia+1]
      if ((alpha>=amin) && (alpha<=amax)) {
        iamin <- ia
        iamax <- ia + 1			
        break
      }
    }
    c1 <- (amax - alpha) / (amax - amin)
    c2 <- (alpha - amin) / (amax - amin)
    res <- c1 * kernels_a[[iamin]] + c2 * kernels_a[[iamax]]
  } else if (alpha < alphas[1]) {
    res <- kernels_a[[1]]
  } else if (alpha > alphas[nalphas]) {
    res <- kernels_a[[nalphas]]
  }  
  return(res)
}

# Kernel combination with linear interpolation
cmf_calc_combined_kernels <- function(kernels, h, alpha_f, alphas) 
{
  mfields <- names(h)
  nfields <- length(mfields)
  for (f in 1:nfields) {
    alpha <- alpha_f[[mfields[f]]]
	if (f == 1) {
      Km <- h[[mfields[f]]] * cmf_kernels_interpolate(kernels[[mfields[f]]], alpha, alphas)
	} else {
      Km <- Km + h[[mfields[f]]] * cmf_kernels_interpolate(kernels[[mfields[f]]], alpha, alphas)
	}
  }
  return(Km)
}

# Kernel combination with linear interpolation for fields with the same value of alpha
cmf_calc_combined_kernels_1alpha <- function(kernels, h, alpha, alphas) 
{
  mfields <- names(h)
  nfields <- length(mfields)
  alpha_f <- list()
  for (f in 1:nfields) {
    alpha_f[[ mfields[f] ]] <- alpha
  }
  cmf_calc_combined_kernels(kernels, h, alpha_f, alphas)
}


