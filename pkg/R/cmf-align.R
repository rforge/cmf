# Molecular alignment in 3D space

#source("cinf-mol2.R")
#source("cinf-isomorph.R")
#source("cinf-mol.R")

mol2xyz <- function(mol) {
  natoms <- length(mol$atoms)
  xyz <- matrix(0, 3, natoms)
  for (iatom in 1:natoms) {
    atom <- mol$atoms[[iatom]]
	xyz[1, iatom] <- atom$x
	xyz[2, iatom] <- atom$y
	xyz[3, iatom] <- atom$z
  }
  xyz
}

xyz2mol <- function(mol, xyz) {
  mol1 <- mol
  natoms <- length(mol$atoms)
  for (iatom in 1:natoms) {
    mol1$atoms[[iatom]]$x <- xyz[1, iatom]
    mol1$atoms[[iatom]]$y <- xyz[2, iatom]
    mol1$atoms[[iatom]]$z <- xyz[3, iatom]
  }
  mol1
}

# Generation of orthogonal rotation matrix for Euler angles phi, teta, psi
euler_orth <- function(phi=0, theta=0, psi=0) {
  r <- matrix(0, 3, 3)
  r[1,1] <- cos(theta)*cos(psi)
  r[1,2] <- -cos(phi)*sin(psi) + sin(phi)*sin(theta)*cos(psi)
  r[1,3] <- sin(phi)*sin(psi) + cos(phi)*sin(theta)*cos(psi)
  r[2,1] <- cos(theta)*sin(psi)
  r[2,2] <- cos(phi)*cos(psi) + sin(phi)*sin(theta)*sin(psi)
  r[2,3] <- -sin(phi)*cos(psi) + cos(phi)*sin(theta)*sin(psi)
  r[3,1] <- -sin(theta)
  r[3,2] <- sin(phi)*cos(theta)
  r[3,3] <- cos(phi)*cos(theta)
  r
}

# Generation of rotation matrix around x,y,z axes:
# alpha_(yaw), beta_(pitch), gamma_(roll)
rotmat_xyz <- function(alpha_=0, beta_=0, gamma_=0) {
  rx <- matrix(0, 3, 3)
  rx[1,1] <- 1; rx[1,2] <- 0;           rx[1,3] <- 0;
  rx[2,1] <- 0; rx[2,2] <- cos(gamma_); rx[2,3] <- -sin(gamma_)
  rx[3,1] <- 0; rx[3,2] <- sin(gamma_); rx[3,3] <- cos(gamma_)
  ry <- matrix(0, 3, 3)
  ry[1,1] <- cos(beta_);  ry[1,2] <- 0; ry[1,3] <- sin(beta_)
  ry[2,1] <- 0;           ry[2,2] <- 1; ry[2,3] <- 0
  ry[3,1] <- -sin(beta_); ry[3,2] <- 0; ry[3,3] <- cos(beta_);
  rz <- matrix(0, 3, 3)
  rz[1,1] <- cos(alpha_); rz[1,2] <- -sin(alpha_); rz[1,3] <- 0
  rz[2,1] <- sin(alpha_); rz[2,2] <- cos(alpha_);  rz[2,3] <- 0
  rz[3,1] <- 0;           rz[3,2] <- 0;            rz[3,3] <- 1
  r <- rx %*% ry %*% rz
  r  
}

# Generation of random orthogonal rotation matrix
rnd_euler_orth <- function() {
  phi <- runif(1, -pi, pi)
  theta <- runif(1, -pi, pi)
  psi <- runif(1, -pi, pi)
  euler_orth(phi, theta, psi)
}

# Generation of random rotation matrix
rnd_rotmat_xyz <- function() {
  alpha_ <- runif(1, -pi, pi)
  beta_ <- runif(1, -pi, pi)
  gamma_ <- runif(1, -pi, pi)
  rotmat_xyz(alpha_, beta_, gamma_)
}

# Random translation vector
rnd_trans_vec <- function(coef_=1) {
  (runif(3) - 0.5) * coef_
}

# Perturbate molecule
pert_mol <- function(mol, tcoef=10) {
  xyz <- mol2xyz(mol)
  R <- rnd_euler_orth()
  T <- rnd_trans_vec(tcoef)
  xyz_p <- R %*% xyz + T
  xyz2mol(mol, xyz_p)
}

# Perturbate molecular database
pert_mdb <- function(mdb, tcoef=10) {
  ncomp <- length(mdb)
  mdb_p <- list()
  for (imol in 1:ncomp) {
     mol <- mdb[[imol]]
     mol_p <- pert_mol(mol, tcoef)
	 mdb_p[[imol]] <- mol_p
  }
  mdb_p
}

# Rigin alignment with Arun algorithm
# K.S.Arun, T.S.Huang, S.D.Blostein, "Least Square Fitting of Two 3-D Point Sets",
# IEEE TRANSACTIONS ON PATTERN ANALYSIS AND MACHINE INTELLIGENCE, VOL. PAMI-9, NO. 5, 1987, pp. 698-700
# Two point sets {p_i} and {p1_i}; i=1,2,...,N are related by p1_i = R * p_i + T + Ni, 
# where R is a rotetion matrix, T a translation vector, and N_i a noise vector.
# the algorithm finds the least-squares solution of R and T.
align_arun <- function(p_i, p1_i) {
  # Step 1: From {p_i}, {p1_i} calculate p. p1; and then {q_i}, {q1_i}
  p <- rowMeans(p_i)
  p1 <- rowMeans(p1_i)
  q_i <- p_i - p
  q1_i <- p1_i - p1
  
  # Step 2: Calculate the 3x3 matrix H = <q_i|q1_i>
  H <- q_i %*% t(q1_i)
  
  # Step 3: Find the SVD of H: H=ULV'
  s <- svd(H)
  U <- s$u
  V <- s$v
  
  # Step 4: Calculate X=VU'
  X <- V %*% t(U)
  
  # Step 5: Calculate, det (x), the determinant of X
  detx <- det(X)
  R <- X
  T <- p1 - R %*% p
  
  list(R=R, T=T, detx=detx)
}

transform_xyz <- function(xyz, R, T) {
  n <- dim(xyz)[2]
  ones <- double(n) + 1
  R %*% xyz + T %*% ones
}

transform_mol <- function(mol, R, T) {
  xyz <- mol2xyz(mol)
  xyz1 <- transform_xyz(xyz, R, T)
  xyz2mol(mol, xyz1)
}

# Superposes moving molecule mol_m on templace molecule mol_t
superpose_mol <- function(mol_m, mol_t) {
  xyz_m <- mol2xyz(mol_m)
  xyz_t <- mol2xyz(mol_t)
  align <- align_arun(xyz_m, xyz_t)
  xyz_a <- transform_xyz(xyz_m, align$R, align$T)
  xyz2mol(mol_m, xyz_a)
}

rmse4xyz <- function(xyz1, xyz2) {
  xyz_d2 <- (xyz1 - xyz2)^2
  rss <- sum(colSums(xyz_d2))
  N <- dim(xyz1)[2]
  sqrt(rss / N)
}

rmse4mol <- function(mol1, mol2) {
  xyz1 <- mol2xyz(mol1)
  xyz2 <- mol2xyz(mol2)
  rmse4xyz(xyz1, xyz2)
}

# Aligns molecular database mdb using template templ
align_mdb_template <- function(mdb, templ, iimol=1:length(mdb)) {
  templ_ct <- mol_get_ct(templ)
  templ_lab <- mol_get_chelabs(templ)
  xyz_t <- mol2xyz(templ)
  mdb_a <- list()
  imol1 <- 0
  for (imol in iimol) {
    imol1 <- imol1 + 1
    mol <- mdb[[imol]]
	mol_ct <- mol_get_ct(mol)
	mol_lab <- mol_get_chelabs(mol)
	xyz_m <- mol2xyz(mol)
	isom_list <- find_substr_isomorph(templ_lab, templ_ct, mol_lab, mol_ct)
	nisom <- length(isom_list)	
	substr_mol <- substruct(mol, isom_list[[1]])
	xyz_ss <- mol2xyz(substr_mol)
	align <- align_arun(xyz_ss, xyz_t)
	xyz_ss_a <- transform_xyz(xyz_ss, align$R, align$T)
	rmse_a <- rmse4xyz(xyz_ss_a, xyz_t)
	align_best <- align
	rmse_a_best <- rmse_a
	if (nisom > 1) {
	  for (isom in 2:nisom) {
	    substr_mol <- substruct(mol, isom_list[[isom]])
	    xyz_ss <- mol2xyz(substr_mol)
	    align <- align_arun(xyz_ss, xyz_t)
	    xyz_ss_a <- transform_xyz(xyz_ss, align$R, align$T)
	    rmse_a <- rmse4xyz(xyz_ss_a, xyz_t)
	    if (rmse_a < rmse_a_best) {
		  align_best <- align
		  rmse_a_best <- rmse_a
		}
	  }
	}
    xyz_a <- transform_xyz(xyz_m, align_best$R, align_best$T)
    mol_a <- xyz2mol(mol, xyz_a)
	mdb_a[[imol1]] <- mol_a
	cat(sprintf("imol=%d imol1=%d nisom=%d detx=%g rmse1=%g rmse2=%g\n", 
	  imol, imol1, nisom, align$detx, rmse4mol(substr_mol,templ), rmse_a_best)); flush.console()
  }
  mdb_a
}

#test_align_1 <- function() {
#  mdb <- read_mol2("ligands.mol2")
#  mol <- mdb[[1]]
#  mol_p <- pert_mol(mol)
#  mol_a <- superpose_mol(mol_p, mol)
#}

#test_align_2 <- function() {
#  mdb <- read_mol2("ligands.mol2")
#  mdb_p <- pert_mdb(mdb, tcoef=1)
#  write_mol2(mdb_p, "ligands-perturbed.mol2")
#  mdb_tmp <- read_mol2("ligands-template.mol2")
#  templ <- mdb_tmp[[1]]
#  mdb_a <- align_mdb_template(mdb_p, templ)
#  write_mol2(mdb_a, "ligands-aligned.mol2")
#}



