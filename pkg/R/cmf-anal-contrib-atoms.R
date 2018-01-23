# Predictions with analysis

require(rgl)

#source("cinf-mol2.R")
#source("cinf-molview.R")
#source("cinf-ptable.R")
#source("cmf-triposff.R")
#source("cmf-params.R")
#source("cmf-kernels.R")

load("data/PT.RData")
# Making predictions with analysus of ampacts of each atom
cmf_pred_anal_atoms <- function
(
  train_fname = "ligands-train.mol2",  # Training set file name
  pred_fname = "ligands-pred.mol2",    # Prediction set file name
  imol = 1,                            # Molecule from file pred_fname to be analyzed
  model_fname = "ligands-model.RData", # Model file name
  ...
)
{
  mdb0_train <- read_mol2(train_fname)
  mdb0_pred <- read_mol2(pred_fname)
  mdb_train <- cmf_params_tripos(mdb0_train)
  mdb_pred <- cmf_params_tripos(mdb0_pred)
  ntrain <- length(mdb_train)
  
  load(model_fname)
  mfields <- names(model$h)
  nfields <- length(mfields)

  mol1 <- mdb_pred[[imol]]
  natoms1 <- length(mol1$atoms)
  atom_contr <<- array(0.0, c(natoms1, nfields))
  for (iatom1 in 1:natoms1) {
    atom1 <- mol1$atoms[[iatom1]]
	val_all <- 0.0
	for (f in 1:nfields) {
	  field <- mfields[f]
	  val_fld <- 0.0
	  for (imol2 in 1:ntrain) {
	    mol2 <- mdb_train[[imol2]]
		natoms2 <- length(mol2$atoms)
		for (iatom2 in 1:natoms2) {
		  atom2 <- mol2$atoms[[iatom2]]
		  if (field %in% tripos_atom_types) {
		    kern <- cmf_aa_ind_kernel(field, atom1, atom2, model$alpha[[field]])
		  } else {
		    kern <- cmf_aa_kernel(field, atom1, atom2, model$alpha[[field]])
		  }
		  val_fld <- val_fld + model$h[[field]] * model$a[imol2] * kern
		}
	  }
	  val_all <- val_all + val_fld
	  atom_contr[iatom1, f] <<- val_fld
	}
	cat(sprintf("%d\t%s\t%g\n", iatom1, atom1$syb, val_all)); flush.console()
  }
  y_pred <- sum(atom_contr) + model$b
  cat(sprintf("bias=%g\n", model$b))
  cat(sprintf("y_pred=%g\n", y_pred))
  
  # Visualization
  mol_view_cpk(mol1)
  mol_view_cylindres(mol1)
  for (iatom1 in 1:natoms1) {
    atom1 <- mol1$atoms[[iatom1]]
	contr <- sum(atom_contr[iatom1,])
	spheres3d(atom1$x, atom1$y, atom1$z, color = if (contr>0) "red" else "blue", radius=abs(contr), alpha=0.5)
	if (abs(contr) > 0.2) {
	  label <- sprintf("%d", iatom1)
	  shiftx <- 0.3
	  shifty <- 0.3
	  shiftz <- 0.3
	  text3d(atom1$x+shiftx, atom1$y+shifty, atom1$z+shiftz, color="black", text=label)
	}
  }
  
}
