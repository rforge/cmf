# Generation of grid for continuous molecular co-fields

#source("cinf-mol2.R")
#source("cmf-grid.R")
#source("cmf-coef.R")
#source("cmf-params.R")
#source("cmf-triposff.R")

# Generation of grid for continuous molecular co-fields
cmf_gen_grid <- function
(
  train_fname = "ligands-train.mol2",       # training set file name
  kernels_fname = "ligands-kernels.RData",  # Computed kernels file name
  model_fname = "ligands-model.RData",      # Model file name
  grid_fname = "ligands-grid-krr.RData",    # Grid with regression coefficients - file name 
  verbose = TRUE,                           # Verbose output
  ...
)
{
  load(kernels_fname)
  load(model_fname)
  
  # Molecular fields
  mfields <- names(model$h)
  nfields <- length(mfields)

  mdb <- read_mol2(train_fname)
  mdb <- cmf_params_tripos(mdb)

  grid <- cmf_init_grid(mdb)

  grids <- list()
  for (f in 1:nfields) {
    field <- mfields[f]
    if (verbose) {
      cat(sprintf("Generating grid for field %s...\n", field))
      flush.console()
    }
    grids[[field]] <- cmf_coef_grid(mdb, model$a, model$alpha[[field]], grid, field)
  }
  
  save(grids, file=grid_fname)
}

# Generation of grid for continuous molecular co-fields
# Extended version with field families
cmf_gen_grid_ex <- function
(
  train_fname = "ligands-train.mol2",       # training set file name
  kernels_fname = "ligands-kernels.RData",  # Computed kernels file name
  model_fname = "ligands-model.RData",      # Model file name
  grid_fname = "ligands-grid-krr.RData",    # Grid with regression coefficients - file name 
  field_family = "PHCH",                    # Field family: 
                                            # ("PHCH" - physico-chemical, 
                                            #  "IND"  - indicator continuous fields,
											#  "MOPAC" - MOPAC fields)
  verbose = TRUE,                           # Verbose output
  ...
)
{
  load(kernels_fname)
  load(model_fname)
  
  # Molecular fields
  mfields <- names(model$h)
  nfields <- length(mfields)

  mdb <- read_mol2(train_fname)
  if (field_family == "PHCH") mdb <- cmf_params_tripos(mdb)

  grid <- cmf_init_grid(mdb)

  grids <- list()
  for (f in 1:nfields) {
    field <- mfields[f]
    if (verbose) {
      cat(sprintf("Generating grid for field %s...\n", field))
      flush.console()
    }
    grids[[field]] <- cmf_coef_grid_ex(mdb, model$a, model$alpha[[field]], grid, field, field_family)
  }
  
  save(grids, file=grid_fname)
}

