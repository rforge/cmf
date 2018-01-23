# Computes CMF kernel matrices and saves them to file

#source("cinf-mol2.R")
#source("cmf-triposff.R")
#source("cmf-params.R")
#source("cmf-kernels.R")
 
alphas <- c(0.01,0.05,0.1,0.25,0.5,0.75,1.0,2.0,3.0,5.0,10.0)

# Computes CMF kernel matrices for training and saves to file
comp_kernels_train <- function
(
  train_fname = "ligands-train.mol2",                    # Training set file name
  kernels_train_fname = "ligands-kernels-train.RData", # The name of the file with kernels for training
  mfields = c("q","vdw","logp","abra","abrb"),         # Molecuar fields
  print_comp_kernels = TRUE,                           # Verbose computation of kernels
  ...
)
{
  mdb0 <- read_mol2(train_fname)
  mdb <- cmf_params_tripos(mdb0)

  nfields <- length(mfields)
  
  syb_types <- get_syb_types_list(mdb)
  
  kernels <- list()
  kernels$alphas <- alphas
  for (f in 1:nfields) {
    kernels[[ mfields[f] ]] <- list()
  }

  for (ialpha in 1:length(alphas)) {
    alpha <- alphas[ialpha]

    for (f in 1:nfields) {
      field <- mfields[f]
      if (print_comp_kernels) {
        cat(sprintf("computing kernel_%s for alpha=%g\n", field, alpha))
        flush.console()
      }
      if (field == "ind") {
        Km <- 0.0
        for(type in syb_types){
          if (print_comp_kernels) cat(type)
          Km <- Km + cmf_indicator_kernel_matrix(mdb, alpha, type, verbose=print_comp_kernels)
        }
      } else {
        Km <- cmf_kernel_matrix_tt(field, mdb, alpha, verbose=print_comp_kernels)
      }
      kernels[[field]][[ialpha]] <- Km
    }
    
  }
  save(kernels, file=kernels_train_fname)
}

# Computes CMF kernel matrices for prediction and saves them to file
comp_kernels_pred <- function
(
  train_fname = "ligands-train.mol2",                # Training set file name
  pred_fname = "ligands-pred.mol2",                  # Prediction set file name
  kernels_pred_fname = "ligands-kernels-pred.RData", # Computed kernels file name
  mfields = c("q","vdw","logp","abra","abrb"),       # Molecular fields
  print_comp_kernels = TRUE,                         # Verbose computation of kerneks
  ...
) 
{
  mdb0_train <- read_mol2(train_fname)
  mdb0_pred <- read_mol2(pred_fname)
  mdb_train <- cmf_params_tripos(mdb0_train)
  mdb_pred <- cmf_params_tripos(mdb0_pred)
  
  nfields <- length(mfields)

  syb_types <- get_syb_types_list(mdb_train)

  kernels_pred <- list()
  kernels_pred$alphas <- alphas
  for (f in 1:nfields) {
    kernels_pred[[ mfields[f] ]] <- list()
  }

  for (ialpha in 1:length(alphas)) {
    alpha <- alphas[ialpha]

    for (f in 1:nfields) {
      field <- mfields[f]
      if (print_comp_kernels) {
        cat(sprintf("computing kernel_%s for alpha=%g\n", field, alpha))
        flush.console()
      }
      if (field == "ind") {
        Km <- 0.0
        for(type in syb_types){
          if (print_comp_kernels) cat(type)
          Km <- Km + cmf_indicator_kernel_matrix_pred(mdb_pred, mdb_train, alpha, type, verbose=print_comp_kernels)
        }
      } else {
        Km <- cmf_kernel_matrix_tp(field, mdb_pred, mdb_train, alpha, verbose=print_comp_kernels)
      }
      kernels_pred[[field]][[ialpha]] <- Km
    }

  }
  save(kernels_pred, file=kernels_pred_fname)
}

# Computes CMF kernel matrices for the combined set of molecules
comp_kernels_all <- function
(
  all_fname = "ligands-all.mol2",                  # The name of the file containing all molecules
  kernels_all_fname = "ligands-kernels-all.RData", # The name of the files containing kernels for all molecules
  ...
)
{
  comp_kernels_train(train_fname=all_fname, kernels_train_fname=kernels_all_fname, ...)
}

# Computes indicator kernel matrices for the training set and saves it to file
comp_ind_kernels_train <- function
(
  train_fname = "ligands-train.mol2",                          # Training set file name
  ind_kernels_train_fname = "ligands-ind-kernels-train.RData", # Computed indicator kernels file name
  print_comp_kernels = TRUE,                                   # Verbose computation of kernels
  ...
)
{
  mdb0 <- read_mol2(train_fname)
  mdb <- cmf_params_tripos(mdb0)

  syb_types <- get_syb_types_list(mdb)
  
  kernels <- list()
  kernels$alphas <- alphas

  for (type in syb_types) {
    kernels[[type]] <- list()
  }

  for (ialpha in 1:length(alphas)) {
    alpha <- alphas[ialpha]

    if (print_comp_kernels) {
      cat(sprintf("computing indicator kernels for alpha=%g\n", alpha))
      flush.console()
    }
    
    for (type in syb_types) {
      cat(type)
      kernels[[type]][[ialpha]] <- cmf_indicator_kernel_matrix(mdb, alpha, type, verbose=print_comp_kernels)
    }
    
  }
  save(kernels, file=ind_kernels_train_fname)
}

# Computes indicator kernel matrices for prediction and saves to file
comp_ind_kernels_pred <- function
(
  train_fname = "ligands-train.mol2",                        # Training set file name
  pred_fname = "ligands-pred.mol2",                          # Prediction set file name
  ind_kernels_pred_fname = "ligands-ind-kernels-pred.RData", # Computed kernels file name
  print_comp_kernels = TRUE,                                 # Verbose computation of kernels
  ...
)
{
  mdb0_train <- read_mol2(train_fname)
  mdb0_pred <- read_mol2(pred_fname)
  mdb_train <- cmf_params_tripos(mdb0_train)
  mdb_pred <- cmf_params_tripos(mdb0_pred)

  syb_types <- get_syb_types_list(mdb_train)
  
  kernels_pred <- list()
  kernels_pred$alphas <- alphas

  for (type in syb_types) {
    kernels_pred[[type]] <- list()
  }

  for (ialpha in 1:length(alphas)) {
    alpha <- alphas[ialpha]

    if (print_comp_kernels) {
      cat(sprintf("computing indicator kernels for alpha=%g\n", alpha))
      flush.console()
    }
    
    for (type in syb_types) {
      if (print_comp_kernels) cat(type)
      kernels_pred[[type]][[ialpha]] <- cmf_indicator_kernel_matrix_pred(mdb_pred, mdb_train, alpha, type, verbose=print_comp_kernels)
    }
    
  }
  save(kernels_pred, file=ind_kernels_pred_fname)
}

# Computes CMF kernel matrices for the combined set of molecules and indicator fields
comp_ind_kernels_all <- function
(
  all_fname = "ligands-all.mol2",                          # The name of the file containing all molecules
  ind_kernels_all_fname = "ligands-ind-kernels-all.RData", # The name of the files containing kernels for all molecules
  ...
)
{
  comp_ind_kernels_train(train_fname=all_fname, ind_kernels_train_fname=ind_kernels_all_fname, ...)
}
