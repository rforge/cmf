# KRR with continuous molecular fields and optimization of  gamma and h.
# alpha is optimized by grid search
# Any combination of molecular fields is supported

#source("cmf-krr.R")
require(conmolfields)

# Activity file name
act_train_fname <- "activity-train.txt"

# Activity column number
act_colnum <- 1

# Separator
sep <- ","

# Kernels file name    
kernels_train_fname <- "ligands-kernels-train.RData"

# Model file name
model_fname <- "ligands-model.RData"

# Types of molecular fields
mfields <- c("q","vdw","logp","abra","abrb")

# Whether to perform grid search for attenuation factor alpha (TRUE/FALSE, 1/0).
alpha_grid_search <- TRUE

# Whether to perform grid search for Tikhonov regularization coefficient gamma (TRUE/FALSE, 1/0).
gamma_grid_search <- FALSE

# Whether to form tconic kernel combination 
conic_kernel_combination <- FALSE

# Whether to optimize h (TRUE/FALSE, 1/0)
optimize_h <- FALSE

# Whether to set b=0 (TRUE/FALSE, 1/0)
set_b_0 <- FALSE

# Print intermediate results in the course of optimiization (TRUE/FALSE, 1/0) 
print_interm_icv <- TRUE                   

# Produce intermediate scatter plots in the course of optimization  (TRUE/FALSE, 1/0)
plot_interm_icv <- TRUE                   

# Print final results on interna; cross-validation (TRUE/FALSE, 1/0)
print_final_icv <- TRUE

# Produce final scatter plot for internal cross-valudation (TRUE/FALSE, 1/0)                     
plot_final_icv <- TRUE                       

cmf_krr_train(
  act_train_fname = act_train_fname, 
  act_colnum = act_colnum, 
  sep = sep, 
  kernels_train_fname = kernels_train_fname, 
  model_fname = model_fname, 
  mfields = mfields, 
  alpha_grid_search = alpha_grid_search,
  gamma_grid_search = gamma_grid_search, 
  conic_kernel_combination = conic_kernel_combination, 
  optimize_h = optimize_h, 
  set_b_0 = set_b_0,
  print_interm_icv = print_interm_icv,
  plot_interm_icv = plot_interm_icv,
  print_final_icv = print_final_icv,
  plot_final_icv = plot_final_icv
)

