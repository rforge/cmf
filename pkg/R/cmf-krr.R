# KRR with continuous molecular fields and optimization of  gamma and h.
# alpha is optimized by grid search
# This script works with any set of molecular fields

#source("cinf-mol2.R")
#source("cinf-regression.R")
#source("cinf-plots.R")

#source("cmf-triposff.R")
#source("cmf-params.R")
#source("cmf-kernels.R")

gamma_list <- c(0.01,0.05,0.1,0.5,0.6,0.7,0.8,0.9,1.0,2.0,2.5,3.0,3.5,4.0,5.0,
  10.0,15.0,20.0,25.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0,100,150,200,300,400,500)

build_krr_model <- function(K_train, y_train, gamma, set_b_0=FALSE)
{  
  ntrain <- length(y_train)
  if (set_b_0) {
	if (rcond(K_train) < 1.e-15) return(NULL)
    a <- solve(K_train + gamma * diag(ntrain), y_train)
    b <- 0
  } else {
    K1 <- K_train + gamma * diag(ntrain)
    K2 <- cbind(K1, matrix(1, nrow=ntrain))
    K3 <- rbind(K2, matrix(1, ncol=ntrain+1))
    K3[ntrain+1, ntrain+1] <- 0
	if (rcond(K3) < 1.e-15) return(NULL)
    y0 <- c(y_train, 0)
    ab <- solve(K3, y0)
    a <- ab[1:ntrain]
    b <- ab[ntrain+1]
  }
  list(a=a, b=b)
}

# KRR with cross-validation
cv_krr <- function(nfolds, K, y, gamma) {
  y_pred_all <- double()
  y_exp_all <- double()
  new_ind <- double()
  for (ifold in 1:nfolds) {
    ind <- 1:length(y)
    ind_train <- ind[ ind%%nfolds != ifold-1 ]
    ind_test  <- ind[ ind%%nfolds == ifold-1 ]
    K_train <- K[ind_train, ind_train]
    K_test  <- K[ind_test, ind_train]
    y_train <- y[ind_train]
    y_test <- y[ind_test]
	  
	m <- build_krr_model(K_train, y_train, gamma, set_b_0)
	if (is.null(m)) next
	  
    y_exp_all <- c(y_exp_all, y_test)
    y_pred <- K_test %*% m$a + m$b
    y_pred_all <- c(y_pred_all, y_pred)
    new_ind <- c(new_ind, ind_test)
  }
  old_ind <- double(length(y))
  for (i in 1:length(y)) old_ind[new_ind[i]] <- i
  regr <- regr_param(y_pred_all, y_exp_all)
  list(R2=regr$R2, RMSE=regr$RMSE, y_pred_cv=y_pred_all[old_ind])
}

# Building model in memory
cmf_krr_train_mem <- function
(
  y,                                           # Activity values for the training set
  kernels,                                     # Kernels for different values of alphas
  alpha_grid_search = TRUE,                    # Whether to perform grid search for alpha (TRUE/FALSE)
  gamma_grid_search = FALSE,                   # Whether to perform grid search for gamma (TRUE/FALSE)
  conic_kernel_combination = FALSE,            # Whether to form tconic kernel combination (TRUE/FALSE)
  optimize_h = FALSE,                          # Whether to optimize h (TRUE/FALSE)
  mfields = c("q","vdw","logp","abra","abrb"), # Molecular fields
  set_b_0 = FALSE,                             # Whether to set b=0
  print_interm_icv = TRUE,                     # Print intermediate results in the course of optimiization
  plot_interm_icv = TRUE,                      # Produce intermediate scatter plots in the course of optimization 
  print_final_icv = TRUE,                      # Print final results on interna; cross-validation
  plot_final_icv = TRUE,                       # Produce final scatter plot for internal cross-valudation
  ...
)
{
  var_y <- var(y)
  ncomp <- length(y) 

  alphas <- kernels$alphas
  nalphas <- length(alphas)

  nfields <- length(mfields)

  Q2_best_of_best <- -1000

  model <- list()

  # Optimization of hyperparameters
  fr <- function(par_list) {
  
    # Try current set of hyperparameters
	try_current_hyper_params <- function() {
      m <- build_krr_model(Km, y, gamma, set_b_0)
      if (is.null(m)) return()
	  y_pred <- Km %*% m$a + m$b
      regr <- regr_param(y_pred, y)
      RMSE <- regr$RMSE
      R2 <- regr$R2
      cv <- cv_krr(10, Km, y, gamma)
      RMSEcv <- cv$RMSE
      Q2 <- cv$R2
      y_pred_cv <- cv$y_pred_cv
      minQ2R2 <- min(Q2, R2)
      if (minQ2R2 > minQ2R2_best) {
        minQ2R2_best <<- minQ2R2
        RMSE_best <<- RMSE
        R2_best <<- R2
        RMSEcv_best <<- RMSEcv
        Q2_best <<- Q2
		if (alpha_grid_search) {
          ialpha_best <<- ialpha
		}
		alpha_best <<- alpha
        gamma_best <<- gamma
        a_best <<- m$a
        b_best <<- m$b
        y_pred_best <<- y_pred
        y_pred_cv_best <<- y_pred_cv
	  }
    }
	  
    R2_best <- -1000
    RMSE_best <- -1
    Q2_best <- -1000
    RMSEcv_best <- -1
    minQ2R2_best <- -1000
    alpha_best <- -1
    gamma_best <- -1
    a_best <- NULL
    b_best <- NULL
    y_pred_best <- double()
    y_pred_cv_best <- double()
    
    # Unpack parameters
	h <- list()
	pos <- 1
    if (optimize_h) {
      if (conic_kernel_combination) {
        for (f in 1:nfields) h[[mfields[f]]] <- abs(par_list[f]) 
      } else {
        for (f in 1:nfields) h[[mfields[f]]] <- par_list[f] 
      }
	  pos <- pos + nfields
	  if (!alpha_grid_search) {
	    alpha <- par_list[pos]
		pos <- pos + 1
	  }
      if (!gamma_grid_search) gamma  <- par_list[pos]
    } else {
      for (f in 1:nfields) h[[mfields[f]]] <- 1
 	  if (!alpha_grid_search) {
	    alpha <- par_list[pos]
		pos <- pos + 1
	  }
      if (!gamma_grid_search) gamma  <- par_list[pos]
    }  
  
    # Optimization of hyper-parameters
	if (alpha_grid_search) { 
	  # grid search for alpha
	  for (ialpha in 1:length(alphas)) {
        alpha <- alphas[[ialpha]]
	    Km <<- matrix(0, nrow=ncomp, ncol=ncomp)
        for (f in 1:nfields) {
          Km <<- Km + h[[mfields[f]]] * kernels[[mfields[f]]][[ialpha]]
        }
        if (gamma_grid_search) {
          for (gamma in gamma_list) {
            try_current_hyper_params()
          } 
        } else {
	      try_current_hyper_params()
        }
      }
      alpha_best <- alphas[ialpha_best]	  
	} else { 
	  # linear interpolation of Km
	  Km <<- cmf_calc_combined_kernels_1alpha(kernels, h, alpha, alphas) 
      if (gamma_grid_search) {
        for (gamma in gamma_list) {
          try_current_hyper_params()
        } 
      } else {
        try_current_hyper_params()
      }
	}
	
    if (Q2_best > Q2_best_of_best) {
      Q2_best_of_best <<- Q2_best
	  if (print_interm_icv) {
	    for (f in 1:nfields) cat(sprintf("h_%s=%g ", mfields[f], h[[ mfields[f] ]] ))
	    cat(sprintf("\n"))
        cat(sprintf("best: alpha=%g gamma=%g RMSE=%g R2=%g RMSEcv=%g Q2=%g \n",
          alpha_best,gamma_best,RMSE_best,R2_best,RMSEcv_best,Q2_best))
        flush.console()
	  }
	  if (plot_interm_icv) {
        cinf_plotxy(y_pred_cv_best, y, xlab="Predicted", ylab="Experiment",
		  main = "Scatter Plot for Cross-Validation (Internal)")
        abline(coef=c(0,1))
	  }
	  model$gamma  <<- gamma_best
	  for (f in 1:nfields) {
	    model$h[[mfields[f]]] <<- h[[mfields[f]]]
	    model$alpha[[mfields[f]]] <<- alpha_best
		if (alpha_best < alphas[1]) model$alpha[[mfields[f]]] <<- alphas[1]
		if (alpha_best > alphas[nalphas]) model$alpha[[mfields[f]]] <<- alphas[nalphas]
	  }
	  model$R2     <<- R2_best
	  model$RMSE   <<- RMSE_best
	  model$y_pred <<- y_pred_best
	  model$y_exp  <<- y
	  model$Q2     <<- Q2_best
	  model$RMSEcv <<- RMSEcv_best
	  model$y_pred_cv <<- y_pred_cv_best
      model$a <<- a_best
      model$b <<- b_best
    }

    RMSEcv_best
  }

  # Pack parameters
  par_list <- list()
  if (optimize_h) par_list <- c(par_list, rep(1,nfields))
  if (!alpha_grid_search) par_list <- c(par_list, 0.25)	
  if (!gamma_grid_search) par_list <- c(par_list, 5)
  npars <- length(par_list)
  if (npars > 1) {
    res <- optim(par_list, fr)
  } else if (npars == 1) {
    res <- optimize(fr, c(0.01, 20))
  } else {
    res <- fr()
  }

  model$set_b_0 <- set_b_0
  
  if (print_final_icv) {
    for (f in 1:nfields) cat(sprintf("h_%s=%g ", mfields[f], model$h[[ mfields[f] ]] ))
    cat(sprintf("\n"))
    cat(sprintf("final: alpha=%g gamma=%g RMSE=%g R2=%g RMSEcv=%g Q2=%g \n",
      model$alpha[1], model$gamma, model$RMSE, model$R2, model$RMSEcv, model$Q2))
    flush.console()
  }
  if (plot_final_icv) {
    cinf_plotxy(model$y_pred_cv, y, xlab="Predicted", ylab="Experiment",
	  main="Scatter Plot for Cross-Validation (Internal)")
    abline(coef=c(0,1))
  }
  
  model
}

# Building model 
cmf_krr_train <- function
(
  act_train_fname = "activity-train.txt",              # Activity file name
  act_colnum = 2,                                      # Activity column number
  sep = ",",                                           # Separator
  kernels_train_fname = "ligands-kernels-train.RData", # Kernels file name
  model_fname = "ligands-model.RData",                 # Model file name
  ...
)
{
  act <- read.table(act_train_fname, header=TRUE, sep=sep)
  y <- act[,act_colnum]
  load(kernels_train_fname)

  model <- cmf_krr_train_mem(y=y, kernels=kernels, ...)
	
  save(model, file=model_fname)
}

# Making predictions in memory
cmf_krr_pred_mem <- function
(
  model,             # Model
  kernels_pred,      # Kernels for prediction
  y_exp,             # Experimental activity values if it is needed to compare with predictionq
  print_pred = TRUE, # Print statistical parameters for prediction
  plot_pred = TRUE,  # Produce scatter plot prediction on external database
  ...
)
{
#  alphas_train <- kernels$alphas
  alphas_pred <- kernels_pred$alphas
  y_train <- model$y_exp

  K_pred <- cmf_calc_combined_kernels(kernels_pred, model$h, model$alpha, alphas_pred)  

  y_pred <- K_pred %*% model$a + model$b

  if (!is.na(y_exp[1])) {
    if (print_pred) {
      regr <- regr_param(y_pred, y_exp)
	  r2ex <- regr_param_ex(y_pred, y_exp, model$y_exp)
      cat(sprintf("R2pred=%g RMSEpred=%g (%g%%) R2pred_ex=%g\n", regr$R2, regr$RMSE, regr$RMSE_pc, r2ex))
      flush.console()
	}
	if (plot_pred) {
      cinf_plotxy(y_pred, y_exp, xlab="Predicted", ylab="Experiment",
	    main = "Scatter Plot for External Prediction")
      abline(coef=c(0,1))
	}
  }

  y_pred
}

# Making predictions 
cmf_krr_pred <- function
(
  model_fname = "ligands-model.RData",                 # Model file name
  kernels_train_fname = "ligands-kernels-train.RData", # Kernels for training file name
  kernels_pred_fname = "ligands-kernels-pred.RData",   # Kernels for prediction file name
  act_colnum = 2,                                      # Column name or activity (if specified)
  sep = ",",                                           # Separator
  act_pred_fname = "activity-pred.txt",                # Activity for prediction set (if specified)
  pred_fname = "ligands-pred.RData",                   # File with predictions
  ...
)
{
  load(model_fname)
  
  load(kernels_train_fname)
  load(kernels_pred_fname)  
  
  iprop <- act_colnum
  if (iprop > 0) {
    act <- read.table(act_pred_fname,header=TRUE,sep=sep)
    y_exp <- act[,iprop]
  } else {
    y_exp <- NA
  }

  y_pred <- cmf_krr_pred_mem(model=model, kernels=kernels, kernels_pred=kernels_pred, y_exp=y_exp, ...)
  
  pred <- list(y_pred=y_pred, y_exp=y_exp)
  save(pred, file=pred_fname);
  pred
}

cmf_extract_subkernels <- function(kernels, ind_rows, ind_cols, mfields)
{
  alphas <- kernels$alphas
  nfields <- length(mfields)
  subkernels <- list()
  subkernels$alphas <- alphas
  for (f in 1:nfields) {
    field <- mfields[f]
    subkernels[[field]] <- list()
    for (ialpha in 1:length(alphas)) {
      subkernels[[field]][[ialpha]] <- kernels[[field]][[ialpha]][ind_rows, ind_cols]
    }
  }
  subkernels
}

# External n-fold cross-validation in memory
cmf_krr_ecv_mem <- function
(
  nfolds = 5,                                  # The number of folds
  y,                                           # Activity values for the training set
  kernels,                                     # Kernels for different values of alphas
  mfields = c("q","vdw","logp","abra","abrb"), # Molecular fields
  print_ecv = TRUE,                            # Produce text output
  plot_ecv = TRUE,                             # Produce scatter plot
...  
)
{
  y_pred_all <- double()
  y_exp_all <- double()
  new_ind <- double()
  models <- list()
  indexes <- list()
  for (ifold in 1:nfolds) {
    if (print_ecv) {
      cat(sprintf("fold = %d of %d\n", ifold, nfolds))
      flush.console()
    }
    ind <- 1:length(y)
    ind_train <- ind[ ind%%nfolds != ifold-1 ]
    ind_test  <- ind[ ind%%nfolds == ifold-1 ]
    kernels_train <- cmf_extract_subkernels(kernels, ind_train, ind_train, mfields)
    kernels_test <- cmf_extract_subkernels(kernels, ind_test, ind_train, mfields)
    y_train <- y[ind_train]
    y_test <- y[ind_test]
    y_exp_all <- c(y_exp_all, y_test)
	  
    model <- cmf_krr_train_mem(
	  y = y_train, 
	  kernels = kernels_train, 
	  mfields = mfields,
	  ...
	)
    y_pred <- cmf_krr_pred_mem(
	  model = model, 
	  kernels = kernels_train, 
	  kernels_pred = kernels_test, 
	  y_exp = y_test,
	  ...
	)
	  
    y_pred_all <- c(y_pred_all, y_pred)
    new_ind <- c(new_ind, ind_test)
	
	models[[ifold]] <- model
    indexes[[ifold]] <- ind_train	
  }
  old_ind <- double(length(y))
  for (i in 1:length(y)) old_ind[new_ind[i]] <- i
  regr <- regr_param(y_pred_all[old_ind], y_exp_all[old_ind])
  Q2ecv <- regr$R2
  RMSEecv <- regr$RMSE
  RMSEecv_pc <- regr$RMSE_pc
	
  if (print_ecv) {
    cat(sprintf("Q2ecv=%.6f RMSEecv=%.6f (%g%%)\n", Q2ecv, RMSEecv, RMSEecv_pc))
    flush.console()
  }
 
  if (plot_ecv) {
    cinf_plotxy(y_pred_all[old_ind], y_exp_all[old_ind], xlab="Predicted", ylab="Experiment",
	  main = "Scatter Plot for External Cross-Validation")
    abline(coef=c(0,1))
  }
  
  list(Q2ecv=Q2ecv, RMSEecv=RMSEecv, RMSEecv_pc=RMSEecv_pc, 
    y_pred_ecv=y_pred_all[old_ind], y_exp=y_exp_all[old_ind], 
	models=models, indexes=indexes)
}

# External n-fold cross-validation
cmf_krr_ecv <- function
(
  nfolds = 5,                                      # The number of folds
  act_all_fname = "activity-all.txt",              # Activity file name
  act_colnum = 2,                                  # Activity column number
  sep = ",",                                       # Separator
  kernels_all_fname = "ligands-kernels-all.RData", # Kernels file name
  ecv_fname = "ligands-ecv.RData",                 # File with external cross-validation results
  ...
)
{
  act <- read.table(act_all_fname, header=TRUE, sep=sep)
  y <- act[,act_colnum]
  load(kernels_all_fname)

  ecv <- cmf_krr_ecv_mem(nfolds=nfolds, y=y, kernels=kernels, ...)
  
  save(ecv, file=ecv_fname)
  
  ecv
}

# Permute kernel matrix using given permutation
cmf_permute_kernels <- function(kernels, permutation, mfields, for_pred=FALSE)
{
  alphas <- kernels$alphas
  nfields <- length(mfields)
  permuted_kernels <- list()
  permuted_kernels$alphas <- alphas
  for (f in 1:nfields) {
    field <- mfields[f]
    permuted_kernels[[field]] <- list()
    for (ialpha in 1:length(alphas)) {
	  if (for_pred)
        permuted_kernels[[field]][[ialpha]] <- kernels[[field]][[ialpha]][, permutation]
	  else
        permuted_kernels[[field]][[ialpha]] <- kernels[[field]][[ialpha]][permutation, permutation]
    }
  }
  permuted_kernels
}

# External n-fold cross-validation with reshuffling in memory
cmf_krr_ecvr_mem <- function
(
  nreshuffles = 10,                            # The number of reshuffles
  y,                                           # Activity values for the training set
  kernels,                                     # Kernels for different values of alphas
  mfields = c("q","vdw","logp","abra","abrb"), # Molecular fields
  print_ecvr = TRUE,                           # Print statistical parameters
  plot_ecvr = TRUE,                            # Produce scatter plot
  seed = -1,                                   # Seed for ranfom number generator
  ...
)
{
  if (seed >= 0) set.seed(seed)
  y_init <- y
  kernels_init <- kernels
  ncomp <- length(y)
  oldnum <- integer(10)
  YPred <- array(dim=c(ncomp, nreshuffles))
  Q2ecv_array <- double(nreshuffles)
  RMSEecv_array <- double(nreshuffles)
  permutations <- list()
  models <- list()
  indexes <- list()
  for (p in 1:nreshuffles) {
    if (print_ecvr) {
      cat(sprintf("reshuffle %d of %d\n", p, nreshuffles))
      flush.console()
	}
    permutation <- sample(1:ncomp, ncomp)
	for (i in 1:ncomp) oldnum[permutation[i]] <- i
	y_perm <- y_init[permutation]
	kernels_perm <- cmf_permute_kernels(kernels_init, permutation, mfields)
	res_perm <- cmf_krr_ecv_mem(y=y_perm, kernels=kernels_perm, mfields=mfields, ...)
	YPred[,p] <- res_perm$y_pred_ecv[oldnum]
	Q2ecv_array[p] <- res_perm$Q2ecv
    RMSEecv_array[p] <- res_perm$RMSEecv	
	models <- c(models, res_perm$models)
	indexes <- c(indexes, res_perm$indexes)
	last <- length(permutations)
	for (ifold in 1:nfolds) permutations[[last+ifold]] <- permutation
  }
  y_pred_mean <- rowMeans(YPred)
  y_pred_sd <- double(ncomp)
  for (i in 1:ncomp) y_pred_sd[i] <- sd(YPred[i,])
  Q2ecv_mean <- mean(Q2ecv_array)
  Q2ecv_sd <- sd(Q2ecv_array)
  RMSEecv_mean <- mean(RMSEecv_array)
  RMSEecv_sd <- sd(RMSEecv_array)
  regr <- regr_param(y_pred_mean, y)
  Q2ecv_aggr <- regr$R2
  RMSEecv_aggr <- regr$RMSE
  RMSEecv_aggr_pc <- regr$RMSE_pc
  nmodels <- length(models)
 
  if (print_ecvr) {
    cat(sprintf("Q2ecv_aggr=%.6f RMSEecv_aggr=%.6f (%g%%)\n", 
	  Q2ecv_aggr, RMSEecv_aggr, RMSEecv_aggr_pc))
    cat(sprintf("Q2ecv_mean=%.6f Q2ecv_sd=%.6f RMSEecv_mean=%.6f RMSEecv_sd=%.6f\n", 
	  Q2ecv_mean, Q2ecv_sd, RMSEecv_mean, RMSEecv_sd))
    flush.console()
  }
 
  if (plot_ecvr) {
    cinf_plotxy(y_pred_mean, y, xlab="Predicted", ylab="Experiment",
	  main = "Scatter Plot for External Cross-Validations with Reshuffles")
    abline(coef=c(0,1))
  }
  
  list(Q2ecv_aggr=Q2ecv_aggr, RMSEecv_aggr=RMSEecv_aggr, RMSEecv_aggr_pc=RMSEecv_aggr_pc,
    YPred=YPred, y_exp=y, y_pred_mean=y_pred_mean, y_pred_sd=y_pred_sd,
	Q2ecv_array=Q2ecv_array, Q2ecv_mean=Q2ecv_mean, Q2ecv_sd=Q2ecv_sd,
	RMSEecv_array=RMSEecv_array, RMSEecv_mean=RMSEecv_mean, RMSEecv_sd=RMSEecv_sd,
	nmodels=nmodels, permutations=permutations, models=models, indexes=indexes, mfields=mfields)  
}

# External n-fold cross-validation with reshuffings
cmf_krr_ecvr <- function
(
  act_all_fname = "activity-all.txt",              # Activity file name
  act_colnum = 2,                                  # Activity column number
  sep = ",",                                       # Separator
  kernels_all_fname = "ligands-kernels-all.RData", # Kernels file name
  ecvr_fname = "ligands-ecvr.RData",               # File name for reshuffled external cross-validation
  ...
)
{
  act <- read.table(act_all_fname, header=TRUE, sep=sep)
  y <- act[,act_colnum]
  load(kernels_all_fname)

  ecvr <- cmf_krr_ecvr_mem(y=y, kernels=kernels, ...)

  save(ecvr, file=ecvr_fname)
  
  ecvr
}

# Making predictions using ecvr results in memory
cmf_krr_ecvr_pred_mem <- function
(
  ecvr,              # Results for reshuffled external cross-validation with models
  kernels,           # Kernels for the training set
  kernels_pred,      # Kernels for prediction
  y_exp,             # Experimental activity values if it is needed to compare with predictionq
  print_pred = TRUE, # Print statistical parameters for prediction
  plot_pred = TRUE,  # Produce scatter plot prediction on external database
  ...
)
{
  ntrain <- dim(kernels_pred[[2]][[1]])[2]
  npred  <- dim(kernels_pred[[2]][[1]])[1]
  nmodels <- ecvr$nmodels
  models <- ecvr$models
  permutations <- ecvr$permutations
  indexes <- ecvr$indexes
  mfields <- ecvr$mfields
  YPred <- matrix(0, npred, nmodels)
  for (imod in 1:nmodels) {
    model <- models[[imod]]
	permutation <- permutations[[imod]]
	ind_train <- indexes[[imod]]
	kernels_perm <- cmf_permute_kernels(kernels, permutation, mfields)
	kernels_pred_perm <- cmf_permute_kernels(kernels_pred, permutation, mfields, for_pred=TRUE)
    kernels_perm_ind <- cmf_extract_subkernels(kernels_perm, ind_train, ind_train, mfields)
    kernels_pred_perm_ind <- cmf_extract_subkernels(kernels_pred_perm, 1:npred, ind_train, mfields)
	y_pred <- cmf_krr_pred_mem(
	  model=model, kernels=kernels_perm_ind, kernels_pred=kernels_pred_perm_ind, 
      y_exp=y_exp, print_pred = FALSE, plot_pred = FALSE, ...
	)
	YPred[,imod] <- y_pred
  }
  YPredMean <- rowMeans(YPred)
  YPredSD <- double(npred)
  for (i in 1:npred) YPredSD[i] <- sd(YPred[i,])
  
  if (!is.na(y_exp[1])) {
    if (print_pred) {
	  dif <- y_exp - YPredMean
	  cat("No.   Prediction   Experiment  Difference\n")
	  for (ipred in 1:npred) {
	    cat(sprintf("%3d  %.3f +- %.3f   %.3f    %6.3f\n", ipred, YPredMean[ipred], YPredSD[ipred]*2, y_exp[ipred], dif[ipred]))
	  }
      if (npred > 1) {
	    regr <- regr_param(YPredMean, y_exp)
        cat(sprintf("R2pred=%g RMSEpred=%g (%g%%)\n", regr$R2, regr$RMSE, regr$RMSE_pc))
	  }
      flush.console()
	}
	if (plot_pred) {
      cinf_plotxy(YPredMean, y_exp, xlab="Predicted", ylab="Experiment",
	    main = "Scatter Plot for Prediction")
      abline(coef=c(0,1))
	}
  }
  
  if (plot_pred) {
    cinf_plotxy(YPredMean, y_exp, xlab="Predicted", ylab="Experiment",
	  main = "Scatter Plot for Prediction")
    abline(coef=c(0,1))
  }
  
  list(YPred=YPred, YPredMean=YPredMean, YPredSD=YPredSD, y_exp=y_exp)
}

# Making predictions using ecvr results 
cmf_krr_ecvr_pred <- function
(
  ecvr_fname = "ligands-ecvr.RData",                   # File name for reshuffled external cross-validation
  kernels_train_fname = "ligands-kernels-train.RData", # Kernels for training file name
  kernels_pred_fname = "ligands-kernels-pred.RData",   # Kernels for prediction file name
  act_colnum = 2,                                      # Column name or activity (if specified)
  sep = ",",                                           # Separator
  act_pred_fname = "activity-pred.txt",                # Activity for prediction set (if specified)
  pred_fname = "ligands-pred.RData",                   # File with predictions
  ...
)
{
  load(ecvr_fname)
  
  load(kernels_train_fname)
  load(kernels_pred_fname)  
  
  iprop <- act_colnum
  if (iprop > 0) {
    act <- read.table(act_pred_fname,header=TRUE,sep=sep)
    y_exp <- act[,iprop]
  } else {
    y_exp <- NA
  }

  ecvr_pred <- cmf_krr_ecvr_pred_mem(ecvr=ecvr, kernels=kernels, kernels_pred=kernels_pred, y_exp=y_exp, ...)
  
  save(ecvr_pred, file=pred_fname);
}
