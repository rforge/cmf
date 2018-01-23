# Analysis of contributions

#source("cinf-regression.R")

#source("cmf-kernels.R")

# Making predictions with analysis
cmf_pred_anal <- function
(
  model_fname = "ligands-model-pred.RData",          # Model file name
  kernels_pred_fname = "ligands-kernels-pred.RData", # Kernels for prediction file name
  act_colnum = 2,                                    # Column name or activity (if specified)
  sep = ",",                                         # Separator
  act_pred_fname = "activity-pred.txt",              # Activity for prediction set (if specified)
  is_train = FALSE,                                  # Whether this analysis is performed for training set set set
  ...
)
{
  # Number of known property in activity-pred.txt
  iprop <- act_colnum

  load(kernels_pred_fname)
  load(model_fname)
  if (is_train) kernels_pred <- kernels 
  alphas_pred <- kernels_pred$alphas
  if (iprop > 0) {
    act <- read.table(act_pred_fname,header=TRUE,sep=sep)
    y_exp <- act[,iprop]
  } else {
    y_exp <- NA
  }
  
  mfields <- names(model$h)
  nfields <- length(mfields)

  K_pred <- cmf_calc_combined_kernels(kernels_pred, model$h, model$alpha, alphas_pred)  
  
  npred <- dim(K_pred)[1]
  ntrain <- dim(K_pred)[2]
			 
  y_pred <- K_pred %*% model$a + model$b

  if (iprop > 0) {
    regr <- regr_param(y_pred, y_exp)
    cat(sprintf("R2=%g RMSE=%g\n", regr$R2, regr$RMSE))
    flush.console()
    plot(y_pred, y_exp, xlab="Predicted", ylab="Experiment")
    abline(coef=c(0,1))
  }
  
  # Analysis of contributions
  contrib <- array(0.0, c(nfields, npred, ntrain))
  for (f in 1:nfields) {
    fname <- mfields[f]
    kernels_interp <- cmf_kernels_interpolate(kernels_pred[[fname]], model$alpha[[fname]], alphas_pred)
    for (p in 1:npred) {
	  for (t in 1:ntrain) {
	    contrib[f,p,t] <- model$h[[fname]] * model$a[t] * kernels_interp[p,t]
	  }
	}
  }
  
  anal <- list()
  anal$contrib <- contrib
  
  # Field contributions
  anal$fields <- mfields
  anal$fld_contrib_av <- numeric(nfields)
  anal$fld_contrib <- array(0.0, c(npred, nfields))
  for (f in 1:nfields) {
    anal$fld_contrib_av[f] <- sum(contrib[f,,]) / npred
	for (p in 1:npred) {
	  anal$fld_contrib[p,f] <- sum(contrib[f,p,]) 
	}
  }
  
  # Training point contributions
  anal$tp_contrib_av <- numeric(ntrain)
  anal$tp_contrib <- array(0.0, c(npred, ntrain))
  for (t in 1:ntrain) {
    anal$tp_contrib_av[t] <- sum(contrib[,,t]) / npred
	for (p in 1:npred) {
	  anal$tp_contrib[p,t] <- sum(contrib[,p,t]) 
	}
  }  

  anal
}
