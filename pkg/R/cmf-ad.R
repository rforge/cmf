# Applicability domain studies

# Produce coverage - mean error plot for applicability domain studies
cinf_plot_coverage_mean_error <- function
(
  y_err,             # Prediction errors 
  d2m,               # Distance to model
  cmep_type = "in",  # Coverage - mean error plot type ("in", "out", "diff")
  smoothing = TRUE,  # Whether to draw smoothing curve
  ...
) 
{
  npoints <- length(y_err)
  
  oldnum <- order(d2m)
  y_err <- y_err[oldnum]
  d2m <- d2m[oldnum]

  y_err_ave_in <- y_err
  for (i in 2:npoints) y_err_ave_in[i] <- mean(y_err[1:i])
  y_err_ave_out <- y_err
  for (i in 1:(npoints-1)) y_err_ave_out[i] <- mean(y_err[i:npoints])

  coverage <- 1:npoints / npoints

  if (cmep_type == "in") error <- y_err_ave_in
  else if (cmep_type == "out") error <- y_err_ave_out
  else if (cmep_type == "diff") error <- y_err_ave_out - y_err_ave_in

  plot(coverage, error, main="Coverage - Mean Absolute Error Plot", xlab="Coverage", ylab="Mean Absolute Error")
  
  if (smoothing) {
    spl <- smooth.spline(coverage, error, df=10)
    lines(spl)
  }
}

# Produce distance to model - error plot for applicability domain studies
cinf_plot_d2m_error <- function
(
  y_err,               # Prediction errors 
  d2m,                 # Distance to model
  mean_parts = TRUE,   # Whether to calcukate mean error values for parts of the points (TRUE/FALSE, 1/0)
  nparts = 3,          # The number of parts for averaging
  color_parts = "red", # Color to draw mean errors of parts
  ...
) 
{
  npoints <- length(y_err)
  
  oldnum <- order(d2m)
  y_err <- y_err[oldnum]
  d2m <- d2m[oldnum]

  plot(d2m, y_err, main="Distance to Model - Error Plot", xlab="Distance to model", ylab="Error")

  if (mean_parts) {
    p_first <- integer(nparts)
    p_last <- integer(nparts)
    p_size <- floor(npoints / nparts)
    for (ip in 1:nparts) p_first[ip] <- (ip-1)*p_size + 1
    if (nparts > 1) for (ip in 1:(nparts-1)) p_last[ip] <- p_first[ip] + p_size
    p_last[nparts] <- npoints

    p_mean <- integer(nparts)
    for (ip in 1:nparts) p_mean[ip] <- mean(y_err[p_first[ip]:p_last[ip]])

    for (ip in 1:nparts) {
      x <- c(d2m[p_first[ip]], d2m[p_last[ip]])
      y <- c(p_mean[ip],p_mean[ip])
      lines(list(x=x,y=y), col=color_parts)
    }
  }
}

# Producing plot "coverage - mean error" for applicability domain studies
# Prediction variance is taken as distance to model
cmf_ecvr_plot_coverage_mean_error <- function
(
  ecvr_fname = "ligands-ecvr.RData", # File name for resukts of external cross-validation with reshuffles
  cmep_type = "in",                  # The type of coverage - mean error plot
  smoothing = TRUE,                  # Whether to draw smoothing curve
  ...  
) 
{
  load(ecvr_fname)
  y_err <- abs(ecvr$y_exp - ecvr$y_pred_mean)
  d2m <- ecvr$y_pred_sd
#ncomp <- length(d2m)
#perm <- sample(1:ncomp, ncomp)
#d2m <- d2m[perm]
  cinf_plot_coverage_mean_error(y_err, d2m, cmep_type=cmep_type, smoothing=smoothing,...)
}

# Producing plot "distance to model - error" for applicability domain studies
# Prediction variance is taken as distance to model
cmf_ecvr_plot_d2m_error <- function
(
  ecvr_fname = "ligands-ecvr.RData", # File name for resukts of external cross-validation with reshuffles
  mean_parts = TRUE,                 # Whether to calcukate mean error values for parts of the points (TRUE/FALSE, 1/0)
  nparts = 3,                        # The number of parts for averaging
  color_parts = "red",               # Color to draw mean errors of parts
  ...  
) 
{
  load(ecvr_fname)
  y_err <- abs(ecvr$y_exp - ecvr$y_pred_mean)
  d2m <- ecvr$y_pred_sd
#ncomp <- length(d2m)
#perm <- sample(1:ncomp, ncomp)
#d2m <- d2m[perm]
  cinf_plot_d2m_error(y_err, d2m, mean_parts=mean_parts, nparts=nparts, color_parts=color_parts, ...)
}
