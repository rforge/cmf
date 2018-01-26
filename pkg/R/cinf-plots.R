# Plotting graphs

# Scatter plot
cinf_plotxy <- function(x, y, xlab="Predicted", ylab="Experimental", margin=0.5, ...) {
  xmin <- min(x)
  xmax <- max(x)
  ymin <- min(y)
  ymax <- max(y)
  xymin <- min(xmin, ymin) - margin
  xymax <- max(xmax, ymax) + margin
  plot(x, y, xlim=c(xymin, xymax), ylim=c(xymin, xymax), xlab=xlab, ylab=ylab, ...)
  abline(coef=c(0,1))
}
