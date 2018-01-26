#
# Build multiple linear regression model
# df - frame with data (output is in the last column)
# stepwise - whether to perform stepwise descriptor selection
#
mlr <- function(df, stepwise=TRUE) {
  names <- names(df)
  ncols <- length(names)
  propname <- names[ncols]
  formula <- paste(propname, "~ .")
  model <- lm(formula, df)
  if (stepwise) model <- step(model)
  plot_regr(fitted.values(model), fitted.values(model) + residuals(model))
  model
}

#
# Plots regression
# x,y - vectors of data
# xlab - label of axis X
# ylab - label of axis Y
#
plot_regr <- function(x, y, xlab="Predicted", ylab="Experimental",
  rparams=TRUE) {
  plot(x, y, xlab=xlab, ylab=ylab)
  abline(coef=c(0,1))
  if (rparams) {
    rp <- regr_param(x, y)
    text <- sprintf("R2=%6.4f RMSE=%.3g MAE=%.3g", rp$R2, rp$RMSE, rp$MAE)
	title(sub=text)
  }
}

#
# Makes scatter plot for MLR model
# model - MLR model
#
plot_mlr_model <- function(model) {
  plot_regr(fitted.values(model), fitted.values(model) + residuals(model))
}

#
# Computes regression parameters: R2, RMSE, RMSE_pc, MAE
# predval - array of predicted values
# expval - array of experimental values
#
regr_param <- function(predval, expval) {
  n <- length(predval)
  ss <- var(expval) * (n - 1)
  rss <- 0
  rsa <- 0
  for (i in 1:n) {
    rss <- rss + (predval[i] - expval[i]) * (predval[i] - expval[i])
	rsa <- rsa + abs(predval[i] - expval[i])
  }
  r2 <- (ss - rss) / ss
  rmse <- sqrt(rss / n)
  rmse_pc <- 100 * rmse / sd(expval)
  mae <- rsa / n
  list(R2=r2, RMSE=rmse, RMSE_pc=rmse_pc, MAE=mae)
}

#
# Computes extended set of regression parameters: R2ex
# predval - array of predicted values
# expval - array of experimental values
# trainval - variance of y in training set
#
regr_param_ex <- function(predval, expval, trainval) {
  n_p <- length(predval)
  rss <- 0
  for (i in 1:n_p) {
    rss <- rss + (predval[i] - expval[i]) * (predval[i] - expval[i])
  }
  r2ex <- 1 - ((rss/(n_p-1)) / var(trainval))
  r2ex
}


