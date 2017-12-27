# Filters out columns from data frame with low variance
rm_low_var <- function(df, threshold=0.0001) {
  ind <- logical()
  ncol <- length(df)
  for (col in 1:ncol) {
    if (var(df[[col]]) < threshold) {
	  ind[col] <- FALSE
	} else {
	  ind[col] <- TRUE
	}
  }
  df[ind]
}