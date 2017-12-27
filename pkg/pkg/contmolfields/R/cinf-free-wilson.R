#
# Processes data table using the Free-Wilson approach
# Returns frame with descriptors and properties
# t - data frame with table
# activity_column - whether the frame contains a column with activity values
#
process_free_wilson <- function(t, activity_column=TRUE) {

  # get the number of columnds in data frame
  ncols <- length(t)

  # get the number of substituent positions
  if (activity_column) {
    nsubpos <- ncols - 1
  } else {
    nsubpos <- ncols
  }

  # get column names
  colnames <- attr(t, "names")

  if (activity_column) {
    # get the activity name
    actname <- colnames[ncols]
  }

  # get the number of compounds
  ncomp <- length(t[[1]])

  # get substituent position names
  subposnames <- array()
  for (i in 1:nsubpos) {
    subposnames[i] <- colnames[i]
  }

  # find substituent list
  sublist <- as.character(levels(t[[1]]))
  nsub <- length(sublist)

  # find the number of descriptors
  ndes = nsubpos * nsub

  # generate descriptors and descriptor names
  desnames <- array()
  descr <- array(0, c(ncomp, ndes))
  ipos <- 0
  for (r in 1:nsubpos) {
    for (s in 1:nsub) {
      ipos <- ipos + 1
	  desnames[ipos] <- sprintf("%s-%s", subposnames[r], sublist[s])
	  for (c in 1:ncomp) {
	    if (as.integer(t[[r]][c]) == as.integer(s)) {
	      descr[c,ipos] <- 1
	    }
	  }
    }
  }

  # Form data frame
  df <- data.frame(descr)
  for (i in 1:ndes) {
    names(df)[i] <- desnames[i]
  }
  if (activity_column) {
    df <- cbind(df, t[ncols])
  }
  df
}

#
# Makes predictions for MLR model
# t - initial table used to build the model
# model - MLR model
# substlist - list of substitutions, example: list(R1="H", R2="Br")
#
predict_free_wilson_mlr <- function(t, model, substlist) {
  # Analyze initial table

  # get the number of columnds in data frame
  ncols <- length(t)

  # get the number of substituent positions
  nsubpos <- ncols - 1

  # get column names
  colnames <- attr(t, "names")

  # get substituent position names
  subposnames <- array()
  for (i in 1:nsubpos) {
    subposnames[i] <- colnames[i]
  }

  # find substituent list
  sublist <- as.character(levels(t[[1]]))
  nsub <- length(sublist)

  # find the number of descriptors
  ndes = nsubpos * nsub
 
  # generate descriptors and descriptor names
  desnames <- array()
  descr <- array(0, c(1, ndes))
  ipos <- 0
  for (r in 1:nsubpos) {
    for (s in 1:nsub) {
      ipos <- ipos + 1
	  desnames[ipos] <- sprintf("%s-%s", subposnames[r], sublist[s])
	  for (c in 1:1) {
	    if (substlist[[r]] == sublist[s]) {
	      descr[c,ipos] <- 1
	    }
	  }
    }
  }

  # Form data frame
  dfpred <- data.frame(descr)
  for (i in 1:ndes) {
    names(dfpred)[i] <- desnames[i]
  }
  dfpred
  pred <- predict(model, dfpred, interval="prediction")
  pred
}

#
# Produces data frame with full design
# using substituents and substitution positions 
# from table t
# t - input data frame containing table
#
full_design <- function(t) {
  # Analyze initial table

  # get the number of columnds in data frame
  ncols <- length(t)

  # get the number of substituent positions
  nsubpos <- ncols - 1

  # get column names
  colnames <- attr(t, "names")

  # get substituent position names
  subposnames <- array()
  for (i in 1:nsubpos) {
    subposnames[i] <- colnames[i]
  }

  # find substituent list
  sublist <- as.character(levels(t[[1]]))
  nsub <- length(sublist)

  form_full_frame <- function(level, curr_frame) {
    if (level <= nsubpos) {
	  for (i in 1:nsub) {
	    curr_frame[1, level] <- sublist[i]
		form_full_frame(level + 1, curr_frame)
	  }
	} else {
	  full_frame <<- rbind(full_frame, curr_frame)
	}
  }
  
  full_frame <- t[0, seq(1, nsubpos)]
  curr_frame <- t[1, seq(1, nsubpos)]
  form_full_frame(1, curr_frame)
  full_frame
}

#
# Makes predictions for data frame using MLR model
# model - MLR model
# table - table with data for prediction
#
predict_table_free_wilson_mlr <- function(model, table) {
  df <- process_free_wilson(table, activity_column=FALSE)
  pred_df <- predict(model, df, interval="prediction")
  pred_table <- cbind(table, pred_df)
  pred_table  
}

#
# Extracts rows from new data frame corresponding to compounds
# not contained in old data frame
# old_df - old data frame
# new_df - new data frame
#
extract_new <- function(old_df, new_df) {
  # number of columns in old data frame 
  ncols_old_df <- length(old_df)
  
  # number of rows in old data frame
  nrows_old_df <- length(old_df[[1]])

  # number of columns in new data frame
  ncols_new_df <- length(new_df)
  
  # number of rows in new data frame
  nrows_new_df <- length(new_df[[1]])
  
  # number of substitution positions
  nsubpos <- ncols_old_df - 1
  
  # logical array indicating whether to extract each row from new_df
  to_extract <- logical(nrows_new_df)
  
  # for each row in new_df
  for (i in 1:nrows_new_df) {
    unique <- TRUE
	for (j in 1:nrows_old_df) {
	  match <- TRUE
	  for (k in 1:nsubpos) {
	    if (old_df[j,k] != new_df[i,k]) {
		  match <- FALSE
		  break
		}
	  }
	  if (match) {
	    unique <- FALSE
		break
	  }
	}
	to_extract[i] <- unique
  }
  
  new_df[to_extract,]
}

#
# Produces data frame with samples of substitutions
# taken from table t
# from multinomial distributions
# t - input data frame containing table
# ssize - sample size
#
sample_subst <- function(t, ssize) {
  # Analyze initial table

  # get the number of columnds in data frame
  ncols <- length(t)

  # get the number of substituent positions
  nsubpos <- ncols - 1

  # get column names
  colnames <- attr(t, "names")

  # get substituent position names
  subposnames <- array()
  for (i in 1:nsubpos) {
    subposnames[i] <- colnames[i]
  }

  # find substituent list
  sublist <- as.character(levels(t[[1]]))
  nsub <- length(sublist)

  # get the number of compounds
  ncomp <- length(t[[1]])

  ##
  ## MAP parameter assessment for binomial distributions
  ## 
  # count substituent occurrences
  counts <- matrix(0, nrow=nsub, ncol=nsubpos)
  for (isub in 1:nsub) {
    for (isubpos in 1:nsubpos) {
	  for (icomp in 1:ncomp) {
	    if (sublist[isub] == t[icomp,isubpos]) {
		  counts[isub,isubpos] <- counts[isub,isubpos] + 1
		}
	  }
	}
  }
  # compute MAP estimations of alpha parameters
  # of multiple multinomial distributions
  alphas_map <- matrix(0.0, nrow=nsub, ncol=nsubpos)
  for (isub in 1:nsub) {
    for (isubpos in 1:nsubpos) {
	  alphas_map[isub,isubpos] <- (counts[isub,isubpos] + 1) / (ncomp + nsub)
	}
  }
  
  full_frame <- t[0, seq(1, nsubpos)]
  curr_frame <- t[1, seq(1, nsubpos)]
  for (is in 1:ssize) {
    for (isubpos in 1:nsubpos) {
	  curr_frame[1,isubpos] <- sample(sublist, 1, replace=TRUE, prob=alphas_map[,isubpos])
	}
    full_frame <- rbind(full_frame, curr_frame)
  }

  full_frame
}

#
# Solves the inverse problem in Free-Wilson analysis
# t - initial table for Free-Wilson analysis
# model - regression model
# y - desired value of y
# ssize - sample size
#
inv_free_wilson <- function(t, model, y, ssize=1000) {
  sample_x <- sample_subst(t, ssize)
  sample_x_pred <- predict_table_free_wilson_mlr(model, sample_x)

#  to_be_selected <- sample_x_pred$lwr < y & y < sample_x_pred$upr
#  sample_x_pred[to_be_selected,]
  
  counts <- make_freq_table(t, sample_x)
  pred <- predict_table_free_wilson_mlr(model, counts)
  
  # statistical weighting
  nrows <- length(pred[[1]])
  factor <- real(nrows)
  sum_factor <- 0
  for (irow in 1:nrows) {
    # adaptation of function normal.select from package LearnBayes
	p1 <- 0.05
	x1 <- pred$lwr[irow]
	p2 <- 0.95
	x2 <- pred$upr[irow]
	sigma <- (x1 - x2)/diff(qnorm(c(p2, p1)))
	mu <- x1 - sigma * qnorm(p1)
	weight <- dnorm(y, mu, sigma) 
	factor[irow] <- pred$count[irow] * weight
	sum_factor <- sum_factor + factor[irow]
  }
  prob <- factor / sum_factor
  
  pred <- cbind(pred, factor=factor)
  pred <- cbind(pred, prob=prob)
  plot(pred$factor, type="h")
  pred
}

#
# Makes frequency table
# t - initial table
# sample - generated sample 
#
make_freq_table <- function(t, sample) {
  fdt <- full_design(t)
  nsubpos <- length(t) - 1
  ncomb <- length(fdt[[1]])
  ssize <- length(sample[[1]])
  counts <- integer(ncomb)
  for (is in 1:ssize) {
    for (irow in 1:ncomb) {
	  identical = TRUE
	  for (icol in 1:nsubpos) {
	    if (sample[is, icol] != fdt[irow, icol]) {
		  identical <- FALSE
		  break
		}
	  }
	  if (identical) {
	    counts[irow] <- counts[irow] + 1
		break
	  }
	}
  }
  cbind(fdt, count=counts)
}

test_free_wilson <- function() {
  t <- read.table("free-wilson-boehm.txt", header=TRUE)
  df <- process_free_wilson(t)
  model <- mlr(df, stepwise=TRUE)
  plot_mlr_model(model)
  print(summary(model))
}
