# Script for computing different classification performance measures

# Calculate kappa-statistics
calc_kappa <- function(tp, tn, fp, fn) {
  n <- tp + tn + fp + fn
  pra <- (tp + tn) / n
  pre <- ((tp+fp)*(tp+fn) + (fn+tn)*(fp+tn)) / (n*n)
  result <- (pra - pre) / (1 - pre)
}

# Calculate classification accuracy
calc_accuracy <- function(tp, tn, fp, fn) {
  n <- tp + tn + fp + fn
  result <- (tp + tn) / n
}

# Calculate classification sensitivity
calc_sensitivity <- function(tp, tn, fp, fn) {
  result <- tp / (tp + fn)
}

# Calculate classification specificity
calc_specificity <- function(tp, tn, fp, fn) {
  result <- tn / (tn + fp)
}

# Calculate classification balanced accuracy
calc_balaccuracy <- function(tp, tn, fp, fn) {
  result <- (calc_sensitivity(tp,tn,fp,fn) + calc_specificity(tp,tn,fp,fn)) / 2
}

# Calculate classification recall
calc_recall <- function(tp, tn, fp, fn) {
  result <- tp / (tp + fn)
}

# Calculate classification precision
calc_precision <- function(tp, tn, fp, fn) {
  result <- tp / (tp + fp)
}

# Calculate classification F1-score
calc_f1 <- function(tp,tn,fp,fn) {
  recall <- calc_recall(tp,tn,fp,fn)
  precision <- calc_precision(tp,tn,fp,fn)
  result <- 2*recall*precision/(recall+precision)
}

# Calculate Matthews correlation coefficient
calc_matcorcoef <- function(tp,tn,fp,fn) {
  result <- (tp*tn-fp*fn)/sqrt((tp+fn)*(fp+tn)*(tp+fp)*(fn+tn))
}

