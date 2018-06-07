# Grid view

require(rgl)
require(misc3d)

grid_view_level <- function(grid, level, alpha=1.0, color_p=PT$Color[["O"]], color_n=PT$Color[["N"]], ...) {
  positive_present <- FALSE
  negative_present <- FALSE
  for (igridx in 1:grid$ngridx) {
    for (igridy in 1:grid$ngridy) {
      for (igridz in 1:grid$ngridz) {
        v <- grid$val[igridx,igridy,igridz]
		if (v > level) {
		  positive_present <- TRUE
		}
		if (-v > level) {
		  negative_present <- TRUE
		}
		if (positive_present && negative_present) break
      }
    }
  }
  val <- grid$val
  if (positive_present) {
    contour3d(grid$val, level, grid$gridx, grid$gridy, grid$gridz, color=color_p, add=TRUE, alpha=alpha, ...)
  }
  if (negative_present) {
    contour3d(-val, level, grid$gridx, grid$gridy, grid$gridz, color=color_n, add=TRUE, alpha=alpha, ...)
  }
}

grid_view_part <- function(grid, part=0.01, alpha=1.0, color_p=PT$Color[["O"]], color_n=PT$Color[["N"]], ...) {
  npoints <- grid$ngridx * grid$ngridy * grid$ngridz
  values <- rep(0, npoints)
  ivalue <- 1
  for (igridx in 1:grid$ngridx) {
    for (igridy in 1:grid$ngridy) {
      for (igridz in 1:grid$ngridz) {
        values[ivalue] <- abs(grid$val[igridx,igridy,igridz])
        ivalue <- ivalue + 1
      }
    }
  }
  svalues <- sort(values, decreasing=TRUE)
  level <- svalues[ceiling(npoints * part)]
  cat(sprintf("level=%g\n", level))
  grid_view_level(grid, level, alpha, color_p, color_n, ...)
}

grid_view_rlevel <- function(grid, rlevel=0.5, alpha=1.0, color_p=PT$Color[["O"]], color_n=PT$Color[["N"]], ...) {
  max_level <- 0
  for (igridx in 1:grid$ngridx) {
    for (igridy in 1:grid$ngridy) {
      for (igridz in 1:grid$ngridz) {
        value <- abs(grid$val[igridx,igridy,igridz])
        if (value > max_level)
		  max_level <- value
      }
    }
  }
  level <- max_level * rlevel
  grid_view_level(grid, level, alpha, color_p, color_n, ...)
}
