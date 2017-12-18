# Descriptor block CONNECT

connect <- function(moldbase) {

  if(!exists("PT"))data(PT)
  
  T <- 0.00001
  
  Chi0V_vec <- numeric()
  Chi1V_vec <- numeric()
  Chi2V_vec <- numeric()
  Chi3cV_vec <- numeric()
  Chi3pV_vec <- numeric()
  Chi4pV_vec <- numeric()
  Chi4pcV_vec <- numeric()
  Chi5cV_vec <- numeric()
  Chi5pV_vec <- numeric()
  Chi6pV_vec <- numeric()
  
  nmol <- length(moldbase)
  for (imol in 1:nmol) {
    natoms <- length(moldbase[[imol]]$atoms)
	nbonds <- length(moldbase[[imol]]$bonds)
	mol <- moldbase[[imol]]
	
	delta_v <- numeric()
	for (iatom in 1:natoms) {
	  el <- mol$atoms[[iatom]]$el
	  nh <- mol$atoms[[iatom]]$nh
	  delta_v[iatom] <- (PT$NumValEl[[el]] - nh) / (PT$NumEl[[el]] - PT$NumValEl[[el]] - 1)
	}
	
    # Compute Chi0V
	Chi0V <- 0
	for (iatom in 1:natoms) {
	  if (abs(delta_v[iatom]) > T) {
	    Chi0V <- Chi0V + 1 / sqrt(delta_v[iatom])
	  }
	}
	Chi0V_vec[imol] <- Chi0V
	
	# Compute Chi1V
	Chi1V <- 0
	if (nbonds > 0) {
	  for (ibond in 1:nbonds) {
	    at1 <- mol$bonds[[ibond]]$at1
	    at2 <- mol$bonds[[ibond]]$at2
		if (abs(delta_v[at1] * delta_v[at2]) > T) {
		  Chi1V <- Chi1V + 1 / sqrt(delta_v[at1] * delta_v[at2])
		}
	  }
	}
	Chi1V_vec[imol] <- Chi1V
	
	# Compute Chi2V, Chi4pV
	Chi2V <- 0
	Chi4pV <- 0
	Chi6pV <- 0
	if (natoms >= 3) {
	  for (iatom in 1:natoms) {
	    atomi <- mol$atoms[[iatom]]
		if (atomi$vd_ > 1) {
		  for (j in 1:(atomi$vd_-1)) {
		    a1 <- atomi$ne_[j]
			for (k in (j+1):atomi$vd_) {
			  a2 <- atomi$ne_[k]
			  if (abs(delta_v[a1]*delta_v[iatom]*delta_v[a2]) > T) {
			    Chi2V <- Chi2V + 1/sqrt(delta_v[a1]*delta_v[iatom]*delta_v[a2])
			  }
			  if (natoms >= 5) {
			    atom1 <- mol$atoms[[a1]]
				atom2 <- mol$atoms[[a2]]
				if ((atom1$vd_ >= 2) && (atom2$vd_ >= 2)) {
				  for (l in 1:atom1$vd_) {
				    a3 <- atom1$ne_[l]
					for (m in 1:atom2$vd_) {
					  a4 <- atom2$ne_[m]
					  if ((a3 != iatom) && (a3 != a2) && (a3 != a4) && (a4 != iatom) && (a4 != a1)) {
					    if (abs(delta_v[a3]*delta_v[a1]*delta_v[iatom]*delta_v[a2]*delta_v[a4]) > T) {
						  Chi4pV <- Chi4pV + 1/sqrt(delta_v[a3]*delta_v[a1]*delta_v[iatom]*delta_v[a2]*delta_v[a4])
						}
						if (natoms >= 7) {
						  atom3 <- mol$atoms[[a3]]
						  atom4 <- mol$atoms[[a4]]
						  if ((atom3$vd_ >= 2) && (atom4$vd_ >= 2)) {
						    for (n in 1:atom3$vd_) {
							  a5 <- atom3$ne_[n]
							  for (o in 1:atom4$vd_) {
							    a6 <- atom4$ne_[o]
								if ((a5 != iatom) && (a5 != a1) && (a5 != a2) && (a5 != a4) && (a5 != a6) &&
								    (a6 != iatom) && (a6 != a1) && (a6 != a2) && (a6 != a3)) {
								  if (abs(delta_v[a5]*delta_v[a3]*delta_v[a1]*delta_v[iatom]*delta_v[a2]*delta_v[a4]*delta_v[a6]) > T) {
								    Chi6pV <- Chi6pV + 1/sqrt(delta_v[a5]*delta_v[a3]*delta_v[a1]*
									                   delta_v[iatom]*delta_v[a2]*delta_v[a4]*delta_v[a6])
								  }
								}
							  }
							}
						  }
						}
					  }
					}
				  }
				}
			  }
			}
		  }
		}
	  }
	}
	Chi2V_vec[imol] <- Chi2V
	Chi4pV_vec[imol] <- Chi4pV
	Chi6pV_vec[imol] <- Chi6pV
	
	# Compute Chi3cV, Chi4pcV
	Chi3cV <- 0
	Chi4pcV <- 0
	for (iatom in 1:natoms) {
	  atomi <- mol$atoms[[iatom]]
	  if (atomi$vd_ >= 3) {
	    for (i in 1:(atomi$vd_-2)) {
		  for (j in (i+1):(atomi$vd_-1)) {
		    for (k in (j+1):atomi$vd_) {
			  a1 <- atomi$ne_[i]
			  a2 <- atomi$ne_[j]
			  a3 <- atomi$ne_[k]
			  if (abs(delta_v[iatom]*delta_v[a1]*delta_v[a2]*delta_v[a3]) > T) {
			    Chi3cV <- Chi3cV + 1/sqrt(delta_v[iatom]*delta_v[a1]*delta_v[a2]*delta_v[a3])
			  }
			}
		  }
		}
	  }
	}
	Chi3cV_vec[imol] <- Chi3cV
	Chi4pcV_vec[imol] <- Chi4pcV
	
	# Compute Chi3pV, Chi4pcV, Chi5pV, Chi5cV
	Chi3pV <- 0
	Chi4pcV <- 0
	Chi5pV <- 0
	Chi5cV <- 0
	if (natoms >= 2) for (ibond in 1:nbonds) {
	  a1 <- mol$bonds[[ibond]]$at1
	  a2 <- mol$bonds[[ibond]]$at2
	  atom1 <- mol$atoms[[a1]]
	  atom2 <- mol$atoms[[a2]]
	  if (natoms >= 4) if ((atom1$vd_ > 1) && (atom2$vd_ > 1)) {
	    for (j in 1:atom1$vd_) for (k in 1:atom2$vd_) {
		  a3 <- atom1$ne_[j]
		  a4 <- atom2$ne_[k]
		  if ((a3 != a2) && (a3 != a4) && (a4 != a1)) {
		    if (abs(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]) > T) {
			  Chi3pV <- Chi3pV + 1 /sqrt(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4])
			}
			if (natoms >= 5) if ((atom1$vd_ > 2) && (j < atom1$vd_)) for (l in (j+1):atom1$vd_) {
			  a5 <- atom1$ne_[l]
			  if ((a5 != a2) && (a5 != a4)) {
		        if (abs(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5]) > T) {
			      Chi4pcV <- Chi4pcV + 1 /sqrt(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5])
			    }
			  }
			}
			if (natoms >= 5) if ((atom2$vd_ > 2) && (k < atom2$vd_)) for (l in (k+1):atom2$vd_) {
			  a5 <- atom2$ne_[l]
			  if ((a5 != a1) && (a5 != a3)) {
		        if (abs(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5]) > T) {
			      Chi4pcV <- Chi4pcV + 1 /sqrt(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5])
			    }
			  }
			}		
		    atom3 <- mol$atoms[[a3]]
		    atom4 <- mol$atoms[[a4]]
            if (natoms >= 6) if ((atom3$vd_ > 1) && (atom4$vd_ > 1)) {
		      for (l in 1:atom3$vd_) for (m in 1:atom4$vd_) {
			    a5 <- atom3$ne_[l]
			    a6 <- atom4$ne_[m]
			    if ((a5 != a1) && (a5 != a2) && (a5 != a4) && (a5 != a6) && (a6 != a1) && (a6 != a2) && (a6 != a3)) {
		          if (abs(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5]*delta_v[a6]) > T) {
			        Chi5pV <- Chi5pV + 1 /sqrt(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5]*delta_v[a6])
			      }
			    }
			  }
            }
            if (natoms >= 6) if ((atom1$vd_ > 2) && (j < atom1$vd_) && (atom2$vd_ > 2) && (k < atom2$vd_)) {
		      for (l in (j+1):atom1$vd_) for (m in (k+1):atom2$vd_) {
			    a5 <- atom1$ne_[l]
			    a6 <- atom2$ne_[m]
			    if ((a5 != a2) && (a5 != a4) && (a5 != a6) && (a6 != a1) && (a6 != a3)) {
		          if (abs(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5]*delta_v[a6]) > T) {
			        Chi5cV <- Chi5cV + 1 /sqrt(delta_v[a1]*delta_v[a2]*delta_v[a3]*delta_v[a4]*delta_v[a5]*delta_v[a6])
			      }
			    }
			  }
			}
          }		  
		}
	  }
	}
	Chi3pV_vec[imol] <- Chi3pV
	Chi4pcV_vec[imol] <- Chi4pcV
	Chi5pV_vec[imol] <- Chi5pV
	Chi5cV_vec[imol] <- Chi5cV
  }
  
  df <- data.frame(
    Chi0V = Chi0V_vec,
	Chi1V = Chi1V_vec,
	Chi2V = Chi2V_vec,
	Chi3cV = Chi3cV_vec,
	Chi3pV = Chi3pV_vec,
	Chi4pV = Chi4pV_vec,
	Chi4pcV = Chi4pcV_vec,
	Chi5cV = Chi5cV_vec,
	Chi5pV = Chi5pV_vec,
	Chi6pV = Chi6pV_vec
  )
  df
}
