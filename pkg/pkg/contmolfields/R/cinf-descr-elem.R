# Descriptor block ELEM

elem <- function(moldbase) {
  nC_vec <- integer()
  nH_vec <- integer()
  nN_vec <- integer()
  nO_vec <- integer()
  nS_vec <- integer()
  nP_vec <- integer()
  nSe_vec <- integer()
  nF_vec <- integer()
  nCl_vec <- integer()
  nBr_vec <- integer()
  nI_vec <- integer()
  nmol <- length(moldbase)
  for (imol in 1:nmol) {
    nC <- 0
	nH <- 0
	nN <- 0
	nO <- 0
	nS <- 0
	nP <- 0
	nSe <- 0
	nF <- 0
	nCl <- 0
	nBr <- 0
	nI <- 0
	natoms <- length(moldbase[[imol]]$atoms)
	for (iatom in 1:natoms) {
	  atom <- moldbase[[imol]]$atoms[[iatom]]
	  nH <- nH + atom$nh
	  if (atom$el == "C") nC <- nC + 1
	  if (atom$el == "H") nH <- nH + 1
	  if (atom$el == "N") nN <- nN + 1
	  if (atom$el == "O") nO <- nO + 1
	  if (atom$el == "S") nS <- nS + 1
	  if (atom$el == "P") nP <- nP + 1
	  if (atom$el == "Se") nSe <- nSe + 1
	  if (atom$el == "F") nF <- nF + 1
	  if (atom$el == "Cl") nCl <- nCl + 1
	  if (atom$el == "Br") nBr <- nBr + 1
	  if (atom$el == "I") nI <- nI + 1
	}
	nC_vec[imol] <- nC
	nH_vec[imol] <- nH
	nN_vec[imol] <- nN
	nO_vec[imol] <- nO
	nS_vec[imol] <- nS
	nP_vec[imol] <- nP
	nSe_vec[imol] <- nSe
	nF_vec[imol] <- nF
	nCl_vec[imol] <- nCl
	nBr_vec[imol] <- nBr
	nI_vec[imol] <- nI
  }
  df <- data.frame(
    nC=nC_vec, 
	nH=nH_vec, 
	nN=nN_vec,
	nO=nO_vec,
	nS=nS_vec,
	nP=nP_vec,
	nSe=nSe_vec,
	nF=nF_vec,
	nCl=nCl_vec,
	nBr=nBr_vec,
	nI=nI_vec
  )
  df
}