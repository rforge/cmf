# View molecule, co-field and field

require(rgl)

#source("R/cinf-mol2.R")
#source("cinf-molview.R")
#source("cmf-gridview.R")
#source("cmf-grid.R")
#source("cmf-fields.R")
#source("cmf-params.R")

cmf_view_mol_field_cofield <- function
(
  mdb_fname = "ligands-train.mol2",          # File name for molecular database 
  imol = 1,                                  # Molecule to visualize
  ft = "q",                                  # Field type to visualize
  grid_fname = "ligands-grid-krr.RData",     # File name for grid with co-fields
  alpha_from_model = FALSE,                  # Take alpha from model? (TRUE/FALSE, 1/0)
  model_fname = "ligands-model.RData",       # Model file name
  alpha = 0.3,                               # Alpha value (if not taken from model)
  rlevel = 0.5,                              # Isosurface level
  alpha_g = 0.7,                             # Alpha (non-transperancy) level
  draw_field = TRUE,                         # Whether to draw field
  draw_cofield = TRUE,                       # Whether to draw co-field
  draw_overlap = FALSE,                      # Whether to draw overlap between fields and co-fields
  ...
)
{
  if (alpha_from_model) {
    load(model_fname)
    alpha <- model$alpha[[ft]]
  }

  mdb0 <- read_mol2(mdb_fname)
  mdb <- cmf_params_tripos(mdb0)
  mol <- mdb[[imol]]
  mol_view_cpk(mol)
  mol_view_cylindres(mol)
  
  # Get default colors of physico-chemical fields
  fcolors <- def_fcolors()
  if (ft %in% tripos_atom_types) {
    # Colors of indixator fields
    color_p <- "red"
	color_n <- "blue"
  } else {
    # Colors of physico-chemical fields
    color_p <- fcolors$p[[ft]]
	color_n <- fcolors$n[[ft]]
  }
  
  # Draw co-field
  if (draw_cofield) {
    load(grid_fname)
    grid_view_rlevel(grids[[ft]], rlevel=rlevel, alpha=alpha_g, color_p=color_p, color_n=color_n, fill=!draw_overlap, ...)
  }
  
  # Draw field
  if (draw_field) {
    if (draw_cofield) {
	  grid_f <- grids[[ft]]
	} else {
      grid_f <- cmf_init_grid(mdb, step=1.0)
	}
    grid_f[[ft]] <- cmf_fval_grid(ft, mol, alpha, grid_f)
    grid_view_rlevel(grid_f[[ft]], rlevel=rlevel, alpha=alpha_g, color_p=color_p, color_n=color_n, fill=!draw_cofield, ...)
  }

  # Draw overlap between fields and co-fields
  if (draw_overlap) {
    grid_p <- cmf_multiply_fields(grids[[ft]], grid_f[[ft]])  
    grid_view_rlevel(grid_p, rlevel=rlevel, alpha=alpha_g, color_p=color_p, color_n=color_n, fill=TRUE, ...)
  }
  
}

# Define default colors for physico-chemical molecular fields
def_fcolors <- function() {
  fcolors <- list()
  fcolors$p[["q"]]        <- "red";    fcolors$n[["q"]]        <- "blue"
  fcolors$p[["vdw"]]      <- "green";  fcolors$n[["vdw"]]      <- "yellow"
  fcolors$p[["logp"]]     <- "yellow"; fcolors$n[["logp"]]     <- "violet"
  fcolors$p[["abra"]]     <- "red";    fcolors$n[["abra"]]     <- "blue"
  fcolors$p[["abrb"]]     <- "red";    fcolors$n[["abrb"]]     <- "blue"
  fcolors$p[["ind"]]      <- "green";  fcolors$n[["ind"]]      <- "yellow"
  fcolors
}
