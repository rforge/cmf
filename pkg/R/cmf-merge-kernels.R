# Merging files with computed kernels
cmf_merge_kernels <- function(kernels_train_fname_list, kernels_train_merged_fname)
{
  kernels_merged <- list()
  for (kernels_fname in kernels_train_fname_list) {
    load(kernels_fname)
    field_names <- names(kernels)
    for (el in field_names) kernels_merged[[el]] <- kernels[[el]]
  }
  kernels <- kernels_merged
  save(kernels, file=kernels_train_merged_fname)
}

