#' Get Unique Chromosome Names and Check Consistency
#'
#' This function extracts the unique chromosome names from the first sample in a `Meth_object`
#' (e.g., a list of methylation data) and checks if all samples have the same chromosomes.
#' If inconsistencies are found, it raises an error and lists the problematic samples.
#'
#' @param Meth_object A named list of data frames or objects (e.g., `meth_raw`), where each element
#'   represents a sample and contains a `chr` column with chromosome names.
#' @param num_samples An integer specifying the number of samples in `Meth_object` to check.
#'
#' @return A character vector of unique chromosome names from the first sample, if all samples are consistent.
#'   If inconsistencies are found, the function stops execution and raises an error.
#'
#'
#' @note
#' - The function assumes `Meth_object` is a **named list** (so `names(Meth_object)[i]` returns the sample name).
#' - If `Meth_object` is unnamed, replace `names(meth_raw)[i]` with a custom identifier (e.g., `paste("Sample", i)`).
#' - The function uses `setequal()` to ensure exact matching of chromosome sets.
#'
#' @export
get_unique_chr_names <- function(
    Meth_object,
    num_samples)
{
  unique_chromosomes <- unique(Meth_object[[1]]$chr)
  inconsistent_samples <- character(0)  # Store names of inconsistent samples
  
  # Check chromosomes are consistent across all samples
  for (i in 1:num_samples) {
    current_chromosomes <- unique(Meth_object[[i]]$chr)
    
    # Check if current sample's chromosomes match the initial set
    if (!setequal(current_chromosomes, unique_chromosomes)) {
      inconsistent_samples <- c(inconsistent_samples, names(Meth_object)[i])
    }
  }
  
  # If any samples are inconsistent, raise an error with details
  if (length(inconsistent_samples) > 0) {
    stop(
      "Not all samples have the same chromosomes. ",
      "The following samples have mismatched chromosomes: ",
      paste(inconsistent_samples, collapse = ", ")
    )
  }
  
  print("The list of unique chromosomes")
  print(unique_chromosomes)
  
  # Return a list of all unique chromosomes
  return(unique_chromosomes)
}