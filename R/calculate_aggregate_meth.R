#' Computes the aggregate methylation for a specified genomic region across 
#' samples, grouped by their respective conditions.
#'
#' @param region_res A nested list where each top-level element is a sample 
#' (e.g., `region_res$SRR28431752$18S$numCs`).
#' @param region A string specifying the region of interest (e.g., `"18S"`).
#' @param conditions A named list where each name is a condition 
#' (e.g., `CONTROL`, `CONDITION 1`) and each element is a vector of sample IDs 
#' belonging to that condition.
#' @param med A logical value indicating whether to compute the median 
#' (if `TRUE`) or the mean (if `FALSE`) for the aggregate methylation values. 
#' Default is `FALSE`.
#' 
#' @return A dataframe with columns:
#' \describe{
#'   \item{sample_ID:}{The sample identifier.}
#'   \item{condition:}{The condition to which the sample belongs.}
#'   \item{aggregate_meth:}{The aggregate methylation value for the specified 
#'   region.}
#' }
#' 
#' @examples
#' # Load the example mouse WGBS dataset
#' file_path <- system.file(
#'   "extdata",
#'   "example_methylation_small.rds",
#'   package = "myoRDNA",
#'   mustWork = TRUE
#' )
#'
#' loaded_data <- readRDS(file_path)
#'
#' # Define the experimental groups
#' samples <- list(
#'   OLD = c(
#'     "SRR28431752", "SRR28431753", "SRR28431754",
#'     "SRR28431755", "SRR28431756", "SRR28431757"
#'   ),
#'   YOUNG = c(
#'     "SRR28431748", "SRR28431749",
#'     "SRR28431750", "SRR28431751"
#'   )
#' )
#'
#' treatment <- c(
#'   0, 0, 0, 0, 0, 0,
#'   1, 1, 1, 1
#' )
#'
#' num_samples <- length(treatment)
#'
#' # Construct a methylRawList object
#' loaded_data <- methods::new(
#'   "methylRawList",
#'   loaded_data,
#'   treatment = treatment
#' )
#'
#' # Extract rDNA methylation data
#' meth_rDNA <- get_rDNA_data(
#'   loaded_data,
#'   num_samples,
#'   rDNA_chr_name = "chrRDNAm"
#' )
#'
#' # Adjust genomic coordinates
#' meth_rDNA_shifted <- shifting(
#'   meth_rdna_raw = meth_rDNA,
#'   cutoff = 35000,
#'   treatment = treatment
#' )
#'
#' # Separate methylation data into rDNA regions
#' region_res <- sort_by_rDNA_region(
#'   meth_rDNA_shifted,
#'   num_samples,
#'   species = "mm39",
#'   min_cov = 10
#' )
#'
#' # Calculate mean methylation for the 18S region
#' mean_meth <- calculate_aggregate_meth(
#'   region_res,
#'   region = "18S",
#'   conditions = samples,
#'   med = FALSE
#' )
#'
#' # Calculate median methylation for the 18S region
#' median_meth <- calculate_aggregate_meth(
#'   region_res,
#'   region = "18S",
#'   conditions = samples,
#'   med = TRUE
#' )
#'
#' # Inspect the results
#' print(mean_meth)
#' print(median_meth)
#'
#' @export
calculate_aggregate_meth <- function(region_res, region, conditions,
                                     med = FALSE) {
  # Initialize an empty dataframe to store results
  result_df <- data.frame(
    Sample_ID = character(),
    Condition = character(),
    aggregate_meth = numeric(),
    stringsAsFactors = FALSE
  )

  # Iterate over each condition and its samples
  for (condition in names(conditions)) {
    for (sam in conditions[[condition]]) {
      # Extract numCs and coverage for the specified region
      numCs <- region_res[[sam]][[region]]$numCs
      cov <- region_res[[sam]][[region]]$coverage

      # Calculate aggregate_meth
      if (med == TRUE){
        aggregate_meth <- stats::median(numCs / cov, na.rm = TRUE)
      } else {
        aggregate_meth <- mean(numCs / cov, na.rm = TRUE)
      }
      
      # Append the results to the dataframe
      result_df <- rbind(
        result_df,
        data.frame(
          Sample_ID = sam,
          Condition = condition,
          aggregate_meth = aggregate_meth,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Return the dataframe
  return(result_df)
}