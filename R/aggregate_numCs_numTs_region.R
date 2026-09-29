#' Computes the aggregate methylation for a specified genomic region across 
#' samples, grouped by their respective conditions.
#'
#' @param region_res A nested list where each top-level element is a sample 
#' (e.g., `region_res$SRR28431752$18S$numCs`).
#' @param region A string specifying the region of interest (e.g., `"18S"`).
#' @param conditions A named list where each name is a condition 
#' (e.g., `OLD`, `YOUNG`) and each element is a vector of sample IDs belonging 
#' to that condition.
#' @param med A logical value indicating whether to compute the median 
#' (if `TRUE`) or the mean (if `FALSE`) for the aggregate methylation values. 
#' Default is `FALSE`.
#' 
#' @return A dataframe with columns:
#' \describe{
#'   \item{sample_ID:}{The sample identifier.}
#'   \item{condition:}{The condition to which the sample belongs.}
#'   \item{aggregate_meth:}{The aggregate methylation value for the 
#'   specified region.}
#' }
#'
aggregate_numCs_numTs_region <- function(region_res, region, conditions,
                                     med = FALSE) {
  # Initialize an empty dataframe to store results
  result_df <- data.frame(
    Sample_ID = character(),
    Condition = character(),
    aggregate_numCs = numeric(),
    aggregate_numTs = numeric(),
    stringsAsFactors = FALSE
  )

  # Iterate over each condition and its samples
  for (condition in names(conditions)) {
    for (sam in conditions[[condition]]) {
      # Extract numCs and coverage for the specified region
      numCs <- region_res[[sam]][[region]]$numCs
      numTs <- region_res[[sam]][[region]]$numTs

      aggregate_numCs <- sum(numCs, na.rm = TRUE)
      aggregate_numTs <- sum(numTs, na.rm = TRUE)
      
      # Append the results to the dataframe
      result_df <- rbind(
        result_df,
        data.frame(
          Sample_ID = sam,
          Condition = condition,
          aggregate_numCs = aggregate_numCs,
          aggregate_numTs = aggregate_numTs,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Return the dataframe
  return(result_df)
}