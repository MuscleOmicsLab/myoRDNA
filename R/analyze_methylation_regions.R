#' Analyze Methylation Across Genomic Regions
#'
#' This function analyzes differential methylation across specified genomic 
#' regions for a given set of samples. It uses an external script 
#' (`diff_meth_region.R`) to compute methylation differences, p-values, and 
#' adjusted p-values for each region. Results are combined into a single data 
#' frame with multiple-testing correction applied.
#'
#' @param region_res A nested list where each top-level element corresponds to 
#' a sample, and each sub-element corresponds to a genomic region 
#'  (e.g., `region_res[[1]]$18S`). Each region must contain methylation data 
#'  (e.g., `numCs`, `coverage`).
#' @param samples A vector of sample IDs to include in the analysis.
#' @param model_type A character string specifying the statistical model to 
#' use for differential methylation analysis. Default is `"beta-binomial"`.
#' @param p_adjust_method A character string specifying the method for p-value 
#' adjustment. Default is `"BH"` (Benjamini-Hochberg).
#'
#' @return A data frame with columns:
#'   - `region`: Name of the genomic region.
#'   - `beta_0`: Intercept term from the model.
#'   - `beta_1`: Coefficient for the methylation difference.
#'   - `p_value`: Raw p-value for the methylation difference.
#'   - `p_adjust`: Adjusted p-value after multiple-testing correction.
#'   - `meth_diff`: Estimated methylation difference between conditions.
#'
analyze_methylation_regions <- function(
    region_res,
    samples,
    model_type = "beta-binomial",
    p_adjust_method = "BH"
) {
  # Define regions from the names of region_res[[1]]
  regions <- names(region_res[[1]])

  # Initialize an empty data frame to store results
  out_results <- data.frame(
    region = character(),
    beta_0 = numeric(),
    beta_1 = numeric(),
    p_value = numeric(),
    p_adjust = numeric(),
    meth_diff = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop over each region
  for (reg in regions) {
    results <- diff_meth_region(
      region_res = region_res,
      region = reg,
      samples = samples,
      model_type = model_type,
      p_adjust_method = p_adjust_method
    )

    # Add the region name to the results
    results$region <- reg

    # Append results to out_results
    out_results <- rbind(out_results, results)
  }

  # Multiple-testing correction
  out_results$p_adjust <- stats::p.adjust(
    out_results$p_value,
    method = p_adjust_method
  )

  # Return the combined results
  return(out_results)
}