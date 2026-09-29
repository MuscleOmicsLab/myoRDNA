#' Calculate Log-Odds of Methylation for a Chromosome
#'
#' Computes the log-odds of methylation (log(numCs/numTs)) for each CpG site
#' in a specified chromosome across samples. Useful for assessing the 
#' suitability of logistic regression models (e.g., in `methylKit` or custom 
#' analyses).
#'
#' @param meth_object A list of `methylRaw` objects (from `methylKit`),
#'   where each object contains methylation data for a sample.
#' @param treatment A vector of treatment/condition labels for each sample.
#'   Must have the same length as `meth_object`.
#' @param chromo A character string specifying the chromosome to analyze.
#'
#' @return A data frame with columns:
#'   - `log_odds`: Log-odds of methylation (log(numCs/numTs)) for each CpG 
#'   site.
#'   - `sample_id`: Identifier for the sample.
#'   
#' @examples
#' # Load the example mouse WGBS dataset
#' file_path <- system.file(
#'   "extdata",
#'   "filtered_meth_raw_chr1_chrRDNAm.rds",
#'   package = "myoRDNA",
#'   mustWork = TRUE
#' )
#'
#' loaded_data <- readRDS(file_path)
#'
#' # Define experimental conditions
#' # 0 = OLD, 1 = YOUNG
#' treatment <- c(
#'   0, 0, 0, 0, 0, 0,
#'   1, 1, 1, 1
#' )
#'
#' # Calculate log-odds for CpG sites on the rDNA chromosome
#' log_odds <- distribution_log_odd(
#'   meth_object = loaded_data,
#'   treatment = treatment,
#'   chromo = "chrRDNAm"
#' )
#'
#' # Display the first few results
#' head(log_odds)
#'
#' @export
distribution_log_odd <- function(
    meth_object,
    treatment,
    chromo
) {

  num_samples <- length(treatment)

  # Check that the treatment vector and methylation object agree
  if (length(meth_object) != num_samples) {
    stop("Length of meth_object must equal length of treatment.")
  }

  # Store output for each sample
  log_odds_list <- vector("list", num_samples)

  for (i in seq_len(num_samples)) {

    # Extract chromosome data for sample i
    chr_data <- meth_object[[i]][meth_object[[i]]$chr == chromo, ]

    # Calculate log-odds for every CpG
    lo <- ifelse(
      chr_data$numCs == 0 | chr_data$numTs == 0,
      NA_real_,
      log(chr_data$numCs / chr_data$numTs)
    )

    # Store one row per CpG
    log_odds_list[[i]] <- data.frame(
      log_odds = lo,
      sample_id = rep(meth_object[[i]]@sample.id, length(lo)),
      stringsAsFactors = FALSE
    )
  }

  # Combine all samples
  log_odds_out <- do.call(rbind, log_odds_list)

  return(log_odds_out)
}