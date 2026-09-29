#' Extract Methylation Site Data for a Genomic Region
#'
#' Extracts methylation counts (`numCs`, `numTs`) and positions for a 
#' specified genomic region
#' across all samples. Assigns a binary treatment variable (`Group`) based on 
#' the first group in `samples`.
#' Skips samples with mismatched vector lengths and removes zero-coverage 
#' sites.
#'
#' @param region_res A nested list where each top-level element corresponds to 
#' a sample, and each sub-element corresponds to a genomic region 
#'  (e.g., `region_res[[1]]$18S`).
#' @param samples A named list of sample IDs grouped by condition 
#' (e.g., `list(OLD = c("sample1", "sample2"), YOUNG = c("sample3"))`).
#' @param region A character string specifying the genomic region to extract. 
#' Default is `"18S"`.
#'
#' @return A data frame with columns:
#'   - `sample`: Sample identifier.
#'   - `position`: Genomic position.
#'   - `numCs`: Number of methylated reads.
#'   - `numTs`: Number of unmethylated reads.
#'   - `Group`: Binary treatment variable (0 for control group, 1 otherwise).
#' 
extract_site_data <- function(
    region_res,
    samples,
    region = "18S"
) {
  control_group <- names(samples[1])

  # Extract numCs, numTs, and positions for all samples
  numCs_reg <- lapply(region_res, function(x) x[[region]]$numCs)
  numTs_reg <- lapply(region_res, function(x) x[[region]]$numTs)
  position_reg <- lapply(region_res, function(x) x[[region]]$position)

  all_data <- data.frame()

  for (sam in names(numCs_reg)) {

    numCs <- as.numeric(numCs_reg[[sam]])
    numTs <- as.numeric(numTs_reg[[sam]])
    position <- as.numeric(position_reg[[sam]])

    # Assign T based on control_group
    Group <- ifelse(sam %in% samples[[control_group]], 0, 1)

    # Check if vectors have the same length
    if (length(numCs) != length(numTs) || length(numCs) != length(position)) {
      warning(
        "numCs, numTs, and position have different lengths for sample ",
        sam,
        ". Skipping this sample.",
        call. = FALSE
      )
      next
    }

    # Create data for this sample
    df <- data.frame(
      sample = sam,
      position = position,
      numCs = numCs,
      numTs = numTs,
      Group = Group
    )

    # Remove zero-coverage sites
    df <- df[(df$numCs + df$numTs) > 0, ]

    # Add to combined data frame
    all_data <- rbind(all_data, df)
  }

  return(all_data)
}