#' Calculate Absolute rDNA Copy Number
#'
#' This function calculates the **absolute rDNA copy number** for each sample in a `Meth_object`
#' (e.g., a list of methylation data or coverage data). It compares the coverage of the rDNA region
#' to the coverage of either the autosomes (chromosomes 1-19) or the entire genome (excluding rDNA).
#' The absolute copy number is adjusted for ploidy (diploid or haploid).
#'
#' @param Meth_object - A named list of `methylRaw` objects (from `methylKit`) or data frames,
#'   where each element represents a sample and contains columns `chr` (chromosome name) and
#'   `coverage` (coverage depth).
#' @param num_samples - An integer specifying the number of samples in `Meth_object` to process.
#' @param autosomes - A logical (`TRUE`/`FALSE`) indicating whether to use **autosomes (chromosomes 1-19)**
#'   as the reference for genome coverage. If `FALSE`, the reference is the entire genome **excluding rDNA**.
#'   Default: `TRUE`.
#' @param rDNA_chr - A character string specifying the name of the rDNA chromosome (e.g., `"chrRDNAm"`).
#' @param diploid - A logical (`TRUE`/`FALSE`) indicating whether the genome is **diploid** (2 copies of each chromosome).
#'   If `TRUE`, the absolute rDNA copy number is multiplied by 2. Default: `TRUE`.
#'
#' @return
#' This function **prints** a tab-separated table to the console with the following columns:
#' \itemize{
#'   \item{Sample ID:}{ The sample identifier (extracted from `Meth_object[[i]]@sample.id`).}
#'   \item{Genome Depth:}{ The mean coverage depth for the reference region (autosomes or genome excluding rDNA).}
#'   \item{rDNA Depth:}{ The mean coverage depth for the rDNA region.}
#'   \item{absolute rDNA CN:}{ The absolute rDNA copy number, calculated as `(cov_rDNA / cov_wgenome) * 2` (if diploid) or `(cov_rDNA / cov_wgenome)` (if haploid).}
#' }
#' The function does **not** return a value but prints results directly.
#'
#'
#' @note
#' - If `autosomes = TRUE`, the reference coverage (`cov_wgenome`) is calculated from chromosomes 1-19.
#'   Ensure these chromosomes exist in your data.
#' - If `autosomes = FALSE`, the reference coverage is calculated from **all chromosomes except `rDNA_chr`**.
#' - The function assumes `Meth_object` is a **named list** and that each element has a `@sample.id` slot (if using `methylRaw` objects) or can be accessed via `names(Meth_object)[i]`.
#' - If `cov_wgenome` is `0`, the function will produce `Inf` or `NA` for `abs_rDNA_CN`. Add a check for `cov_wgenome > 0` if this is a concern.
#' - The function uses `cat` to print tab-separated output. To capture the output, redirect it to a file or use `capture.output().
#'
#' @export
absolute_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
{
  # Initialize a dataframe to store results
  results_df <- data.frame(
    Sample_ID = character(),
    Genome_Depth = numeric(),
    rDNA_Depth = numeric(),
    Absolute_rDNA_CN = numeric(),
    stringsAsFactors = FALSE
  )

  # Print header
  cat("Sample ID", "Genome Depth", "rDNA Depth", "absolute rDNA CN",
      sep = "\t",
      fill = TRUE)

  if (autosomes == TRUE)
  {
    for (i in 1:num_samples)
    {
      chr_data <- Meth_object[[i]][sub("^chr", "", Meth_object[[i]]$chr) %in% as.character(1:19),]
      rDNA_data <- Meth_object[[i]][Meth_object[[i]]$chr == rDNA_chr, ]
      cov_wgenome <- mean(chr_data$coverage, na.rm = TRUE)
      cov_rDNA <- mean(rDNA_data$coverage, na.rm = TRUE)

      if (diploid == TRUE)
      {
        abs_rDNA_CN <- 2.0 * (cov_rDNA / cov_wgenome)
      } else
      {
        abs_rDNA_CN <- (cov_rDNA / cov_wgenome)
      }

      # Print the output
      cat(
        Meth_object[[i]]@sample.id,
        cov_wgenome,
        cov_rDNA,
        abs_rDNA_CN,
        sep = "\t",
        fill = TRUE
      )

      # Append results to the dataframe
      results_df <- rbind(
        results_df,
        data.frame(
          Sample_ID = Meth_object[[i]]@sample.id,
          Genome_Depth = cov_wgenome,
          rDNA_Depth = cov_rDNA,
          Absolute_rDNA_CN = abs_rDNA_CN,
          stringsAsFactors = FALSE
        )
      )
    }
  } else
  {
    for (i in 1:num_samples)
    {
      chr_data <- Meth_object[[i]]
      rDNA_data <- Meth_object[[i]][Meth_object[[i]]$chr == rDNA_chr, ]
      cov_wgenome <- mean(chr_data$coverage, na.rm = TRUE)
      cov_rDNA <- mean(rDNA_data$coverage, na.rm = TRUE)

      if (diploid == TRUE)
      {
        abs_rDNA_CN <- 2.0 * (cov_rDNA / cov_wgenome)
      } else
      {
        abs_rDNA_CN <- (cov_rDNA / cov_wgenome)
      }

      # Print the output
      cat(
        Meth_object[[i]]@sample.id,
        cov_wgenome,
        cov_rDNA,
        abs_rDNA_CN,
        sep = "\t",
        fill = TRUE
      )

      # Append results to the dataframe
      results_df <- rbind(
        results_df,
        data.frame(
          Sample_ID = Meth_object[[i]]@sample.id,
          Genome_Depth = cov_wgenome,
          rDNA_Depth = cov_rDNA,
          Absolute_rDNA_CN = abs_rDNA_CN,
          stringsAsFactors = FALSE
        )
      )
    }
  }

  # Return the dataframe
  return(results_df)
}
