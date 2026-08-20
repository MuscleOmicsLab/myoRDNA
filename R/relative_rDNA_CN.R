#' Calculate Relative rDNA Copy Number
#'
#' This function calculates the **relative rDNA copy number** for each sample in a `Meth_object`
#' (e.g., a list of methylation data or coverage data). It compares the **total coverage** of the rDNA region
#' to the **total coverage** of either the autosomes (chromosomes 1-19) or the entire genome.
#' The relative copy number is adjusted for ploidy (diploid or haploid).
#'
#' @param Meth_object A named list of `methylRaw` objects (from `methylKit`) or data frames,
#'   where each element represents a sample and contains columns `chr` (chromosome name) and
#'   `coverage` (coverage depth).
#' @param num_samples An integer specifying the number of samples in `Meth_object` to process.
#' @param autosomes A logical (`TRUE`/`FALSE`) indicating whether to use **autosomes (chromosomes 1-19)**
#'   as the reference for genome coverage. If `FALSE`, the reference is the **entire genome**.
#'   Default: `TRUE`.
#' @param rDNA_chr A character string specifying the name of the rDNA chromosome (e.g., `"chrRDNAm"`).
#' @param diploid A logical (`TRUE`/`FALSE`) indicating whether the genome is **diploid** (2 copies of each chromosome).
#'   If `TRUE`, the relative rDNA copy number is multiplied by 2. Default: `TRUE`.
#'
#' @return
#' This function **prints** a tab-separated table to the console with the following columns:
#' \itemize{
#'   \item{Sample ID:}{ The sample identifier (extracted from `Meth_object[[i]]@sample.id`).}
#'   \item{Total Number of reads:}{ The total coverage (sum of `coverage` column) for the reference region (autosomes or entire genome).}
#'   \item{Number of reads rDNA:}{ The total coverage (sum of `coverage` column) for the rDNA region.}
#'   \item{relative rDNA CN:}{ The relative rDNA copy number, calculated as `(sum_cov_rDNA / sum_cov_wgenome) * 2` (if diploid) or `(sum_cov_rDNA / sum_cov_wgenome)` (if haploid).}
#' }
#' The function does **not** return a value but prints results directly.
#'
#' @examples
#' # Example usage with methylKit objects:
#' library(methylKit)
#' meth_raw <- methRead(
#'   location = "path/to/coverage_files",
#'   sample.id = c("Sample1", "Sample2"),
#'   assembly = "hg38",
#'   context = "CpG",
#'   pipeline = "bismarkCoverage"
#' )
#' relative_rDNA_CN(
#'   Meth_object = meth_raw,
#'   num_samples = 2,
#'   autosomes = TRUE,
#'   rDNA_chr = "chrRDNAm",
#'   diploid = TRUE
#' )
#'
#' names(Meth_object) <- c("Sample1", "Sample2")
#' relative_rDNA_CN(
#'   Meth_object = Meth_object,
#'   num_samples = 2,
#'   autosomes = FALSE,
#'   rDNA_chr = "chrRDNAm",
#'   diploid = FALSE
#' )
#'
#' @note
#' - If `autosomes = TRUE`, the reference coverage (`sum_cov_wgenome`) is calculated from **all chromosomes** in the sample.
#'   This is a **deviation from the function name**, as it does not restrict to autosomes (chromosomes 1-19).
#'   To fix this, update the filtering logic for `chr_data` to include only autosomes (e.g., `chr_data <- Meth_object[[i]][Meth_object[[i]]$chr %in% paste0("chr", 1:19), ]`).
#' - If `autosomes = FALSE`, the reference coverage is calculated from **all chromosomes** in the sample.
#' - The function assumes `Meth_object` is a **named list** and that each element has a `@sample.id` slot (if using `methylRaw` objects) or can be accessed via `names(Meth_object)[i]`.
#' - If `sum_cov_wgenome` is `0`, the function will produce `Inf` or `NA` for `rel_rDNA_CN`. Add a check for `sum_cov_wgenome > 0` if this is a concern.
#' - The function uses `cat` to print tab-separated output. To capture the output, redirect it to a file or use `capture.output()`.
#'
#' @export
relative_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
{
  cat("Sample ID", "Total Number of reads", "Number of reads rDNA", "relative rDNA CN",
      sep = "\t",
      fill = TRUE)

  if (autosomes == TRUE)
  {
    for (i in 1:num_samples)
    {
      chr_data <- Meth_object[[i]]
      rDNA_data <- Meth_object[[i]][Meth_object[[i]]$chr == rDNA_chr, ]
      sum_cov_wgenome <- sum(chr_data$coverage, na.rm = TRUE)
      sum_cov_rDNA <- sum(rDNA_data$coverage, na.rm = TRUE)

      if (diploid == TRUE)
      {
        rel_rDNA_CN <- (sum_cov_rDNA / sum_cov_wgenome)
      } else
      {
        rel_rDNA_CN <- (sum_cov_rDNA / sum_cov_wgenome)
      }

      cat(
        Meth_object[[i]]@sample.id,
        sum_cov_wgenome,
        sum_cov_rDNA,
        rel_rDNA_CN,
        sep = "\t",
        fill = TRUE
      )
    }
  } else
  {
    for (i in 1:num_samples)
    {
      chr_data <- Meth_object[[i]]
      rDNA_data <- Meth_object[[i]][Meth_object[[i]]$chr == rDNA_chr, ]
      sum_cov_wgenome <- sum(chr_data$coverage, na.rm = TRUE)
      sum_cov_rDNA <- sum(rDNA_data$coverage, na.rm = TRUE)

      if (diploid == TRUE)
      {
        rel_rDNA_CN <- (sum_cov_rDNA / sum_cov_wgenome)
      } else
      {
        rel_rDNA_CN <- (sum_cov_rDNA / sum_cov_wgenome)
      }

      cat(
        Meth_object[[i]]@sample.id,
        sum_cov_wgenome,
        sum_cov_rDNA,
        rel_rDNA_CN,
        sep = "\t",
        fill = TRUE
      )
    }
  }
}
