#' Extract rDNA Chromosome Data
#'
#' Extracts all data associated with the rDNA chromosome from a `meth_object`.
#'
#' @name get_rDNA_data
#' @param meth_object A methylation object (e.g., from `methylKit`).
#' @param num_samples Number of samples in `meth_object` to process.
#' @param rDNA_chr_name Name of the rDNA chromosome (default: `"rDNAm"`).
#'
#' @return A subset of `meth_object` containing only the rDNA chromosome data.
#' 
#' @examples
#' # Load the example mouse WGBS dataset
#' file_path <- system.file(
#'   "extdata",
#'   "filtered_meth_raw_chr1_chrRDNAm_small.rds",
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
#' num_samples <- length(treatment)
#'
#' # Construct a methylRawList object
#' loaded_data <- methylKit::methylRawList(
#'   loaded_data,
#'   treatment = treatment
#' )
#'
#' # Extract methylation data for the rDNA chromosome
#' meth_rDNA <- get_rDNA_data(
#'   meth_object = loaded_data,
#'   num_samples = num_samples,
#'   rDNA_chr_name = "chrRDNAm"
#' )
#'
#' # Inspect the extracted data
#' print(meth_rDNA)
#' 
#' @export
get_rDNA_data <- function(
    meth_object,
    num_samples,
    rDNA_chr_name = "rDNAm"
) {
    rdna_gr <- GenomicRanges::GRanges(
        seqnames = rDNA_chr_name,
        ranges = IRanges::IRanges(start = 1, end = 1e9)
    )
    meth_rdna_raw <- methylKit::selectByOverlap(meth_object, rdna_gr)
    return(meth_rdna_raw)
}