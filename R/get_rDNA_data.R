#' Extract rDNA Chromosome Data
#'
#' Extracts all data associated with the rDNA chromosome from a `meth_object`.
#'
#' @param meth_object A methylation object (e.g., from `methylKit`).
#' @param num_samples Number of samples in `meth_object` to process.
#' @param rDNA_chr_name Name of the rDNA chromosome (default: `"rDNAm"`).
#'
#' @return A subset of `meth_object` containing only the rDNA chromosome data.
#' @examples
#' get_rDNA_data(meth_object, num_samples = 10, rDNA_chr_name = "rDNAm")
library("GenomicRanges")
get_rDNA_data <- function(
    meth_object,
    num_samples,
    rDNA_chr_name = "rDNAm"
) {
    rdna_gr <- GRanges(
        seqnames = rDNA_chr_name,
        ranges = IRanges(start = 1, end = 1e9)
    )
    meth_rdna_raw <- selectByOverlap(meth_object, rdna_gr)
    return(meth_rdna_raw)
}