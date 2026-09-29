#' Shift rDNA Methylation Data
#'
#' This function shifts the start and end positions of rDNA methylation data
#' by a calculated shift value, ensuring all coordinates remain non-negative.
#' The shift is computed as the difference between the maximum position in the 
#' data and a user-specified cutoff. Invalid coordinates (negative after 
#' shifting)are filtered out.
#'
#' @param meth_rdna_raw A list of `methylRaw` objects (from `methylKit`),
#'   where each object contains methylation data with `start` and `end` 
#'   positions.
#' @param cutoff A numeric value used to compute the shift.
#'   The shift is calculated as `max_position - cutoff`.
#' @param treatment A character vector or factor specifying the treatment 
#' group for each sample. Used to label the output `methylRawList`.
#'
#' @return A `methylRawList` object containing the shifted methylation data,
#'   with invalid coordinates (negative after shifting) removed.
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
#' # Define sample groups
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
#' # Shift rDNA genomic coordinates
#' meth_rDNA_shifted <- shifting(
#'   meth_rdna_raw = meth_rDNA,
#'   cutoff = 35000,
#'   treatment = treatment
#' )
#' 
#' head(meth_rDNA_shifted)
#'
#' @export
shifting <- function(
    meth_rdna_raw,
    cutoff,
    treatment
) {
  # 1. Extract all positions from meth_raw
  all_positions <- unlist(lapply(meth_rdna_raw, 
                                 function(x) methylKit::getData(x)$start))
  all_positions

  # 2. Calculate the maximum position
  max_position <- max(all_positions)
  max_position

  # 3. Compute the shift
  shift <- max_position - cutoff
  shift

  meth_rdna_shifted <- lapply(meth_rdna_raw, function(x) {
    df <- methylKit::getData(x)

    df$start <- df$start - shift
    df$end   <- df$end - shift

    # Filter invalid coordinates
    df <- df[df$start >= 0 & df$end >= 0, ]

    # Rebuild methylRaw object
    new_x <- methods::new("methylRaw",
                 df,
                 sample.id = x@sample.id,
                 assembly  = x@assembly,
                 context   = x@context,
                 resolution = x@resolution)

    return(new_x)
  })
  meth_rdna_shifted <- methylKit::methylRawList(
    meth_rdna_shifted,
    treatment = treatment
  )
  return(meth_rdna_shifted)
}