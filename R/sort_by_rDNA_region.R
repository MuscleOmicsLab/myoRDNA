#' Sort Methylation Data by rDNA Regions
#'
#' @description
#' Extracts and organizes methylation data for specified rDNA regions across 
#' samples.
#' Supports predefined species (`mm39`, `rno7`). Validates sample bounds and 
#' required columns.
#'
#' @param rDNA_meth_object A list of `methylRaw` objects (from `methylKit`) or 
#' data frames.
#' @param num_samples Integer. Number of samples to process (must be >= 1 and 
#' <= `length(rDNA_meth_object)`).
#' @param species Character. Species reference for rDNA regions: `"mm39"` 
#' (default), `"rno7"`, or `"hg38"` (unsupported).
#' @param min_cov Numeric. Minimum coverage threshold to filter data 
#' (default: `10`).
#'
#' @return A named list of lists, where each sublist contains:
#'   - `numCs`: Number of methylated Cs in the region.
#'   - `numTs`: Number of unmethylated Ts in the region.
#'   - `coverage`: Coverage depth for the region.
#'   - `position`: Start positions of the region.
#'   Structured as `sample_id -> region -> data`.
#'   
#'   
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
#' loaded_data <- methods::new(
#'   "methylRawList",
#'   loaded_data,
#'   treatment = treatment
#' )
#'
#' # Extract methylation data from the rDNA chromosome
#' meth_rDNA <- get_rDNA_data(
#'   loaded_data,
#'   num_samples,
#'   rDNA_chr_name = "chrRDNAm"
#' )
#'
#' # Separate rDNA methylation data into annotated regions
#' region_res <- sort_by_rDNA_region(
#'   meth_rDNA,
#'   num_samples,
#'   species = "mm39",
#'   min_cov = 10
#' )
#'
#' # Inspect the regions available for the first sample
#' names(region_res[[1]])
#'
#' @note
#' - For `species = "hg38"`, the function stops with a warning 
#' (unsupported).
#' - Missing columns in input data will cause an error.
#'
#' @export
sort_by_rDNA_region <- function(
    rDNA_meth_object,
    num_samples,
    species = "mm39",
    min_cov = 10
) {
  if (species == "mm39") {
    message("Species selected: mm39")
    rdna_regions <- data.frame(
    region = c("18S", "5.8S", "28S", 
            "45S", "5_ETS", "ITS_1", 
            "ITS_2","3_ETS","r"),
    start = c(4008, 6878, 8123, 
            1, 1, 5878, 
            7035, 12852, 13404),
    end = c(5877, 7034, 12852, 
            13403, 4007, 6877, 
            8122, 13403, 30000)
    )
  } else if (species == "rno7") {
      message("Species selected: rno7")
      rdna_regions <- data.frame(
          region = c("18S", "5.8S", "28S"),
          start = c(1, 1873, 2030),
          end = c(1872, 2029, 6836)
      )
  } else if (species == "hg38") {
      stop(
        "Species selected: hg38. ",
        "This species is not supported. Please select another species.",
        call. = FALSE
      )
      rdna_regions <- data.frame(
        region = c("18S", "5.8S", "28S", 
              "45S", "5_ETS", "ITS_1", 
              "ITS_2","3_ETS","r"),
        start = c(4008, 6878, 8123, 
              1, 1, 5878, 
              7035, 12852, 13404),
        end = c(5877, 7034, 12852, 
              13403, 4007, 6877, 
              8122, 13403, 30000)
      )
  }
  if (num_samples < 1 || num_samples > length(rDNA_meth_object)) {
    
    stop(
      "num_samples must be between 1 and ",
      length(rDNA_meth_object)
    )
    
  }
  
  region_results <- list()
  for (i in seq_len(num_samples)) {
    
    # Get methylation object for this sample
    sample <- rDNA_meth_object[[i]]
    
    # Get sample ID
    sample_id <- sample@sample.id
    
    #print(paste("Processing sample:", sample_id))
    
    sample_data <- methylKit::getData(sample)
    
    required_columns <- c(
      "start",
      "numCs",
      "numTs",
      "coverage"
    )
    
    missing_columns <- setdiff(
      required_columns,
      colnames(sample_data)
    )
    
    if (length(missing_columns) > 0) {
      
      stop(
        "Sample ", sample_id,
        " is missing the following columns: ",
        paste(missing_columns, collapse = ", ")
      )
      
    }
    
    sample_data <- sample_data[
      sample_data$coverage >= min_cov,
    ]
    
    region_results[[sample_id]] <- stats::setNames(
      vector("list", nrow(rdna_regions)),
      rdna_regions$region
    )
    
    for (r in seq_len(nrow(rdna_regions))) {
      
      region_name <- rdna_regions$region[r]
      
      region_start <- rdna_regions$start[r]
      
      region_end <- rdna_regions$end[r]
      
      region_data <- sample_data[
        sample_data$start >= region_start &
          sample_data$start <= region_end,
      ]
      
      region_results[[sample_id]][[region_name]] <- list(
        
        numCs = region_data$numCs,
        
        numTs = region_data$numTs,
        
        coverage = region_data$coverage,
        
        position = region_data$start
        
      )
      
    }
    
  }
  return(region_results)
}