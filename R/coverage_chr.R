#' Plot Average Coverage per Chromosome
#'
#' @description
#' Calculates the average coverage for each chromosome across all samples and 
#' plots it as a bar chart.
#'
#' @param meth_object A list of data frames or `methylRaw` objects containing 
#' methylation data.
#' @param num_samples Integer. Number of samples to process.
#' @param chromosome_names Character vector. Names of chromosomes to include 
#' in the analysis.
#'
#' @return A `ggplot` bar chart visualizing the average coverage per 
#' chromosome.
#' 
#' @examples
#' # Load the example methylation dataset
#' file_path <- system.file(
#'   "extdata",
#'   "example_methylation_small.rds",
#'   package = "myoRDNA",
#'   mustWork = TRUE
#' )
#'
#' loaded_data <- readRDS(file_path)
#'
#' # Define the number of samples
#' treatment <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1)
#' num_samples <- length(treatment)
#'
#' # Extract unique chromosome names
#' chromosome_names <- get_unique_chr_names(
#'   loaded_data,
#'   num_samples
#' )
#'
#' # Display the chromosome names
#' print(chromosome_names)
#'
#' @note
#' - The function assumes the input data frames contain a `coverage` column.
#' - The plot uses `ggplot2` for visualization.
#' 
#' @export
coverage_chr <- function(
    meth_object,
    num_samples,
    chromosome_names
) {
    results <- data.frame(
        chromosome = character(),
        sample = numeric(),
        avg_coverage = numeric(),
        stringsAsFactors = FALSE
    )
    for (chromo in chromosome_names) {
      for (i in seq_len(num_samples)) {
        chr_cov <- meth_object[[i]][meth_object[[i]]$chr == chromo, ]
        # Filter for coverage >= 10
        chr_cov_filtered <- chr_cov[chr_cov$coverage >= 10, ]
        # Calculate average coverage for filtered data
        avg_chr_cov <- base::mean(chr_cov_filtered$coverage, na.rm = TRUE)
        # Append the results to the data frame
        results <- rbind(results, data.frame(
          chromosome = chromo,
          sample = i,
          avg_coverage = avg_chr_cov
        ))
      }
    }
    
    # Calculate the average coverage for each chromosome across all samples
    avg_coverage_per_chromosome <- stats::aggregate(avg_coverage ~ chromosome,
                                          data = results,
                                          FUN = base::mean)

    return(avg_coverage_per_chromosome)
}