attach_by_id <- function(
    df1,
    df2,
    id_column = "Sample_ID",
    samples,
    condition_name = "Category",
    cor_columns = NULL
) {
  # Create a named vector to map ID values to their condition
  id_conditions <- unlist(
    lapply(names(samples), function(condition) {
      setNames(rep(condition, length(samples[[condition]])), samples[[condition]])
    })
  )
  
  # Merge the dataframes by the ID column
  merged_df <- merge(
    df1,
    df2,
    by = id_column,
    suffixes = c("_df1", "_df2")
  )
  
  # Add the condition column to the merged dataframe
  merged_df[[condition_name]] <- id_conditions[merged_df[[id_column]]]
  
  # Compute Pearson correlation if columns are specified
  if (!is.null(cor_columns)) {
    if (length(cor_columns) != 2) {
      stop("cor_columns must be a vector of exactly two column names.")
    }
    
    # Check if the specified columns exist in the merged dataframe
    if (!all(cor_columns %in% names(merged_df))) {
      stop("One or both specified columns for correlation do not exist in the merged dataframe.")
    }
    
    # Pearson correlation for all samples
    cor_test_all <- cor.test(
      merged_df[[cor_columns[1]]],
      merged_df[[cor_columns[2]]],
      method = "pearson"
    )
    print(cor_test_all)
  }
  
  # Return the merged dataframe with conditions
  return(merged_df)
}
