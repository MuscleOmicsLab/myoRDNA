file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm.rds",
  package = "myoRDNA"
)

loaded_data <- readRDS(file_path)

# Samples
samples <- list(
  'OLD' = c(
    "SRR28431752", "SRR28431753", "SRR28431754", "SRR28431755", 
    "SRR28431756", "SRR28431757"
  ),
  'YOUNG' = c(
    "SRR28431748", "SRR28431749", "SRR28431750", "SRR28431751"
  )
)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

rel_out <- relative_rDNA_CN(
    Meth_object = loaded_data,
    num_samples = num_samples,
    autosomes = TRUE,
    rDNA_chr = "chrRDNAm",
    diploid = TRUE
)

abs_out <- absolute_rDNA_CN(
    Meth_object = loaded_data,
    num_samples = num_samples,
    autosomes = TRUE,
    rDNA_chr = "chrRDNAm",
    diploid = TRUE
)

# Test
test_that("test calculate correlation coefficient", {

  file_path <- system.file(
    "extdata",
    "merged_df.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- attach_by_id(
    df1 = abs_out,
    df2 = rel_out,
    id_column = "Sample_ID",
    samples = samples,
    condition_name = "Condition",
    cor_columns = c("Absolute_rDNA_CN", "Relative_rDNA_CN")
  )
  
  # Compare results
  expect_equal(actual_output, expected_output)
})
