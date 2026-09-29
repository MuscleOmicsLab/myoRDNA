file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm_xz.rds",
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

# Test
test_that("test log odds for numCs and numTs", {

  file_path <- system.file(
    "extdata",
    "log_odds.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- distribution_log_odd(
    loaded_data, 
    treatment, 
    chromo="chrRDNAm")
  
  
  
  # Compare results
  expect_equal(actual_output, expected_output)
})