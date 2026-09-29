file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm_xz.rds",
  package = "myoRDNA"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

# Test
test_that("test get unique chromosome names", {

  file_path <- system.file(
    "extdata",
    "unique_names.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- get_unique_chr_names(
    loaded_data,
    num_samples)
  
  # Compare results
  expect_equal(actual_output, expected_output)
})