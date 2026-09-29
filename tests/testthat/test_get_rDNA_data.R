file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm.rds",
  package = "myoRDNA"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

loaded_data <- new(
  "methylRawList",
  loaded_data,
  treatment = treatment
)

# Test
test_that("test get rDNA data", {

  file_path <- system.file(
    "extdata",
    "meth_rDNA.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- get_rDNA_data(
    loaded_data,
    num_samples,
    rDNA_chr_name = "chrRDNAm")
  
  # Compare results
  expect_equal(actual_output, expected_output)
})