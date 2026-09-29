file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds",
  package = "myoRDNA"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

# Test
test_that("test calculate the relative rDNA copy number", {

  file_path <- system.file(
    "extdata",
    "rel_out.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- relative_rDNA_CN(
    Meth_object = loaded_data,
    num_samples = num_samples,
    autosomes = TRUE,
    rDNA_chr = "chrRDNAm",
    diploid = TRUE
  )
  
  #print(actual_output)
  
  #saveRDS(actual_output, "inst/extdata/rel_out.rds")
  
  # Compare results
  expect_equal(actual_output, expected_output)
})