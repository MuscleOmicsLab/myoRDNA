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

num_samples <- length(treatment)

meth_rDNA <- get_rDNA_data(
    loaded_data,
    num_samples,
    rDNA_chr_name = "chrRDNAm")

# Test
test_that("test calculate regions", {

  file_path <- system.file(
    "extdata",
    "region_res.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- sort_by_rDNA_region(
    meth_rDNA,
    num_samples,
    species = "mm39",
    min_cov = 10
  )
  
  # Compare results
  expect_equal(actual_output, expected_output)
})

#region_res <- sort_by_rDNA_region(
#    meth_rDNA,
#    num_samples,
#    species = "rno7",
#    min_cov = 10
#)

#region_res <- sort_by_rDNA_region(
#    meth_rDNA,
#    num_samples,
#    species = "hg38",
#    min_cov = 10
#)