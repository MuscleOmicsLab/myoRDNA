file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds",
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

loaded_data <- new(
  "methylRawList",
  loaded_data,
  treatment = treatment
)

meth_rDNA <- get_rDNA_data(
  loaded_data,
  num_samples,
  rDNA_chr_name = "chrRDNAm")

cutoff = 35000
meth_rDNA_shifted <- shifting(
  meth_rdna_raw=meth_rDNA,
  cutoff=cutoff,
  treatment=treatment
)

region_res <- sort_by_rDNA_region(
  meth_rDNA_shifted,
  num_samples,
  species = "mm39",
  min_cov = 10
)

region <- '18S'
# Test mean aggregation (med = FALSE)
test_that("calculate_aggregate_meth correctly calculates mean methylation", {

  file_path <- system.file(
    "extdata",
    "aggregate_meth_avg.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- calculate_aggregate_meth(
    region_res,
    region,
    conditions = samples,
    med = FALSE
  )
  
  #print(actual_output)
  
  #saveRDS(actual_output, "inst/extdata/aggregate_meth_avg.rds")
  
  # Compare results
  expect_equal(actual_output, expected_output)
})


# Test median aggregation (med = TRUE)
test_that("calculate_aggregate_meth correctly calculates median 
          methylation", {

  file_path <- system.file(
    "extdata",
    "aggregate_meth_med.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- calculate_aggregate_meth(
    region_res,
    region,
    conditions = samples,
    med = TRUE
  )
  
  #print(actual_output)
  
  #saveRDS(actual_output, "inst/extdata/aggregate_meth_med.rds")
  
  # Compare results
  expect_equal(actual_output, expected_output)
})