file_path <- system.file(
  "extdata",
  "filtered_meth_raw_chr1_chrRDNAm_small.rds",
  package = "myoRDNA"
)

loaded_data <- readRDS(file_path)

treatment = c(0,0,0,0,0,0,
              1,1,1,1)

num_samples <- length(treatment)

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

# Test
test_that("test site by site differential methylation", {

  file_path <- system.file(
    "extdata",
    "differential_meth_site.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- diff_meth(
    meth_data=meth_rDNA_shifted,
    samples=samples,
    species="mm39",
    num_samples=num_samples,
    min_cov=10,
    region="18S",
    model_type="binomial",
    p_adjust_method="BH",
    site = TRUE)
  
  #print(actual_output)
  
  #saveRDS(actual_output, "inst/extdata/differential_meth_site.rds")
  
  # Compare results
  expect_equal(actual_output, expected_output)
})

# Test
test_that("test region by region differential methylation", {

  file_path <- system.file(
    "extdata",
    "differential_meth_region.rds",
    package = "myoRDNA"
  )

  expect_true(nzchar(file_path))
  expect_true(file.exists(file_path))

  # Load expected test data
  expected_output <- readRDS(file_path)
  
  # Run the function
  actual_output <- diff_meth(
    meth_data=meth_rDNA_shifted,
    samples=samples,
    species="mm39",
    num_samples=num_samples,
    min_cov=10,
    region="18S",
    model_type="binomial", #beta-binomial, quasi-binomial
    p_adjust_method="BH",
    site = FALSE)
  
  #print(actual_output)
  
  #saveRDS(actual_output, "inst/extdata/differential_meth_region.rds")
  
  # Compare results
  expect_equal(actual_output, expected_output)
})