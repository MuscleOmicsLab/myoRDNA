relative_rDNA_CN <- function(
    Meth_object,
    num_samples,
    autosomes = TRUE,
    rDNA_chr,
    diploid = TRUE)
{
  cat("Sample ID", "Total Number of reads", "Number of reads rDNA", "relative rDNA CN",
      sep = "\t",
      fill = TRUE)
  
  if (autosomes == TRUE)
  {
    for (i in 1:num_samples)
    {
      chr_data <- Meth_object[[i]]
      rDNA_data <- Meth_object[[i]][Meth_object[[i]]$chr == rDNA_chr, ]
      sum_cov_wgenome <- sum(chr_data$coverage, na.rm = TRUE)
      sum_cov_rDNA <- sum(rDNA_data$coverage, na.rm = TRUE)
      
      if (diploid == TRUE)
      {
        rel_rDNA_CN <- 2.0 * (sum_cov_rDNA / sum_cov_wgenome)
      } else
      {
        rel_rDNA_CN <- (sum_cov_rDNA / sum_cov_wgenome)
      }
      
      cat(
        Meth_object[[i]]@sample.id,
        sum_cov_wgenome,
        sum_cov_rDNA,
        rel_rDNA_CN,
        sep = "\t",
        fill = TRUE
      )
    }
  } else
  {
    for (i in 1:num_samples)
    {
      chr_data <- Meth_object[[i]]
      rDNA_data <- Meth_object[[i]][Meth_object[[i]]$chr == rDNA_chr, ]
      sum_cov_wgenome <- sum(chr_data$coverage, na.rm = TRUE)
      sum_cov_rDNA <- sum(rDNA_data$coverage, na.rm = TRUE)
      
      if (diploid == TRUE)
      {
        rel_rDNA_CN <- 2.0 * (sum_cov_rDNA / sum_cov_wgenome)
      } else
      {
        rel_rDNA_CN <- (sum_cov_rDNA / sum_cov_wgenome)
      }
      
      cat(
        Meth_object[[i]]@sample.id,
        sum_cov_wgenome,
        sum_cov_rDNA,
        rel_rDNA_CN,
        sep = "\t",
        fill = TRUE
      )
    }
  }
}