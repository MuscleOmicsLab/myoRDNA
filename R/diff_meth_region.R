#' Perform Differential Methylation Analysis for a Genomic Region
#'
#' Aggregates methylation counts (`numCs` and `numTs`) for a specified genomic 
#' region, fits a logistic regression model (binomial, quasi-binomial, or 
#' beta-binomial), and computes methylation differences, p-values, and 
#' predicted probabilities.
#'
#' @param region_res A nested list where each top-level element corresponds to 
#' a sample, and each sub-element corresponds to a genomic region 
#'  (e.g., `region_res[[1]]$18S`).
#'   Each region must contain methylation data (e.g., `numCs`, `numTs`).
#' @param region A character string specifying the name of the genomic region 
#' to analyze.
#' @param samples A named vector of sample IDs, where names correspond to 
#' conditions
#' (e.g., 
#' `c(Young = c("sample1", "sample2"), Old = c("sample3", "sample4"))`).
#' @param model_type A character string specifying the statistical model to 
#' use. Options: `"binomial"`, `"quasi-binomial"`, or `"beta-binomial"` 
#' (default).
#' @param p_adjust_method A character string specifying the method for p-value 
#' adjustment.
#'   Default is `"BH"` (Benjamini-Hochberg). Note: Adjustment is applied in 
#'   thewrapper function.
#'
#' @return A data frame with columns:
#'   - `beta_0`: Intercept term from the model.
#'   - `beta_1`: Coefficient for the methylation difference.
#'   - `p_value`: Raw p-value for the methylation difference.
#'   - `meth_diff`: Estimated methylation difference between conditions 
#'   (in percentage points).
#' 
diff_meth_region <- function(
    region_res,
    region,
    samples,
    model_type = "beta-binomial",
    p_adjust_method = "BH"
) {
  # Aggregate the data
  data <- aggregate_numCs_numTs_region(
    region_res,
    region,
    conditions = samples
  )

  # Create binary treatment variable
  data$Group <- ifelse(data$Condition == names(samples)[1], 0, 1)

  # Fit model
  if (model_type == "binomial") {

    model <- stats::glm(
      cbind(aggregate_numCs, aggregate_numTs) ~ Group,
      data = data,
      family = stats::binomial(link = "logit")
    )

    beta_0 <- stats::coef(model)["(Intercept)"]
    beta_1 <- stats::coef(model)["Group"]
    p_value <- base::summary(model)$coefficients["Group", "Pr(>|z|)"]

  } else if (model_type == "quasi-binomial") {

    model <- stats::glm(
      cbind(aggregate_numCs, aggregate_numTs) ~ Group,
      data = data,
      family = stats::quasibinomial(link = "logit")
    )

    beta_0 <- stats::coef(model)["(Intercept)"]
    beta_1 <- stats::coef(model)["Group"]
    p_value <- base::summary(model)$coefficients["Group", "Pr(>|z|)"]

  } else if (model_type == "beta-binomial") {

    if (!requireNamespace("glmmTMB", quietly = TRUE)) {
      stop("Package 'glmmTMB' is required but not installed. 
           Please install it.")
    }

    model <- glmmTMB::glmmTMB(
      cbind(aggregate_numCs, aggregate_numTs) ~ Group,
      family = glmmTMB::betabinomial(link = "logit"),
      data = data
    )

    beta_0 <- glmmTMB::fixef(model)$cond["(Intercept)"]
    beta_1 <- glmmTMB::fixef(model)$cond["Group"]
    p_value <- base::summary(model)$coefficients$cond["Group", "Pr(>|z|)"]

  } else {
    stop("Invalid model_type. Choose 'binomial', 'quasi-binomial', or 
         'beta-binomial'.")
  }

  # Extract numeric coefficient values
  beta_0 <- unname(beta_0)
  beta_1 <- unname(beta_1)

  # Calculate predicted methylation probabilities
  p_young <- stats::plogis(beta_0)
  p_old <- stats::plogis(beta_0 + beta_1)

  # Methylation difference in percentage points
  meth_diff <- 100 * (p_old - p_young)

  # Create results data frame
  region_results <- data.frame(
    beta_0 = beta_0,
    beta_1 = beta_1,
    p_value = p_value,
    meth_diff = meth_diff
  )

  return(region_results)
}