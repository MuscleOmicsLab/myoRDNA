#' Perform Differential Methylation Analysis at Individual Sites
#'
#' This function analyzes differential methylation at individual genomic 
#' positions.
#' It fits a logistic regression model (binomial, quasi-binomial, or 
#' beta-binomial) for each position, computes odds ratios, p-values, and 
#' predicted methylation probabilities,and applies multiple-testing correction. 
#' Positions with insufficient data are skipped.
#'
#' @param all_data A data frame containing methylation data with columns:
#'   - `position`: Genomic position.
#'   - `numCs`: Number of methylated reads.
#'   - `numTs`: Number of unmethylated reads.
#'   - `Group`: Binary treatment variable (e.g., 0 for control, 1 for 
#'   treatment).
#' @param model_type A character string specifying the statistical model to 
#' use.
#'   Options: `"beta-binomial"` (default), `"quasi-binomial"`, or 
#'   `"binomial"`.
#' @param p_adjust_method A character string specifying the method for p-value
#'  adjustment. Default is `"BH"` (Benjamini-Hochberg).
#'
#' @return A data frame with columns:
#'   - `position`: Genomic position.
#'   - `beta_0`: Intercept term from the model.
#'   - `beta_1`: Coefficient for the methylation difference.
#'   - `odds_ratio`: Odds ratio for the methylation difference.
#'   - `p_value`: Raw p-value for the methylation difference.
#'   - `p_adjust`: Adjusted p-value after multiple-testing correction.
#'   - `p_cond_1`: Predicted methylation probability for the control group.
#'   - `p_cond_2`: Predicted methylation probability for the treatment group.
#'   - `meth_diff`: Methylation difference between groups (in percentage 
#'   points).
#'   - `log10_p_adjust`: Negative log10 of the adjusted p-value.
#'
diff_meth_site <- function(
    all_data,
    model_type = "beta-binomial",
    p_adjust_method = "BH"
) {
  # Initialize results data frame
  site_results <- data.frame(
    position = numeric(),
    beta_0 = numeric(),
    beta_1 = numeric(),
    odds_ratio = numeric(),
    p_value = numeric(),
    stringsAsFactors = FALSE
  )

  # Loop through unique positions
  for (pos in unique(all_data$position)) {

    # Extract data for this position
    pos_data <- all_data[all_data$position == pos, ]

    # Skip if insufficient data
    if (nrow(pos_data) < 3 || length(unique(pos_data$Group)) < 2) {
      next
    }

    # Fit the model based on model_type
    if (model_type == "beta-binomial") {
      tryCatch({
        model <- glmmTMB::glmmTMB(
          cbind(numCs, numTs) ~ Group,
          family = glmmTMB::betabinomial(link = "logit"),
          data = pos_data,
          control = glmmTMB::glmmTMBControl(
            optimizer = stats::optim,
            optArgs = list(
              method = "Nelder-Mead"
            ),
            optCtrl = list(
              maxit = 5000
            )
          )
        )
        beta_0 <- glmmTMB::fixef(model)$cond["(Intercept)"]
        beta_1 <- glmmTMB::fixef(model)$cond["Group"]
        p_value <- base::summary(model)$coefficients$cond["Group", "Pr(>|z|)"]
      }, warning = function(w) {
        message(
          "At position ", pos, ": ",
          conditionMessage(w)
        )
        NULL
      }, error = function(e) {
        message(
          "At position ", pos, ": ",
          conditionMessage(e)
        )
        NULL
      })
    } else if (model_type == "quasi-binomial") {
      model <- stats::glm(
        cbind(numCs, numTs) ~ Group,
        family = stats::quasibinomial(link = "logit"),
        data = pos_data
      )
      beta_0 <- stats::coef(model)["(Intercept)"]
      beta_1 <- stats::coef(model)["Group"]
      p_value <- base::summary(model)$coefficients["Group", "Pr(>|t|)"]
    } else if (model_type == "binomial") {
      model <- stats::glm(
        cbind(numCs, numTs) ~ Group,
        family = stats::binomial(link = "logit"),
        data = pos_data
      )
      beta_0 <- stats::coef(model)["(Intercept)"]
      beta_1 <- stats::coef(model)["Group"]
      p_value <- base::summary(model)$coefficients["Group", "Pr(>|z|)"]
    } else {
      stop("Invalid model_type. Use 'beta-binomial', 'quasi-binomial', 
           or 'binomial'.")
    }

    # Calculate odds ratio
    odds_ratio <- exp(beta_1)

    # Store results
    site_results <- rbind(
      site_results,
      data.frame(
        position = pos,
        beta_0 = beta_0,
        beta_1 = beta_1,
        odds_ratio = odds_ratio,
        p_value = p_value
      )
    )
  }

  # Remove rows with NA values and warn the user
  na_rows <- sum(is.na(site_results$p_value))
  if (na_rows > 0) {
    warning(
      na_rows,
      " rows with NA values were removed. ",
      "Consider using a different model (e.g., quasi-binomial ",
      "or beta-binomial) to handle overdispersion.",
      call. = FALSE
    )
    site_results <- stats::na.omit(site_results)
  }

  # Multiple-testing correction
  site_results$p_adjust <- stats::p.adjust(
    site_results$p_value,
    method = p_adjust_method
  )

  # Calculate predicted methylation probabilities
  site_results$p_cond_1 <- stats::plogis(site_results$beta_0)
  site_results$p_cond_2 <- stats::plogis(site_results$beta_0 + site_results$beta_1)

  # Methylation difference in percentage points
  site_results$meth_diff <- 100 * (site_results$p_cond_2 - site_results$p_cond_1)

  # Calculate -log10(p_adjust)
  site_results$log10_p_adjust <- -log10(site_results$p_adjust)

  return(site_results)
}