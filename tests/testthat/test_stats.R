context("Testing statistics")

library(notame)
# Summary statistics ----
test_that("summary statistics work without grouping", {

  smry <- summary_statistics(mark_nas(example_set, 0))
  ex <- exprs(mark_nas(example_set, 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(ex, 1, fun)), smry[, gsub("finite_", "", fun)])
  }
  expect_equal(unname(apply(ex, 1, finite_quantile, probs = 0.25)), smry$Q25)
  expect_equal(unname(apply(ex, 1, finite_quantile, probs = 0.75)), smry$Q75)
})


test_that("summary statistics work with grouping", {
  foreach::registerDoSEQ()

  smry <- summary_statistics(mark_nas(example_set, 0), grouping_cols = "Group")
  exa <- exprs(mark_nas(example_set[, example_set$Group == "A"], 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(exa, 1, fun)), smry[, gsub("finite", "A", fun)])
  }
  expect_equal(unname(apply(exa, 1, finite_quantile, probs = 0.25)), smry$A_Q25)
  expect_equal(unname(apply(exa, 1, finite_quantile, probs = 0.75)), smry$A_Q75)

  exb <- exprs(mark_nas(example_set[, example_set$Group == "B"], 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(exb, 1, fun)), smry[, gsub("finite", "B", fun)])
  }
  expect_equal(unname(apply(exb, 1, finite_quantile, probs = 0.25)), smry$B_Q25)
  expect_equal(unname(apply(exb, 1, finite_quantile, probs = 0.75)), smry$B_Q75)
})

test_that("summary statistics work with all NA features", {

  ex_set_na <- mark_nas(example_set, 0)
  exprs(ex_set_na)[1, ] <- NA

  smry <- summary_statistics(ex_set_na)

  expect_equal(nrow(smry), nrow(exprs(ex_set_na)))
  expect_equal(smry$Feature_ID, featureNames(ex_set_na))
  expect_true(all(is.na(smry[1, 2:ncol(smry)])))
})

# Cohen's D ----
test_that("Cohen's d works", {
  ex <- example_set %>%
    drop_qcs() %>%
    mark_nas(0)
  cd <- ex %>%
    combined_data()

  cd1 <- cd[cd$Time == 1, ]
  cd2 <- cd[cd$Time == 2, ]

  d <- c()
  for (feature in featureNames(ex)) {
    tdiff <- data.frame(feature = cd2[, feature] - cd1[, feature],
                        group = cd1$Group)
    mean_a <- finite_mean(tdiff$feature[tdiff$group == "A"])
    sd_a <- finite_sd(tdiff$feature[tdiff$group == "A"])
    mean_b <- finite_mean(tdiff$feature[tdiff$group == "B"])
    sd_b <- finite_sd(tdiff$feature[tdiff$group == "B"])
    d <- c(d, (mean_b - mean_a)/sqrt(mean(c(sd_a^2, sd_b^2))))
  }

  cohd <- cohens_d(ex, id = "Subject_ID", time = "Time")

  df <- data.frame(Feature_ID = featureNames(ex),
                   Cohen_d = d,
                   stringsAsFactors = FALSE)
  rownames(df) <- df$Feature_ID
  expect_equal(cohd, df)
})

test_that("Cohen's d work with all NA features", {

  ex_set_na <- drop_qcs(mark_nas(example_set, 0))
  exprs(ex_set_na)[1:2, ] <- NA

  cohd <- cohens_d(ex_set_na)

  expect_equal(nrow(cohd), nrow(exprs(ex_set_na)))
  expect_equal(cohd$Feature_ID, featureNames(ex_set_na))
  expect_true(all(is.na(cohd[1:2, 2:ncol(cohd)])))
})

test_that("Cohen's d is not run with multiple time or group levels", {

  expect_error(cohens_d(example_set), "Column Group")

  ex <- example_set
  ex$Group <- c(1, 2)
  expect_error(cohens_d(ex, id = "Subject_ID", time = "Time"), "Column Time")

  ex$Time <- c(1,2)
  ex$Group <- 1
  expect_error(cohens_d(ex, id = "Subject_ID", time = "Time"), "Column Group")
})

# Fold change ----
test_that("Fold change works", {
  ex <- example_set %>%
    mark_nas(0)
  cd <- ex %>%
    combined_data()

  cd1 <- cd[cd$Time == 1, ]
  cd2 <- cd[cd$Time == 2, ]

  fc <- data.frame(Feature_ID = featureNames(ex),
                   FC_B_vs_A = 1,
                   FC_QC_vs_A = 1,
                   FC_QC_vs_B = 1,
                   stringsAsFactors = FALSE)
  rownames(fc) <- fc$Feature_ID
  for (i in seq_len(nrow(fc))) {
    feature <- fc$Feature_ID[i]
    mean_a <- finite_mean(cd[cd$Group == "A", feature])
    mean_b <- finite_mean(cd[cd$Group == "B", feature])
    mean_qc <- finite_mean(cd[cd$Group == "QC", feature])

    fc$FC_B_vs_A[i] <- mean_b/mean_a
    fc$FC_QC_vs_A[i] <- mean_qc/mean_a
    fc$FC_QC_vs_B[i] <- mean_qc/mean_b
  }

  foldc <- fold_change(ex)

  expect_equal(foldc, fc)
})

test_that("Fold change works with all NA features", {

  ex_set_na <- drop_qcs(mark_nas(example_set, 0))
  exprs(ex_set_na)[1:2, ] <- NA

  foldc <- fold_change(ex_set_na)

  expect_equal(nrow(foldc), nrow(exprs(ex_set_na)))
  expect_equal(foldc$Feature_ID, featureNames(ex_set_na))
  expect_true(all(is.na(foldc[1:2, 2:ncol(foldc)])))
})

# P-value correction ----
test_that("P-value correction works", {

  ps <- data.frame(x = letters,
                   x_P = runif(26),
                   P_x = runif(26),
                   y_P = runif(26))

  adj <- adjust_p_values(ps, flags = rep(NA, 26))

  expect_equal(adj$x_P_FDR, p.adjust(ps$x_P, method = "BH"))
  expect_equal(adj[colnames(ps)], ps)
  expect_equal(ncol(adj), ncol(ps) + 2)

  adj2 <- adjust_p_values(ps, flags = c("a", "a", "a",
                                        rep(NA, 23)))

  expect_equal(adj2$x_P_FDR,
               c(rep(NA, 3), p.adjust(ps$x_P[4:26], method = "BH")))

})

# Linear model ----
test_that("Linear model works", {

  cd <- combined_data(drop_qcs(example_set))
  lm_fit <- lm(HILIC_pos_259_9623a4_4322 ~ Time,
               data = cd)
  smry <- summary(lm_fit)

  # Works for a simple example
  lm_res <- perform_lm(drop_qcs(example_set),
                       formula_char = "Feature ~ Time")

  expect_equal(lm_res$Time2_Estimate[1], smry$coefficients[2,1])
  expect_equal(lm_res$Time2_Std_Error[1], smry$coefficients[2,2])
  expect_equal(lm_res$Time2_t_value[1], smry$coefficients[2,3])
  expect_equal(lm_res$Time2_P[1], smry$coefficients[2,4])
  expect_equal(lm_res$R2[1], smry$r.squared)
  expect_equal(lm_res$Adj_R2[1], smry$adj.r.squared)

  # Works with column with only NA values
  ex_set_na <- drop_qcs(mark_nas(example_set, 0))
  exprs(ex_set_na)[1:2, ] <- NA

  lm_res <- perform_lm(ex_set_na,
                       formula_char = "Feature ~ Time")
  expect_equal(nrow(lm_res), nrow(exprs(example_set)))
  expect_equal(lm_res$Feature_ID, featureNames(example_set))
  expect_true(all(is.na(lm_res[1:2, 2:ncol(lm_res)])))

  # FDR correction ignored for flagged compounds
  ex_set_na <- flag_quality(ex_set_na)
  lm_res2 <- perform_lm(ex_set_na, formula_char = "Feature ~ Group")
  flag_idx <- !is.na(flag(ex_set_na))
  expect_true(all(is.na(lm_res2$Group_P_FDR[flag_idx])))

})

# Logistic regression ----
test_that("Logistic regression works", {

  cd <- combined_data(drop_qcs(example_set))
  glm_fit <- glm(Group ~ HILIC_pos_259_9623a4_4322,
               data = cd,
               family = binomial())
  smry <- summary(glm_fit)

  glm_res <- perform_logistic(drop_qcs(example_set),
                       formula_char = "Group ~ Feature")


  expect_equal(glm_res$Feature_Estimate[1], smry$coefficients[2,1])
  expect_equal(glm_res$Feature_Std_Error[1], smry$coefficients[2,2])
  expect_equal(glm_res$Feature_z_value[1], smry$coefficients[2,3])
  expect_equal(glm_res$Feature_P[1], smry$coefficients[2,4])
  expect_equal(glm_res$Feature_LCI95[1], confint(glm_fit)[2,1])
  expect_equal(glm_res$Feature_UCI95[1], confint(glm_fit)[2,2])


  ex_set_na <- drop_qcs(mark_nas(example_set, 0))
  exprs(ex_set_na)[1:2, ] <- NA

  glm_res <- perform_logistic(ex_set_na,
                       formula_char = "Group ~ Feature")
  expect_equal(nrow(glm_res), nrow(exprs(example_set)))
  expect_equal(glm_res$Feature_ID, featureNames(example_set))
  expect_true(all(is.na(glm_res[1:2, 2:ncol(glm_res)])))
})


# Paired t-test ----
test_that("Paired t-test works", {
  ex <- drop_qcs(example_set)
  cd <- combined_data(ex)
  t_res <- perform_paired_t_test(ex, group = "Time", id = "Subject_ID")
  feature <- "HILIC_pos_259_9623a4_4322"
  mean1 <- finite_mean(cd[cd$Time == 1, colnames(cd) == feature])
  mean2 <- finite_mean(cd[cd$Time == 2, colnames(cd) == feature])
  # Check comparison order
  expect_equal(t_res[t_res$Feature_ID == feature, 2], mean1 - mean2)
  # Check row names
  expect_identical(rownames(t_res), featureNames(drop_qcs(example_set)))
  # Check column names
  expect_identical(colnames(t_res), c(
    "Feature_ID",
    "1_vs_2_Estimate",
    "1_vs_2_Lower_CI95",
    "1_vs_2_Upper_CI95",
    "1_vs_2_t_test_P",
    "1_vs_2_t_test_P_FDR"
  ))
})

# Pairwise t-tests ----
test_that("Pairwise t-test works", {
  object <- drop_qcs(example_set)
  pData(object)$Group <- factor(rep(c(rep("A", 3), rep("B", 3), rep("C", 2)), 3))
  pData(object)$Subject_ID <- factor(rep(1:8, 3))
  pData(object)$Time <- factor(c(rep(1, 8), rep(2, 8), rep(3, 8)))

  pwt_res <- perform_pairwise_t_test(object, group = "Time")

  expect_identical(rownames(pwt_res), featureNames(drop_qcs(example_set)))
  prefixes <- c("1_vs_2_", "1_vs_3_", "2_vs_3_")
  suffixes <- c("Estimate", "Lower_CI95", "Upper_CI95", "t_test_P", "t_test_P_FDR")
  cols <- expand.grid(prefixes, suffixes)
  # Check column names
  expect_identical(colnames(pwt_res), c(
    "Feature_ID", "1_Mean", "2_Mean",
    do.call(paste0, cols[order(cols$Var1) & cols$Var1 == prefixes[1], ]),
    "3_Mean",
    do.call(paste0, cols[order(cols$Var1), ])[6:15]
  ))
  # These should be identical as no paired mode
  pData(object)$Subject_ID <- factor(rep(1:12, 2))
  expect_identical(perform_pairwise_t_test(object, group = "Time"), pwt_res)
  # These shouldn't match cause paired mode
  # In this case 4 pairs in each
  expect_failure(expect_identical(perform_pairwise_t_test(object,
                                                          group = "Time",
                                                          id = "Subject_ID",
                                                          is_paired = TRUE
  ), pwt_res))
})

test_that("Pairwise paired t-test works", {
  object <- drop_qcs(example_set)
  pData(object)$Group <- factor(rep(c(rep("A", 3), rep("B", 3), rep("C", 2)), 3))
  pData(object)$Subject_ID <- factor(rep(1:8, 3))
  pData(object)$Time <- factor(c(rep(1, 8), rep(2, 8), rep(3, 8)))

  pwpt_res <- perform_pairwise_t_test(object,
                                      group = "Time",
                                      id = "Subject_ID",
                                      is_paired = TRUE
  )

  expect_identical(rownames(pwpt_res), featureNames(drop_qcs(example_set)))
  prefixes <- c("1_vs_2_", "1_vs_3_", "2_vs_3_")
  suffixes <- c("Estimate", "Lower_CI95", "Upper_CI95", "t_test_P", "t_test_P_FDR")
  cols <- expand.grid(prefixes, suffixes)
  expect_identical(colnames(pwpt_res), c("Feature_ID",
                                        do.call(paste0, cols[order(cols$Var1), ])
  ))
  # Change Subject IDs
  pData(object)$Subject_ID <- factor(rep(1:12, 2))
  pwpt_res_2 <- perform_pairwise_t_test(object,
                                        group = "Time",
                                        id = "Subject_ID",
                                        is_paired = TRUE
  )
  # These shouldn't match because means are counted only from paired samples
  # In this case 4 pairs in each
  expect_failure(expect_identical(pwpt_res_2, pwpt_res))
})


