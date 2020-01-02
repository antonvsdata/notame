context("Testing statistics")

library(amp)

test_that("summary statistics work without grouping", {

  smry <- summary_statistics(mark_nas(example_set, 0),
                             grouping_cols = NA)
  ex <- exprs(mark_nas(example_set, 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(ex, 1, fun)), smry[, gsub("finite_", "", fun)])
  }
  expect_equal(unname(apply(ex, 1, finite_quantile, probs = 0.25)), smry$Q25)
  expect_equal(unname(apply(ex, 1, finite_quantile, probs = 0.75)), smry$Q75)
})


test_that("summary statistics work with grouping", {

  smry <- summary_statistics(mark_nas(example_set, 0))

  exa <- exprs(mark_nas(example_set[, example_set$Group == "A"], 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(exa, 1, fun)), smry[, gsub("finite", "Group_A", fun)])
  }
  expect_equal(unname(apply(exa, 1, finite_quantile, probs = 0.25)), smry$Group_A_Q25)
  expect_equal(unname(apply(exa, 1, finite_quantile, probs = 0.75)), smry$Group_A_Q75)

  exb <- exprs(mark_nas(example_set[, example_set$Group == "B"], 0))

  for (fun in c("finite_mean", "finite_sd", "finite_median", "finite_mad")) {
    expect_equal(unname(apply(exb, 1, fun)), smry[, gsub("finite", "Group_B", fun)])
  }
  expect_equal(unname(apply(exb, 1, finite_quantile, probs = 0.25)), smry$Group_B_Q25)
  expect_equal(unname(apply(exb, 1, finite_quantile, probs = 0.75)), smry$Group_B_Q75)
})

test_that("summary statistics work with all NA features", {

  ex_set_na <- mark_nas(example_set, 0)
  exprs(ex_set_na)[1, ] <- NA

  smry <- summary_statistics(ex_set_na)

  expect_equal(nrow(smry), nrow(exprs(ex_set_na)))
  expect_equal(smry$Feature_ID, featureNames(ex_set_na))
  expect_true(all(is.na(smry[1, 2:ncol(smry)])))
})


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
    d <- c(d, (mean_b - mean_a)/mean(c(sd_a, sd_b)))
  }

  cohd <- cohens_d(ex)

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
  ex$Group <- c(1,2)
  expect_error(cohens_d(ex), "Column Time")

  ex$Time <- c(1,2)
  ex$Group <- 1
  expect_error(cohens_d(ex), "Column Group")
})


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



