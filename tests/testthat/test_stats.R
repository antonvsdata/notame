context("Testing statistics")

library(notame)
# Summary statistics ----
test_that("summary statistics work without grouping", {
  # Force foreach to count sequentially
  foreach::registerDoSEQ()
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
                   B_vs_A_2_minus_1_Cohen_d = d,
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

test_that("Cohen's d data checking works", {
  ex <- example_set
  ex$Group <- c(1, 2)
  expect_error(cohens_d(ex, id = "Subject_ID", time = "Time"), "should be a factor")

  ex <- example_set
  ex$Time <- c(1,2)
  expect_error(cohens_d(ex, id = "Subject_ID", time = "Time"), "should be a factor")

  ex <- example_set
  ex$Group <- factor(1)
  expect_error(cohens_d(ex, id = "Subject_ID", time = "Time"), "should have at least two levels")

  ex <- example_set
  ex$Time <- factor(1)
  expect_error(cohens_d(ex, id = "Subject_ID", time = "Time"), "should have at least two levels")
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
                   B_vs_A_FC = 1,
                   QC_vs_A_FC = 1,
                   QC_vs_B_FC = 1,
                   stringsAsFactors = FALSE)
  rownames(fc) <- fc$Feature_ID
  for (i in seq_len(nrow(fc))) {
    feature <- fc$Feature_ID[i]
    mean_a <- finite_mean(cd[cd$Group == "A", feature])
    mean_b <- finite_mean(cd[cd$Group == "B", feature])
    mean_qc <- finite_mean(cd[cd$Group == "QC", feature])

    fc$B_vs_A_FC[i] <- mean_b/mean_a
    fc$QC_vs_A_FC[i] <- mean_qc/mean_a
    fc$QC_vs_B_FC[i] <- mean_qc/mean_b
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

test_that("Cohens D values are counted right", {

  object <- drop_qcs(example_set)
  pData(object)$Group <- factor(c("A", "B", "C"))

  data <- combined_data(object)
  features <- featureNames(object)
  group1 <- data[which(data[, "Group"] == levels(data[,"Group"])[1]), ]
  group2 <- data[which(data[, "Group"] == levels(data[,"Group"])[2]), ]
  group3 <- data[which(data[, "Group"] == levels(data[,"Group"])[3]), ]
  ds <-  foreach::foreach(i = seq_along(features), .combine = rbind) %do% {
    feature <- features[i]
    f1 <- group1[, feature]
    f2 <- group2[, feature]
    f3 <- group3[, feature]
    d <- data.frame(Feature_ID = feature,
                    B_vs_A_Cohen_d = (finite_mean(f2) - finite_mean(f1)) /
                      sqrt((finite_sd(f1)^2 + finite_sd(f2)^2) / 2),
                    C_vs_A_Cohen_d = (finite_mean(f3) - finite_mean(f1)) /
                      sqrt((finite_sd(f1)^2 + finite_sd(f3)^2) / 2),
                    C_vs_B_Cohen_d = (finite_mean(f3) - finite_mean(f2)) /
                      sqrt((finite_sd(f3)^2 + finite_sd(f2)^2) / 2),
                    stringsAsFactors = FALSE)
  }
  rownames(ds) <- ds$Feature_ID
  foreach::registerDoSEQ()
  cohd <- cohens_d(object)
  expect_identical(cohd, ds)
})

test_that("Cohens D values between time points are counted right", {

  object <- drop_qcs(example_set)
  pData(object)$Group <- factor(rep(c(rep("A", 3), rep("B", 3), rep("C", 2)), 3))
  pData(object)$Subject_ID <- factor(rep(1:8, 3))
  pData(object)$Time <- factor(c(rep(1, 8), rep(2, 8), rep(3, 8)))

  # Create results with time points manually
  data <- combined_data(object)
  features <- featureNames(object)

  group_combos <- combn(levels(pData(object)[, "Group"]), 2)
  time_combos <- combn(levels(pData(object)[, "Time"]), 2)
  for (i in seq_len(ncol(group_combos))) {
    group1 <- data[which(data[, "Group"] == group_combos[1, i]), ]
    group2 <- data[which(data[, "Group"] == group_combos[2, i]), ]
    for (j in seq_len(ncol(time_combos))) {
      time1 <- data[which(data[, "Time"] == time_combos[1, i]), ]
      time2 <- data[which(data[, "Time"] == time_combos[2, i]), ]
      common_ids <- intersect(time1[, "Subject_ID"], time2[, "Subject_ID"])
      new_data <- time2[time2[, "Subject_ID"] %in% common_ids, features] -
        time1[time1[, "Subject_ID"] %in% common_ids, features]
      # Split to groups
      group1 <- new_data[which(time1[, "Group"] == levels(time1[,"Group"])[1]), ]
      group2 <- new_data[which(time1[, "Group"] == levels(time1[,"Group"])[2]), ]
    }
  }

  group1 <- data[which(data[, "Group"] == levels(data[,"Group"])[1]), ]
  group2 <- data[which(data[, "Group"] == levels(data[,"Group"])[2]), ]
  group3 <- data[which(data[, "Group"] == levels(data[,"Group"])[3]), ]
  time1 <- data[which(data[, "Time"] == levels(data[, "Time"])[1]), ]
  time2 <- data[which(data[, "Time"] == levels(data[, "Time"])[2]), ]
  time3 <- data[which(data[, "Time"] == levels(data[, "Time"])[3]), ]
  common_ids_21 <- intersect(time1[, "Subject_ID"], time2[, "Subject_ID"])
  common_ids_31 <- intersect(time1[, "Subject_ID"], time3[, "Subject_ID"])
  common_ids_32 <- intersect(time3[, "Subject_ID"], time2[, "Subject_ID"])
  # Change between time points
  new_data_21 <- time2[time2$Subject_ID %in% common_ids_21, features] -
    time1[time1$Subject_ID %in% common_ids_21, features]
  new_data_31 <- time3[time3$Subject_ID %in% common_ids_31, features] -
    time1[time1$Subject_ID %in% common_ids_31, features]
  new_data_32 <- time3[time3$Subject_ID %in% common_ids_32, features] -
    time2[time2$Subject_ID %in% common_ids_32, features]
  new_data <- rbind(new_data_21, new_data_31, new_data_32)
  # Split to groups
  group1 <- new_data_21[which(time1[, "Group"] == levels(time1[,"Group"])[1]), ] # A 2-1
  group2 <- new_data_21[which(time1[, "Group"] == levels(time1[,"Group"])[2]), ] # B 2-1
  group3 <- new_data_21[which(time1[, "Group"] == levels(time1[,"Group"])[3]), ] # C 2-1
  group4 <- new_data_31[which(time2[, "Group"] == levels(time2[,"Group"])[1]), ] # A 3-1
  group5 <- new_data_31[which(time2[, "Group"] == levels(time2[,"Group"])[2]), ] # B 3-1
  group6 <- new_data_31[which(time2[, "Group"] == levels(time2[,"Group"])[3]), ] # C 3-1
  group7 <- new_data_32[which(time3[, "Group"] == levels(time3[,"Group"])[1]), ] # A 3-2
  group8 <- new_data_32[which(time3[, "Group"] == levels(time3[,"Group"])[2]), ] # B 3-2
  group9 <- new_data_32[which(time3[, "Group"] == levels(time3[,"Group"])[3]), ] # C 3-2
  ds <-  foreach::foreach(i = seq_along(features), .combine = rbind) %do% {
    feature <- features[i]
    f1 <- group1[, feature] # A 2-1
    f2 <- group2[, feature] # B 2-1
    f3 <- group3[, feature] # C 2-1
    f4 <- group4[, feature] # A 3-1
    f5 <- group5[, feature] # B 3-1
    f6 <- group6[, feature] # C 3-1
    f7 <- group7[, feature] # A 3-2
    f8 <- group8[, feature] # B 3-2
    f9 <- group9[, feature] # C 3-2
    d <- data.frame(Feature_ID = feature,
                    "B_vs_A_2_minus_1_Cohen_d" = (finite_mean(f2) - finite_mean(f1)) /
                      sqrt((finite_sd(f1)^2 + finite_sd(f2)^2) / 2),
                    "B_vs_A_3_minus_1_Cohen_d" = (finite_mean(f5) - finite_mean(f4)) /
                      sqrt((finite_sd(f4)^2 + finite_sd(f5)^2) / 2),
                    "B_vs_A_3_minus_2_Cohen_d" = (finite_mean(f8) - finite_mean(f7)) /
                      sqrt((finite_sd(f7)^2 + finite_sd(f8)^2) / 2),
                    "C_vs_A_2_minus_1_Cohen_d" = (finite_mean(f3) - finite_mean(f1)) /
                      sqrt((finite_sd(f1)^2 + finite_sd(f3)^2) / 2),
                    "C_vs_A_3_minus_1_Cohen_d" = (finite_mean(f6) - finite_mean(f4)) /
                      sqrt((finite_sd(f4)^2 + finite_sd(f6)^2) / 2),
                    "C_vs_A_3_minus_2_Cohen_d" = (finite_mean(f9) - finite_mean(f7)) /
                      sqrt((finite_sd(f7)^2 + finite_sd(f9)^2) / 2),
                    "C_vs_B_2_minus_1_Cohen_d" = (finite_mean(f3) - finite_mean(f2)) /
                      sqrt((finite_sd(f2)^2 + finite_sd(f3)^2) / 2),
                    "C_vs_B_3_minus_1_Cohen_d" = (finite_mean(f6) - finite_mean(f5)) /
                      sqrt((finite_sd(f5)^2 + finite_sd(f6)^2) / 2),
                    "C_vs_B_3_minus_2_Cohen_d" = (finite_mean(f9) - finite_mean(f8)) /
                      sqrt((finite_sd(f8)^2 + finite_sd(f9)^2) / 2),
                    stringsAsFactors = FALSE)
  }
  rownames(ds) <- ds$Feature_ID
  foreach::registerDoSEQ()
  expect_identical(cohens_d(object, time = "Time", id = "Subject_ID"), ds)
})

test_that("Cohens D warnings work", {

  object <- drop_qcs(example_set)
  pData(object)$Group <- factor(rep(c(rep("A", 3), rep("B", 3), rep("C", 2)), 3))
  pData(object)$Subject_ID <- factor(rep(1:8, 3))
  pData(object)$Time <- factor(c(rep(1, 8), rep(2, 8), rep(3, 8)))
  # Remove one sample
  expect_warning(cohens_d(object[, -1],
                          time = "Time",
                          id = "Subject_ID",
                          ),
                 regexp = "[One or more subject(s) missing time points]"
  )
  # all time points missing in one group
  expect_warning(cohens_d(object[, -(23:24)], time = "Time", id = "Subject_ID"),
                 regexp = "[Groups don't have two observations of at least two subjects]"
  )
  # all but one time points missing in one group
  expect_warning(cohens_d(object[, -24], time = "Time", id = "Subject_ID"),
                 regexp = "[Groups don't have two observations of at least two subjects]"
  )
  # Same subject is in two groups
  tmp <- object
  pData(tmp)$Group[1] <- "B"
  expect_warning(cohens_d(tmp, time = "Time", id = "Subject_ID"),
                 regexp = "[Same subject recorded in two groups]"
  )
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


