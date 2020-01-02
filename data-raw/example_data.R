library(dplyr)
devtools::load_all()
# Create a toy dataset for use in examples

# Feature data
set.seed(38)
n_features <- 20
feature_data <- data.frame(Split = "HILIC_pos",
                           Alignment = seq_len(n_features),
                           Mass = runif(n_features, 100, 500),
                           RetentionTime = runif(n_features, 0.5, 8),
                           Column = "HILIC", Mode = "pos",
                           stringsAsFactors = FALSE)

# Create Feature ID
round_mz <- as.numeric(feature_data$Mass) %>% round(digits = 4) %>%
  as.character() %>% gsub("[.]", "_", .)
round_rt <- as.numeric(feature_data$RetentionTime) %>% round(digits = 4) %>%
  as.character() %>% gsub("[.]", "_", .)
feature_data$Feature_ID <- paste0("HILIC_pos_", round_mz, "a", round_rt)
feature_data <- select(feature_data, "Feature_ID", everything())

rownames(feature_data) <- feature_data$Feature_ID

# Pheno data
n_samples <- 30

qc_idx <- seq(1, n_samples/2, length.out = 3) %>% round()
subject_ids <- as.character(seq_len(n_samples/2))
group <- sample(LETTERS[1:2], n_samples/2, replace = TRUE)
time <- as.character(rep(c(1,2), each = n_samples/2))
group[qc_idx] <- "QC"
subject_ids[qc_idx] <- "QC"
time[c(qc_idx, qc_idx + n_samples/2)] <- "QC"
qc <- ifelse(group == "QC", "QC", "Sample")
pheno_data <- data.frame(Injection_order = seq_len(n_samples),
                         Sample_ID = paste0("Demo_", seq_len(n_samples)),
                         Subject_ID = subject_ids,
                         Group = factor(group), QC = factor(qc),
                         Time = factor(time))

rownames(pheno_data) <- pheno_data$Sample_ID



# Assay data

# Random means for each feature
means <- runif(n_features, 3000, 33000)

# Normally distributed data around the mean
assay_data <- t(sapply(means, function(x) {rnorm(n_samples, x, 0.3*x)}))
# Add drift effect to the data
coefs <- runif(n_samples, 0.4, 0.9) * sample(c(-1, 1), n_samples, replace = TRUE)

# Randomly choose linear or logarithmic trend
for (i in seq_len(nrow(assay_data))) {
  if (rnorm(1) > 0) {
    assay_data[i,] <- assay_data[i,] + means[i] * coefs[i] * log(pheno_data$Injection_order)
  } else {
    assay_data[i,] <- assay_data[i,] + means[i] * coefs[i] * 0.08* pheno_data$Injection_order
  }
}
# Set random indexes to zero (for missing values)
n_missing <- 30
row_zeros <- sample(seq_len(nrow(assay_data)), n_missing, replace = TRUE)
col_zeros <- sample(seq_len(ncol(assay_data)), n_missing, replace = TRUE)
for (i in seq_len(n_missing)) {
  assay_data[row_zeros[i], col_zeros[i]] <- 0
}

assay_data <- abs(assay_data)

# Set dimension names
rownames(assay_data) <- rownames(feature_data)
colnames(assay_data) <- rownames(pheno_data)

# Construct object
example_set <- construct_MetaboSet(exprs = assay_data, pheno_data = pheno_data, feature_data = feature_data,
                                   group_col = "Group", time_col = "Time", subject_col = "Subject_ID")[[1]]



usethis::use_data(example_set, overwrite = TRUE)
