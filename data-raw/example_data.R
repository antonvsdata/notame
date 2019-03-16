library(dplyr)
library(Biobase)

# Create a toy dataset for use in examples

# Feature data
set.seed(38)
n_features <- 20
feature_data <- data.frame(Split = "HILIC_pos",
                           Alignment = seq_len(n_features),
                           Mass = runif(n_features, 100, 500),
                           RetentionTime = runif(n_features, 0.5, 8),
                           Column = "HILIC", Mode = "pos",
                           Flag = NA_character_)

# Create Feature ID
round_mz <- as.numeric(feature_data$Mass) %>% round(digits = 4) %>%
  as.character() %>% gsub("[.]", "_", .)
round_rt <- as.numeric(feature_data$RetentionTime) %>% round(digits = 4) %>%
  as.character() %>% gsub("[.]", "_", .)
feature_data$Feature_ID <- paste0("HILIC_neg_", round_mz, "a", round_rt)
feature_data <- select(feature_data, "Feature_ID", everything())

rownames(feature_data) <- feature_data$Feature_ID

# Pheno data

n_samples <- 30

qc_idx <- seq(1, n_samples, length.out = 5) %>% round()
group <- sample(LETTERS[1:2], n_samples, replace = TRUE)
group[qc_idx] <- "QC"
qc <- ifelse(seq_len(n_samples) %in% qc_idx, "QC", "Sample")
pheno_data <- data.frame(Injection_order = seq_len(n_samples),
                         Sample_ID = paste0("Demo_", seq_len(n_samples)),
                         Group = group, QC = qc,
                         Time = sample(seq_len(2), n_samples, replace = TRUE))

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

rownames(assay_data) <- rownames(feature_data)
colnames(assay_data) <- rownames(pheno_data)

pheno_data <- AnnotatedDataFrame(data = pheno_data)
feature_data <- AnnotatedDataFrame(data = feature_data)

example_set <- MetaboSet(exprs = assay_data,
                         phenoData = pheno_data,
                         featureData = feature_data,
                         group_col = "Group",
                         time_col = "Time",
                         subject_col = NA_character_,
                         predicted = matrix(NA_real_, nrow = nrow(assay_data),
                                            ncol = ncol(assay_data),
                                            dimnames = dimnames(assay_data)))



usethis::use_data(example_set)
