
# Check that objects have same special columns
#
# Used to check that the special columns of pheno data parts of MetaboSet objects
# match when merging, called by check_match
#
# @param x,y MetaboSet objects
# @param fun the function to apply, usually one of group_col, time_col or subject_col
check_column_match <- function(x, y, fun, name) {
  check <- fun(x) == fun(y)
  if (is.na(check)) {
    check <- is.na(fun(x)) & is.na(fun(y))
  }
  if (!check) {
    stop(paste(name, "returns different column names"))
  }
  common <- intersect(sampleNames(x), sampleNames(y))
  if (!is.na(fun(x))) {
    if(!identical(pData(x)[common, fun(x)], pData(y)[common, fun(y)])) {
      stop(paste(name, "columns contain different elements for common samples"))
    }
  }
}

# Check that two MetaboSet object can be combined
#
# Checks many matching criteria, basically pheno data needs to have similar special columns,
# the number of samples needs to be the same and feature data need to have the
# same columns names. Throws an error if any of the criteria is not fulfilled.
#
# @param x,y MetaboSet objects
check_match <- function(x, y) {
  # Lots of checks to ensure that everything goes smoothly

  # Amount of samples must be equal
  if (nrow(pData(x)) != nrow(pData(y))) {
    warning("Unequal amount of samples")
  }
  # Resulting feature ID must be unique
  feature_id <- c(fData(x)$Feature_ID, fData(y)$Feature_ID)
  if (!all_unique(feature_id)) {
    stop("Merge would result in duplicated feature ID")
  }
  # group_col, time_col, subject_col need to match
  funs <- list("group_col" = group_col, "time_col" = time_col, "subject_col" = subject_col)
  for (i in seq_along(funs)) {
    check_column_match(x, y, funs[[i]], names(funs)[i])
  }
  common <- intersect(sampleNames(x), sampleNames(y))
  if (!identical(pData(x)[common, "Injection_order"], pData(y)[common, "Injection_order"])) {
    stop("Injection orders of common samples are not identical")
  }
  if (!identical(pData(x)$Sample_ID, pData(y)$Sample_ID)) {
    warning("Sample IDs are not identical")
    samples_x <- setdiff(sampleNames(x), sampleNames(y))
    samples_y <- setdiff(sampleNames(y), sampleNames(x))
    log_text("Merging objects with unequal amounts of samples.")
    log_text("Samples only in first object:")
    log_text(paste0(paste(samples_x, collapse = ", "), "\n"))
    log_text("Samples only in second object:")
    log_text(paste0(paste(samples_y, collapse = ", "), "\n"))
  }


  overlap_cols <- intersect(colnames(pData(x)), colnames(pData(y))) %>%
    setdiff(c("Sample_ID", "Injection_order", group_col(x), time_col(x), subject_col(x)))

  if (length(overlap_cols)) {
    for (overlap_col in overlap_cols) {
      if (!identical(pData(x)[common, overlap_col], pData(y)[common, overlap_col])) {
        stop(paste("Columns named", overlap_col, "in pheno data have different content"))
      }
    }
  }

  if (!identical(colnames(fData(x)), colnames(fData(y)))) {
    stop("fData have different column names")
  }

}

# Merge two MetaboSet objects together
merge_mode_helper <- function(x, y) {
  # Create dummy injection order if original ones differ
  common <- intersect(sampleNames(x), sampleNames(y))
  if (!identical(pData(x)[common, "Injection_order"], pData(y)[common, "Injection_order"])) {
    log_text("Injection order differs between modes. Creating dummy injection order")
    x_modes <- unique(fData(x)$Split)
    # Save original injection order for first mode
    if (length(x_modes) == 1) {
      pData(x)[, paste0(x_modes, "_Injection_order")] <- x$Injection_order
    }
    dummy_injection <- as.numeric(-seq_along(pData(x)$Sample_ID))
    names(dummy_injection) <- x$Sample_ID
    pData(x)$Injection_order <- dummy_injection
    # Save original injection order in other modes
    pData(y)[, paste0(unique(fData(y)$Split), "_Injection_order")] <- y$Injection_order
    # Update dummy injection
    y_in_x <- y$Sample_ID %in% x$Sample_ID
    new_io <- seq(from = min(dummy_injection) - 1, length.out = sum(!y_in_x))
    names(new_io) <- y$Sample_ID[!y_in_x]
    dummy_injection <- append(dummy_injection, new_io)

    pData(y)$Injection_order <- dummy_injection[match(pData(y)$Sample_ID, names(dummy_injection))]
    log_text("Dummy injection order (row numbers) created")
  }
  # Check that the match is ok
  check_match(x, y)

  merged_pdata <- dplyr::full_join(pData(x), pData(y),
                                   by = intersect(colnames(pData(x)),
                                                  colnames(pData(y))
                                   )
  ) %>%
    Biobase::AnnotatedDataFrame()
  rownames(merged_pdata) <- merged_pdata$Sample_ID
  merged_fdata <- rbind(fData(x), fData(y)) %>%
    Biobase::AnnotatedDataFrame()
  if (identical(colnames(exprs(x)), colnames(exprs(y)))) {
    merged_exprs <- rbind(exprs(x), exprs(y))
  } else {
    merged_exprs <- dplyr::bind_rows(as.data.frame(exprs(x)), as.data.frame(exprs(y))) %>% as.matrix()
    rownames(merged_exprs) <- rownames(merged_fdata)
  }

  merged_group_col <- ifelse(!is.na(group_col(x)), group_col(x), group_col(y))
  merged_time_col <- ifelse(!is.na(time_col(x)), time_col(x), time_col(y))
  merged_subject_col <- ifelse(!is.na(subject_col(x)), subject_col(x), subject_col(y))

  merged_object <- MetaboSet(exprs = merged_exprs,
                             phenoData = merged_pdata,
                             featureData = merged_fdata,
                             group_col = merged_group_col,
                             time_col = merged_time_col,
                             subject_col = merged_subject_col)

  merged_object
}

# Convert metaboset objects in ... to al list
# OR if a list is given in the first place, preserve that list
to_list <- function(...) {
  # Combine the objects to a list
  objects <- list(...)
  # If a list is given in the first place, it should move to top level
  if (length(objects) == 1) {
    if (class(objects[[1]]) == "list") {
      objects <- objects[[1]]
    }
  }

  objects
}

#' Merge MetaboSet objects together
#'
#' Merges two or more MetaboSet objects together. Can be used to merge analytical modes or batches.
#'
#'
#' @param ... MetaboSet objects or a list of Metaboset objects
#' @param merge what to merge? features is used for combining analytical modes,
#' samples is used for batches
#'
#' @return A merged MetaboSet object
#'
#' @details When merging samples, sample IDs that beging with "QC" or "Ref" are combined so that they have
#' running numbers on them. This means that if both bathces have samples called "QC_1", this will not result in an error,
#' but the sample IDs will be adjusted so that they are unique
#'
#' @examples
#' # Merge analytical modes
#' merged <- merge_metabosets(hilic_neg_sample, hilic_pos_sample,
#'                            rp_neg_sample, rp_pos_sample)
#' # Merge batches
#' batch1 <- merged_sample[, merged_sample$Batch == 1]
#' batch2 <- merged_sample[, merged_sample$Batch == 2]
#' merged <- merge_batches(batch1, batch2)
#'
#' @export
merge_metabosets <- function(..., merge = c("features", "samples")) {

  merge <- match.arg(merge)
  # Combine the objects to a list
  objects <- to_list(...)
  # Class check
  if (!all(sapply(objects, class) == "MetaboSet")) {
    stop("The arguments should only contain MetaboSet objects")
  }
  # Choose merging function
  if (merge == "features") {
    merge_fun <- merge_mode_helper
  } else {
    merge_fun <- merge_batch_helper
  }
  # Merge objects together one by one
  merged <- NULL
  for (object in objects) {
    if (is.null(merged)) {
      merged <- object
    } else {
      merged <- merge_fun(merged, object)
    }
  }
  merged
}


fdata_batch_helper <- function(fx, fy) {

  non_identical_cols <- !identical(colnames(fx), colnames(fy))
  if (non_identical_cols) {
    only_x_cols <- setdiff(colnames(fx), colnames(fy))
    only_x <- fx[c("Feature_ID", only_x_cols)]
    fx[only_x_cols] <- NULL
    only_y_cols <- setdiff(colnames(fy), colnames(fx))
    only_y <- fy[c("Feature_ID", only_y_cols)]
    fy[only_y_cols] <- NULL
  }

  # Combine common features: all NAs in fx are replaced by a value from fy
  common_features <- intersect(fx$Feature_ID, fy$Feature_ID)
  for (cf in common_features) {
    na_idx <- is.na(fx[cf, ])
    fx[cf, na_idx] <- fy[cf, na_idx]
  }
  new_features <- setdiff(fy$Feature_ID, fx$Feature_ID)

  merged_fdata <- rbind(fx, fy[new_features, ])

  if (non_identical_cols) {
    merged_fdata <- dplyr::left_join(merged_fdata, only_x, by = "Feature_ID") %>%
      dplyr::left_join(only_y, by = "Feature_ID")
    rownames(merged_fdata) <- merged_fdata$Feature_ID
  }
  merged_fdata
}


merge_batch_helper <- function(x, y) {

  merged_pdata <- rbind(pData(x), pData(y))
  if (anyDuplicated(merged_pdata$Sample_ID)) {
    log_text("Found duplicated sample IDs when merging, renaming QC and Ref samples")
    merged_pdata$Sample_ID[grepl("QC", merged_pdata$Sample_ID)] <- paste0("QC_", seq_len(sum(grepl("QC", merged_pdata$Sample_ID))))
    merged_pdata$Sample_ID[grepl("Ref", merged_pdata$Sample_ID)] <- paste0("Ref_", seq_len(sum(grepl("Ref", merged_pdata$Sample_ID))))
  }

  rownames(merged_pdata) <- merged_pdata$Sample_ID
  merged_pdata <- merged_pdata %>%
    Biobase::AnnotatedDataFrame()

  if (identical(rownames(exprs(x)), rownames(exprs(y)))) {
    merged_exprs <- cbind(exprs(x), exprs(y))
    colnames(merged_exprs) <- rownames(merged_pdata)
  } else {
    merged_exprs <- dplyr::bind_rows(as.data.frame(t(exprs(x))), as.data.frame(t(exprs(y)))) %>% t()
    colnames(merged_exprs) <- rownames(merged_pdata)
  }

  merged_fdata <- fdata_batch_helper(fData(x), fData(y)) %>%
    Biobase::AnnotatedDataFrame()
  merged_group_col <- ifelse(!is.na(group_col(x)), group_col(x), group_col(y))
  merged_time_col <- ifelse(!is.na(time_col(x)), time_col(x), time_col(y))
  merged_subject_col <- ifelse(!is.na(subject_col(x)), subject_col(x), subject_col(y))

  merged_object <- MetaboSet(exprs = merged_exprs,
                             phenoData = merged_pdata,
                             featureData = merged_fdata,
                             group_col = merged_group_col,
                             time_col = merged_time_col,
                             subject_col = merged_subject_col)

  merged_object
}

#' DEPRECATED: Merge MetaboSet objects of batches together
#'
#' NOTE: This function is deprecated, use merge_metabosets(..., merge = "samples") instead.
#'
#' @param ... MetaboSet objects or a list of Metaboset objects
#'
#' @return A merged MetaboSet object
#'
#' @examples
#' batch1 <- merged_sample[, merged_sample$Batch == 1]
#' batch2 <- merged_sample[, merged_sample$Batch == 2]
#' merged <- merge_batches(batch1, batch2)
#'
#' @export
merge_batches <- function(...) {

  warning("merge_batches is deprecated, merge_metabosets can now merge object from different batches as
          well as objects from different modes. merge_batches(...) is equivalent to merge_metabosets(..., merge = 'samples')", .call = FALSE)

  merge_metabosets(..., merge = "samples")
}
