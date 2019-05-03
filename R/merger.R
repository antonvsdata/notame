
# Check that objects have same special columns
#
# Used to check that the special columns of pheno data parts of MetaboSet objects
# match when merging, called by check_match
#
# @param x,y MetaboSet objects
# @param fun the function to apply, usually one of group_col, time_col or subject_col
check_column_match <- function(x, y, fun) {
  check <- fun(x) == fun(y)
  if (is.na(check)) {
    check <- is.na(fun(x)) & is.na(fun(y))
  }
  if (!check) {
    stop(paste(as.character(substitute(fun)), "returns different column names"))
  }
  if (!is.na(fun(x))) {
    if(!identical(pData(x)[, fun(x)], pData(y)[, fun(y)])) {
      stop(paste(as.character(substitute(fun)), "columns contain different elements"))
    }
  }
}

# Check that two MetaboSet object can be combined
#
# Checks many matching criteria, basically pheno data needs to have similar special columns,
# the amount of samples needs to be the same and feature data and results need to have the
# same columns names. Throws an error if any of the criteria is not fulfilled.
#
# @param x,y MetaboSet objects
check_match <- function(x, y) {
  # Lots of checks to ensure that everything goes smoothly

  # Amount of samples must be equal
  if (nrow(pData(x)) != nrow(pData(y))) {
    stop("Unequal amount of samples")
  }
  # Resulting feature ID must be unique
  feature_id <- c(fData(x)$Feature_ID, fData(y)$Feature_ID)
  if (!all_unique(feature_id)) {
    stop("Merge would result in duplicated feature ID")
  }
  # group_col, time_col, subject_col need to match
  funs <- list(group_col, time_col, subject_col)
  for (i in seq_along(funs)) {
    check_column_match(x, y, funs[[i]])
  }

  if (!identical(pData(x)$Injection_order, pData(y)$Injection_order)) {
    stop("Injection orders are not identical")
  }
  if (!identical(pData(x)$Sample_ID, pData(y)$Sample_ID)) {
    stop("Sample IDs are not identical")
  }


  overlap_cols <- intersect(colnames(pData(x)), colnames(pData(y))) %>%
    setdiff(c("Sample_ID", "Injection_order", group_col(x), time_col(x), subject_col(x)))

  if (length(overlap_cols)) {
    for (overlap_col in overlap_cols) {
      if (!identical(pData(x)[overlap_col], pData(y)[overlap_col])) {
        stop(paste("Columns named", overlap_col, "in pheno data have different content"))
      }
    }
  }

  if (!identical(colnames(fData(x)), colnames(fData(y)))) {
    stop("fData have different column names")
  }

  if (!identical(colnames(exprs(x)), colnames(exprs(y)))) {
    stop("exprs have different column names")
  }

  if (!identical(colnames(results(x)), colnames(results(y)))) {
    stop("results have different column names")
  }

}

# Merge two metaboset objects together
merge_helper <- function(x, y) {
  # Check that the match is ok
  check_match(x,y)

  merged_pdata <- dplyr::left_join(pData(x), pData(y), by =
                                     intersect(colnames(pData(x)), colnames(pData(y)))) %>%
    Biobase::AnnotatedDataFrame()
  rownames(merged_pdata) <- rownames(pData(x))
  merged_exprs <- rbind(exprs(x), exprs(y))
  merged_fdata <- rbind(fData(x), fData(y)) %>%
    Biobase::AnnotatedDataFrame()
  merged_predicted <- rbind(predicted(x), predicted(y))
  merged_results <- rbind(results(x), results(y))

  merged_group_col <- ifelse(!is.na(group_col(x)), group_col(x), group_col(y))
  merged_time_col <- ifelse(!is.na(time_col(x)), time_col(x), time_col(y))
  merged_subject_col <- ifelse(!is.na(subject_col(x)), subject_col(x), subject_col(y))

  merged_object <- MetaboSet(exprs = merged_exprs,
                             phenoData = merged_pdata,
                             featureData = merged_fdata,
                             group_col = merged_group_col,
                             time_col = merged_time_col,
                             subject_col = merged_subject_col,
                             predicted = merged_predicted,
                             results = merged_results)

  merged_object
}

#' Merge MetaboSet objects together
#'
#' @param ... MetaboSet objects or a list of Metaboset objects
#'
#' @return A merged MetaboSet object
#'
#' @export
merge_metabosets <- function(...) {

  # Combine the objects to a list
  objects <- list(...)
  # If a list is given in the first place, it should move to top level
  if (length(objects) == 1) {
    if (class(objects[[1]]) == "list") {
      objects <- objects[[1]]
    }
  }
  # Class check
  if (!all(sapply(objects, class) == "MetaboSet")) {
    stop("The arguments should only contain MetaboSet objects")
  }

  # Merge objects together one by one
  merged <- NULL
  for (object in objects) {
    if (is.null(merged)) {
      merged <- object
    } else {
      merged <- merge_helper(merged, object)
    }
  }

  merged
}
