#

#' Extract quality information of features
#'
#' @param object a MetaboSet object
#'
#' @export
quality <- function(object) {
  if (!all(c("RSD", "RSD_r", "D_ratio", "D_ratio_r") %in% colnames(fData(object)))) {
    return(NULL)
  }
  fData(object)[c(
    "Feature_ID", "RSD", "RSD_r", "D_ratio",
    "D_ratio_r"
  )]
}

erase_quality <- function(object) {
  if (!all(c("RSD", "RSD_r", "D_ratio", "D_ratio_r") %in% colnames(fData(object)))) {
    return(NULL)
  }
  fData(object)[c("RSD", "RSD_r", "D_ratio", "D_ratio_r")] <- NULL
  object
}


#' Assess quality information of features
#'
#' @param object a MetaboSet object
#'
#' @export
assess_quality <- function(object) {
  # Remove old quality metrics
  if (!is.null(quality(object))) {
    object <- erase_quality(object)
  }

  qc_data <- exprs(object)[, object$QC == "QC"]
  sample_data <- exprs(object)[, object$QC != "QC"]

  quality_metrics <- foreach::foreach(
    i = seq_len(nrow(sample_data)), .combine = rbind,
    .export = c("finite_sd", "finite_mad", "finite_mean", "finite_median")
  ) %dopar% {
    data.frame(
      Feature_ID = rownames(sample_data)[i],
      RSD = finite_sd(qc_data[i, ]) / abs(finite_mean(qc_data[i, ])),
      RSD_r = finite_mad(qc_data[i, ]) / abs(finite_median(qc_data[i, ])),
      D_ratio = finite_sd(qc_data[i, ]) / finite_sd(sample_data[i, ]),
      D_ratio_r = finite_mad(qc_data[i, ]) / finite_mad(sample_data[i, ]),
      row.names = rownames(sample_data)[i], stringsAsFactors = FALSE
    )
  }

  object <- join_fData(object, quality_metrics)

  object
}


#' Flag low-quality features
#'
#' Flags low-quality features using the quality metrics defined in (Broadhurst 2018). The metrics are described in
#' more detain in Details. A condition for keeping the features is given as a character,
#' which is passed to \code{dplyr::filter}.
#'
#' @param object a MetaboSet object
#' @param condition character, condition for keeping the features, see Details
#'
#' @details The quality metrics measure two things: internal spread of the QCs,
#' and spread of the QCs compared to the spread of the biological samples.
#'   Internal spread is measured with relative standard deviation (RSD), also known as coefficient of variation (CV).
#'   \deqn{RSD = sd(QC) / mean(QC) }
#'   Where \eqn{sd(QC)} is the standard deviation of the QC samples and \eqn{mean(QC)} is
#'   the sample mean of the signal in the QC samples.
#'   RSD can also be replaced by a non-parametric, robust version based on the median and
#'   median absolute deviation (MAD):
#'   \deqn{RSD_r = 1.4826 * MAD(QC) / median(QC)}
#'   The spread of the QC samples compared to the biological samples is measured using a metric called D-ratio:
#'   \deqn{D_ratio = sd(QC) / sd(biological)}
#'   Or, as before, a non-parametric, robust alternative:
#'   \deqn{D_ratio_r = MAD(QC) / MAD(biolofical) }
#' The default condition keeps features that pass either of the two following conditions:
#' \deqn{RSD_r < 0.2 & D_ratio_r < 0.4}
#' \deqn{RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1}
#'
#' @return a MetaboSet object with the features flagged
#'
#' @references Broadhurst, David et al. Guidelines and considerations for the use of system suitability
#' and quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies.
#' Metabolomics : Official journal of the Metabolomic Society vol. 14,6 (2018): 72. doi:10.1007/s11306-018-1367-3
#'
#' @examples
#' ex_set <- flag_quality(merged_sample)
#' fData(ex_set)
#' # Custom condition
#' ex_set <- flag_quality(merged_sample, condition = "RSD_r < 0.3 & D_ratio_r < 0.6")
#' fData(ex_set)
#'
#' @export
flag_quality <- function(object,
                         condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) | (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") {
  if (is.null(quality(object))) {
    object <- assess_quality(object)
  }
  add_citation(
    "Quality metrics were computed according to guidelines in:",
    "Broadhurst, David et al. Guidelines and considerations for the use of system suitability and
    quality control samples in mass spectrometry assays applied in untargeted clinical metabolomic studies.
    Metabolomics : Official journal of the Metabolomic Society vol. 14,6 (2018): 72. doi:10.1007/s11306-018-1367-3"
  )
  good <- paste0("fData(object) %>% dplyr::filter(", condition, ")") %>%
    parse(text = .) %>%
    eval()
  good <- good$Feature_ID

  idx <- is.na(flag(object)) & !fData(object)$Feature_ID %in% good
  flag(object)[idx] <- "Low_quality"

  percentage <- scales::percent(sum(flag(object) == "Low_quality", na.rm = TRUE) / nrow(fData(object)))
  log_text(paste0("\n", percentage, " of features flagged for low quality"))

  object
}

#' Flag features with low detection rate
#'
#' Flags features with too high amount of missing values. There are two detection rate limits, both defined as the
#' minimum proportion of samples that need to have a value (not NA) for the feature to be kept. \code{qc_limit} is
#' the detection rate limit for QC samples, \code{group_limit} is the detection rate limit for the actual study groups.
#' If the group limit is passed for AT LEAST ONE GROUP, then the feature is kept. Features with low detection rate
#' in QCs are flagged as "Low_qc_detection", while low detection rate in the study groups is flagged as
#' "Low_group_detection". The detection rates for all the groups are recorded in \code{fData(object)}.
#'
#' @param object a MetaboSet object
#' @param qc_limit the detection rate limit for QC samples
#' @param group_limit the detection rate limit for study groups
#' @param group the columns name in sample information to use as the grouping variable
#'
#' @return a MetaboSet object with the features flagged
#'
#' @examples
#' ex_set <- flag_detection(merged_sample)
#' fData(ex_set)
#'
#' @export
flag_detection <- function(object, qc_limit = 0.7, group_limit = 0.5, group = group_col(object)) {
  found_qc <- Biobase::esApply(object[, object$QC == "QC"], 1, prop_found)
  bad_qc <- names(which(found_qc < qc_limit))

  found_qc_df <- data.frame(
    Feature_ID = names(found_qc),
    Detection_rate_QC = found_qc,
    stringsAsFactors = FALSE
  )

  idx <- is.na(flag(object)) & fData(object)$Feature_ID %in% bad_qc
  flag(object)[idx] <- "Low_qc_detection"

  # Compute proportions found in each study group
  if (!is.na(group)) {
    proportions <- combined_data(object)[, c("Sample_ID", group, featureNames(object))] %>%
      tidyr::gather(Feature_ID, Intensity, featureNames(object)) %>%
      dplyr::group_by(Feature_ID, !!as.name(group)) %>%
      dplyr::summarise(proportion_found = prop_found(Intensity)) %>%
      tidyr::spread(!!as.name(group), "proportion_found")
    # Remove a possible QC column
    proportions$QC <- NULL
    colnames(proportions)[-1] <- paste0("Detection_rate_", group, colnames(proportions)[-1])
    # Check if any group has enough non-missing entries
    proportions$good <- apply(proportions[-1], 1, function(x) {
      any(x >= group_limit)
    })

    idx <- is.na(flag(object)) & (!fData(object)$Feature_ID %in% proportions$Feature_ID[proportions$good])
    flag(object)[idx] <- "Low_group_detection"
    # Add detection rates to feature data
    proportions <- dplyr::left_join(proportions, found_qc_df, by = "Feature_ID")
    proportions$good <- NULL
  } else {
    proportions <- found_qc_df
  }

  percentage <- scales::percent(
    sum(flag(object) %in% c("Low_qc_detection", "Low_group_detection"), na.rm = TRUE) / nrow(fData(object))
  )
  log_text(paste0("\n", percentage, " of features flagged for low detection rate"))

  object <- join_fData(object, proportions)

  object
}


#' Flag contaminants
#'
#' Flags contaminant features by comparing the median values of blanks and biological samples.
#' Biological sampels are defined as samples that are not marked as blanks and are not QCs.
#' If the median of blanks > the median of biological samples times a set ratio, the feature is
#' flagged as contaminant
#'
#' @param object a MetaboSet object
#' @param blank_col character, the column name in pData with blank labels
#' @param blank_label character, the label for blank samples in blank_col
#' @param flag_thresh numeric, the ratio threshold for flagging contaminants.
#' If the median of blanks > flag_thresh * median of biological samples, the feature gets flagged.
#' @param flag_label character, the label used when flagging contaminants. Can be changed if
#' sample processing contaminants and carryover contaminants are flagged separately.
#'
#' @export
flag_contaminants <- function(object, blank_col, blank_label, flag_thresh = 0.05,
                              flag_label = "Contaminant") {
  blanks <- object[, pData(object)[, blank_col] == blank_label]
  samples <- object[, object$QC != "QC" & pData(object)[, blank_col] != blank_label]

  blank_median <- apply(exprs(blanks), 1, finite_median)
  sample_median <- apply(exprs(samples), 1, finite_median)
  blank_flag <- blank_median / sample_median > flag_thresh

  idx <- is.na(flag(object)) & !is.na(blank_flag)
  idx <- idx & blank_flag
  flag(object)[idx] <- flag_label

  percentage <- scales::percent(sum(flag(object) == flag_label, na.rm = TRUE) / nrow(object))
  log_text(paste0("\n", percentage, " of features flagged as contaminants"))

  blank_ratio <- data.frame(
    Feature_ID = featureNames(object), Blank_ratio = blank_median / sample_median,
    stringsAsFactors = FALSE
  )
  object <- join_fData(object, blank_ratio)

  object
}
