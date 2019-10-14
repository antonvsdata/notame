#' @export
setGeneric("quality", signature = "object",
           function(object) standardGeneric("quality"))

#' @describeIn MetaboSet extract quality information of features
#' @export
setMethod("quality", c(object = "MetaboSet"),
          function(object) {
            if (!all(c("RSD", "RSD_r", "D_ratio","D_ratio_r") %in% colnames(results(object)))) {
              return(NULL)
            }
            results(object)[c("Feature_ID", "RSD", "RSD_r", "D_ratio",
                            "D_ratio_r")]
          })

#' @export
setGeneric("erase_quality", signature = "object",
           function(object) standardGeneric("erase_quality"))

#' @export
setMethod("erase_quality", c(object = "MetaboSet"),
          function(object) {
            if (!all(c("RSD", "RSD_r", "D_ratio","D_ratio_r") %in% colnames(results(object)))) {
              return(NULL)
            }
            results(object)[c("RSD", "RSD_r", "D_ratio", "D_ratio_r")] <- NULL
            object
          })

#' @export
setGeneric("assess_quality", signature = "object",
           function(object) standardGeneric("assess_quality"))

#' @describeIn MetaboSet compute quality metrics
#' @importFrom foreach "%dopar%"
#' @importFrom Biobase exprs
#' @export
setMethod("assess_quality", c(object = "MetaboSet"),
          function(object) {
            # Remove old quality metrics
            if (!is.null(quality(object))) {
              object <- erase_quality(object)
            }

            qc_data <- exprs(object)[, object$QC == "QC"]
            sample_data <- exprs(object)[, object$QC != "QC"]

            quality_metrics <- foreach::foreach(i = seq_len(nrow(sample_data)), .combine = rbind,
                                                .export = c("finite_sd", "finite_mad", "finite_mean", "finite_median")) %dopar% {
              data.frame(Feature_ID = rownames(sample_data)[i],
                         RSD = finite_sd(qc_data[i, ]) / abs(finite_mean(qc_data[i, ])),
                         RSD_r = finite_mad(qc_data[i, ]) / abs(finite_median(qc_data[i, ])),
                         D_ratio = finite_sd(qc_data[i, ]) / finite_sd(sample_data[i, ]),
                         D_ratio_r = finite_mad(qc_data[i, ]) / finite_mad(sample_data[i, ]),
                         row.names = rownames(sample_data)[i], stringsAsFactors = FALSE)
            }

            object <- join_results(object, quality_metrics)

            object
          })


#' Flag low-quality features
#'
#' Flags low-quality features using the quality metrics defined in (Broadhurst 2018). The metrics are described in
#' more detain in Details. A condition for keeping the features is given as a character, which is passed to \code{dplyr::filter}.
#'
#' @param object a MetaboSet object
#' @param condition character, condition for keeping the features, see Details
#'
#' @details The quality metrics measure two things: internal spread of the QCs, and spread of the QCs compared to the spread of the biological samples.
#'   Internal spread is measured with relative standard deviation (RSD), also known as coefficient of variation (CV).
#'   \deqn{RSD = \frac{s_{QC}}{\bar{x}_{QC}}}
#'   Where\eqn{s_{QC}} is the standard deviation of the QC samples and \eqn{\bar{x}_{QC}} is the sample mean of the signal in the QC samples.
#'   RSD can also be replaced by a non-parametric, robust version based on the median and median absolute deviation (MAD):
#'   \deqn{RSD\_r = \frac{1.4826 \cdot MAD_{QC}}{\text{median}(x_{QC}}}
#'   The spread of the QC samples compared to the biological samples is measured using a metric called D-ratio:
#'   \deqn{D\_ratio = \frac{s_{QC}}{s_{biological}}}
#'   Or, as before, a non-parametric, robust alternative:
#'   \deqn{D\_ratio\_r = \frac{MAD_{QC}}{MAD_{biological}}}
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
#' ex_set <- flag_quality(example_set)
#' results(ex_set)
#' # Custom condition
#' ex_set <- flag_quality(example_set, condition = "RSD_r < 0.3 & D_ratio_r < 0.6")
#' results(ex_set)
#'
#' @export
setGeneric("flag_quality", signature = "object",
           function(object,
                    condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) | (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") standardGeneric("flag_quality"))

#' @describeIn MetaboSet flag low-quality features
#' @export
setMethod("flag_quality", c(object = "MetaboSet"),
          function(object,
                   condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) |
                                (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") {
            if (is.null(quality(object))) {
              object <- assess_quality(object)
            }

            good <- paste0("results(object) %>% dplyr::filter(", condition, ")") %>%
              parse(text = .) %>% eval()
            good <- good$Feature_ID

            idx <- is.na(flag(object)) & !results(object)$Feature_ID %in% good
            flag(object)[idx] <- "Low_quality"

            percentage <- scales::percent(sum(flag(object) == "Low_quality", na.rm = TRUE)/nrow(results(object)))
            log_text(paste0("\n", percentage, " of features flagged for low quality"))

            object
          })

#' Flag features with low detection rate
#'
#' Flags features with too high amount of missing values. There are two detection rate limits, both defined as the minimum proportion of samples that need to
#' have a value (not NA) for the feature to be kept. \code{qc_limit} is the detection rate limit for QC samples, \code{group_limit} is the detection rate limit
#' for the actual study groups. If the group limit is passed for AT LEAST ONE GROUP, then the feature is kept. Features with low detection rate in QCs are
#' flagged as "Low_qc_detection", while low detection rate in the study groups is flagged as "Low_group_detection". The detection rates for all the groups are
#' recorded in \code{results(object)}
#'
#' @param object a MetaboSet object
#' @param qc_limit the detection rate limit for QC samples
#' @param group_limit the detection rate limit for study groups
#' @param group the columns name in sample information to use as the grouping variable
#'
#' @return a MetaboSet object with the features flagged
#'
#' @examples
#' ex_set <- flag_detection(example_set)
#' results(ex_set)
#'
#' @export
setGeneric("flag_detection", signature = "object",
           function(object, qc_limit = 0.7, group_limit = 0.8, group = group_col(object)) standardGeneric("flag_detection"))

#' @describeIn MetaboSet flag features with low detection rate
#' @importFrom Biobase featureNames
#' @export
setMethod("flag_detection", c(object = "MetaboSet"),
          function(object, qc_limit, group_limit, group) {

            found_qc <- Biobase::esApply(object[, object$QC == "QC"], 1, prop_found)
            bad_qc <- names(which(found_qc < qc_limit))

            found_qc_df <- data.frame(Feature_ID = names(found_qc),
                                      Detection_rate_QC = found_qc,
                                      stringsAsFactors = FALSE)

            idx <- is.na(flag(object)) & results(object)$Feature_ID %in% bad_qc
            flag(object)[idx] <- "Low_qc_detection"

            # Compute proportions found in each study group
            if (!is.na(group)) {
              proportions <- combined_data(object)[, c("Sample_ID", group, featureNames(object))] %>%
                tidyr::gather_("Feature_ID", "Intensity", featureNames(object)) %>%
                dplyr::group_by_("Feature_ID", group) %>%
                dplyr::summarise(proportion_found = prop_found(Intensity)) %>%
                tidyr::spread_(group, "proportion_found")
              # Remove a possible QC column
              proportions$QC <- NULL
              colnames(proportions)[-1] <- paste0("Detection_rate_", group, "_", colnames(proportions)[-1])
              # Check if any group has enough non-missing entries
              proportions$good <- apply(proportions[-1], 1, function(x){any(x >= group_limit)})

              idx <- is.na(flag(object)) & (!results(object)$Feature_ID %in% proportions$Feature_ID[proportions$good])
              flag(object)[idx] <- "Low_group_detection"
              # Add detection rates to feature data
              proportions <- dplyr::left_join(proportions, found_qc_df, by = "Feature_ID")
              proportions$good <- NULL
            } else {
              proportions <- found_qc_df
            }

            percentage <- scales::percent(sum(flag(object) %in% c("Low_qc_detection", "Low_group_detection"), na.rm = TRUE)/nrow(results(object)))
            log_text(paste0("\n", percentage, " of features flagged for low detection rate"))

            object <- join_results(object, proportions)

            object
          })
