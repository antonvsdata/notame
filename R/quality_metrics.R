#' @export
setGeneric("quality", signature = "object",
           function(object) standardGeneric("quality"))

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


#' @export
setGeneric("flag_quality", signature = "object",
           function(object,
                    condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) |
                                (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") standardGeneric("flag_quality"))

#' @export
setMethod("flag_quality", c(object = "MetaboSet"),
          function(object,
                   condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) |
                                (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") {

            good <- paste0("results(object) %>% dplyr::filter(", condition, ")") %>%
              parse(text = .) %>% eval()
            good <- good$Feature_ID

            idx <- is.na(results(object)$Flag) & !results(object)$Feature_ID %in% good
            results(object)$Flag[idx] <- "Low_quality"

            object
          })


#' @export
setGeneric("flag_detection", signature = "object",
           function(object, qc_limit = 0.7, group_limit = 0.8, group = group_col(object)) standardGeneric("flag_detection"))

#' @importFrom Biobase featureNames
#' @export
setMethod("flag_detection", c(object = "MetaboSet"),
          function(object, qc_limit, group_limit, group) {

            found_qc <- Biobase::esApply(object[, object$QC == "QC"], 1, prop_found)
            bad_qc <- names(which(found_qc < qc_limit))

            found_qc_df <- data.frame(Feature_ID = names(found_qc),
                                      Detection_rate_QC = found_qc,
                                      stringsAsFactors = FALSE)

            idx <- is.na(results(object)$Flag) & results(object)$Feature_ID %in% bad_qc
            results(object)$Flag[idx] <- "Low_qc_detection"

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

              idx <- is.na(results(object)$Flag) & (!results(object)$Feature_ID %in% proportions$Feature_ID[proportions$good])
              results(object)$Flag[idx] <- "Low_group_detection"
              # Add detection rates to feature data
              proportions <- dplyr::left_join(proportions, found_qc_df, by = "Feature_ID")
              proportions$good <- NULL
            } else {
              proportions <- found_qc_df
            }



            object <- join_results(object, proportions)

            object
          })
