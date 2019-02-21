

setGeneric("assess_quality", signature = "object",
           function(object) standardGeneric("assess_quality"))

#' @importFrom foreach "%dopar%"
#' @importFrom Biobase exprs fData "fData<-"
setMethod("assess_quality", c(object = "MetaboSet"),
          function(object) {
            qc_data <- exprs(object)[, object$QC == "QC"]
            sample_data <- exprs(object)[, object$QC != "QC"]

            quality_metrics <- foreach::foreach(i = 1:nrow(sample_data), .combine = rbind) %dopar% {
              data.frame(Feature_ID = rownames(sample_data)[i],
                         RSD = finite_sd(qc_data[i, ]) / abs(finite_mean(qc_data[i, ])),
                         RSD_r = finite_mad(qc_data[i, ]) / abs(finite_median(qc_data[i, ])),
                         D_ratio = finite_sd(qc_data[i, ]) / finite_sd(sample_data[i, ]),
                         D_ratio_r = finite_mad(qc_data[i, ]) / finite_mad(sample_data[i, ]),
                         stringsAsFactors = FALSE)
            }

            object <- join_fdata(object, quality_metrics)

            object
          })


setGeneric("quality", signature = "object",
           function(object) standardGeneric("quality"))

#' @importFrom Biobase fData
setMethod("quality", c(object = "MetaboSet"),
          function(object) {
            fData(object)[c("Feature_ID", "RSD", "RSD_r", "D_ratio",
                                     "D_ratio_r")]
          })



setGeneric("flag_quality", signature = "object",
           function(object,
                    condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) |
                                (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") standardGeneric("flag_quality"))

#' @importFrom Biobase fData "fData<-"
setMethod("flag_quality", c(object = "MetaboSet"),
          function(object,
                   condition = "(RSD_r < 0.2 & D_ratio_r < 0.4) |
                                (RSD < 0.1 & RSD_r < 0.1 & D_ratio < 0.1)") {

            good <- paste0("fData(object) %>% dplyr::filter(", condition, ")") %>%
              parse(text = .) %>% eval()
            good <- good$Feature_ID

            idx <- is.na(fData(object)$Flag) & !fData(object)$Feature_ID %in% good
            fData(object)$Flag[idx] <- "Low_quality"

            object
          })


setGeneric("flag_detection", signature = "object",
           function(object, qc_limit = 0.7, group_limit = 0.8, group = group_col(object)) standardGeneric("flag_detection"))

#' @importFrom Biobase fData "fData<-" featureNames
#' @importFrom magrittr "%>%"
setMethod("flag_detection", c(object = "MetaboSet"),
          function(object, qc_limit, group_limit, group) {

            found_qc <- Biobase::esApply(object[, object$QC == "QC"], 1, prop_found)
            bad_qc <- names(which(found_qc < qc_limit))

            found_qc_df <- data.frame(Feature_ID = names(found_qc),
                                      Detection_rate_QC = found_qc,
                                      stringsAsFactors = FALSE)

            idx <- is.na(fData(object)$Flag) & fData(object)$Feature_ID %in% bad_qc
            fData(object)$Flag[idx] <- "Low_qc_detection"

            # Compute proportions found in each study group
            if (!is.na(group)) {
              proportions <- lcms_data(object)[, c("Sample_ID", group, featureNames(object))] %>%
                tidyr::gather_("Feature_ID", "Intensity", featureNames(object)) %>%
                dplyr::group_by_("Feature_ID", group) %>%
                dplyr::summarise(proportion_found = prop_found(Intensity)) %>%
                tidyr::spread_(group, "proportion_found")
              # Remove a possible QC column
              proportions$QC <- NULL
              colnames(proportions)[-1] <- paste0("Detection_rate_", group, "_", colnames(proportions)[-1])
              # Check if any group has enough non-missing entries
              proportions$good <- apply(proportions[-1], 1, function(x){any(x >= group_limit)})

              idx <- is.na(fData(object)$Flag) & (!fData(object)$Feature_ID %in% proportions$Feature_ID[proportions$good])
              fData(object)$Flag[idx] <- "Low_group_detection"
              # Add detection rates to feature data
              proportions <- dplyr::left_join(proportions, found_qc_df, by = "Feature_ID")
              proportions$good <- NULL
            } else {
              proportions <- found_qc_df
            }



            object <- join_fdata(object, proportions)

            object
          })

