

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
                         Detection_rate = 1 - prop_na(qc_data[i, ]),
                         stringsAsFactors = FALSE)
            }

            fData(object) <- dplyr::left_join(fData(object),
                                                       quality_metrics,
                                                       by = "Feature_ID")
            rownames(fData(object)) <- fData(object)$Feature_ID

            object
          })


setGeneric("quality", signature = "object",
           function(object) standardGeneric("quality"))

#' @importFrom Biobase fData
setMethod("quality", c(object = "MetaboSet"),
          function(object) {
            fData(object)[c("Feature_ID", "RSD", "RSD_r", "D_ratio",
                                     "D_ratio_r", "Detection_rate")]
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


setGeneric("flag_missigness", signature = "object",
           function(object, qc_limit, group_limit, group) standardGeneric("flag_missigness"))

#' @importFrom Biobase fData "fData<-"
#' @importFrom magrittr "%>%"
setMethod("flag_missigness", c(object = "MetaboSet"),
          function(object, qc_limit = 0.7, group_limit = 0.8, group = group(object)) {



            object
            })

