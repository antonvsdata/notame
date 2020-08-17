

#' Set default color scales on load
#'
#' @param libname,pckgname default parameters
#'
#' @importFrom grDevices rgb
.onLoad <- function(libname, pkgname) {
  op <- options()
  op.notame <- list(
    notame.logging = FALSE,
    notame.log_file = NULL,
    notame.citations = list("Preprocessing and analyses were performed using notame package:" = citation("notame"),
                            "notame is built on a class from Biobase package:" = citation("Biobase"),
                            "visualizations in notame are built with ggplot2:" = citation("ggplot2")),
    notame.color_scale_con = ggplot2::scale_color_viridis_c(),
    notame.color_scale_dis = ggplot2::scale_color_brewer(palette = "Set1"),
    notame.fill_scale_con = ggplot2::scale_fill_viridis_c(),
    notame.fill_scale_dis = ggplot2::scale_fill_brewer(palette = "Set1"),
    notame.fill_scale_div_con = ggplot2::scale_fill_distiller(palette = "RdBu"),
    notame.fill_scale_div_dis = ggplot2::scale_fill_brewer(palette = "RdBu"),
    notame.shape_scale = ggplot2::scale_shape_manual(values = c(16, 17, 15, 3, 7, 8, 11, 13))
  )
  toset <- !(names(op.notame) %in% names(op))
  if(any(toset)) options(op.notame[toset])

  invisible()
}

install_helper <- function(cran, bioconductor, github, gitlab, ...) {
  if (!missing(cran)) {
    for (pckg in cran) {
      if (!requireNamespace(pckg, quietly = TRUE)) {
        cat(paste("\nPackage", pckg, "missing, attempting to install from CRAN\n"))
        tryCatch({
          install.packages(pckg, ...)
        }, error = function(e) {cat(e$message)})
      }
    }
  }

  if (!missing(bioconductor)) {
    for (pckg in bioconductor) {
      if (!requireNamespace(pckg, quietly = TRUE)) {
        cat(paste("\nPackage", pckg, "missing, attempting to install from Bioconductor\n"))
        tryCatch({
          BiocManager::install(pckg, ...)
        }, error = function(e) {cat(e$message)})
      }
    }
  }

  if (!missing(github)) {
    for (pckg in github) {
      if (!requireNamespace(strsplit(pckg, split = "/")[[1]][2], quietly = TRUE)) {
        cat(paste("\nPackage", pckg, "missing, attempting to install from GitHub\n"))
        tryCatch({
          devtools::install_github(pckg, ...)
        }, error = function(e) {cat(e$message)})
      }
    }
  }

  if (!missing(gitlab)) {
    for (pckg in gitlab) {
      if (!requireNamespace(strsplit(pckg, split = "/")[[1]][2], quietly = TRUE)) {
        cat(paste("\nPackage", pckg, "missing, attempting to install from GitLab\n"))
        tryCatch({
          devtools::install_gitlab(pckg, ...)
        }, error = function(e) {cat(e$message)})
      }
    }
  }
}

#' Install dependencies
#'
#' Attempt to install dependencies package by package, skipping packages with errors
#' By default, only installs core packages needed for preprocessing. Other packages can
#' be installed if needed (and the preprocessing packages can be ignored)
#'
#' @param preprocessing logical, install core preprocessing and visualization packages?
#' @param extra logical, install extra packages needed for special visualizations and stats?
#' @param batch_corr logical, install packages needed for batch effect correction methods?
#' @param misc logicl, install miscallenous packages needed for running tests, modifying vignettes etc.?
#' @param ... other parameters passed to installing functions, like lib
#'
#' @export
install_dependencies <- function(preprocessing = TRUE, extra = FALSE, batch_corr = FALSE, misc = FALSE, ...) {
  # Core dependencies
  core_cran <- c("BiocManager",
                 "cowplot",
                 "missForest",
                 "openxlsx",
                 "randomForest",
                 "RColorBrewer",
                 "Rtsne")
  core_bioconductor <- "pcaMethods"
  # Extra parts for certain visualizations and statistics
  extra_cran <- c("car",
                  "doParallel",
                  "ggbeeswarm",
                  "ggdendro",
                  "ggrepel",
                  "Hmisc",
                  "hexbin",
                  "igraph",
                  "lme4",
                  "lmerTest",
                  "MuMIn",
                  "PK",
                  "rmcorr")
  extra_bioconductor <- c("mixOmics", "supraHex")
  extra_gitlab <- "CarlBrunius/MUVR"

  batch_cran <- "fpc"
  batch_bioconductor <- "RUVSeq"
  batch_github <- "rwehrens/BatchCorrMetabolomics"
  batch_gitlab <- "CarlBrunius/batchCorr"

  misc_cran <- c("knitr",
                 "rmarkdown",
                 "testthat")

  if (preprocessing){
    install_helper(cran = core_cran, bioconductor = core_bioconductor, ...)
  }

  if (extra) {
    install_helper(cran = extra_cran, bioconductor = extra_bioconductor,
                   gitlab = extra_gitlab, ...)
  }
  if (batch_corr) {
    install_helper(cran = batch_cran, bioconductor = batch_bioconductor,
                   github = batch_github, gitlab = batch_gitlab, ...)
  }
  if (misc) {
    install_helper(cran = misc_cran, ...)
  }

}

add_citation <- function(name, ref) {
  cites <- getOption("notame.citations")
  if (!name %in% names(cites)) {
    cites[[name]] <- ref
    options(notame.citations = cites)
  }
}

#' Show citations
#'
#' This function lists citations for all the major packages used by the notame functions that
#' have been called during the session. All notame functions update the list automatically.
#' The citations are taken from the call to \code{citation("package")}, and complemented with
#' a brief description of what the package was used for.
#' NOTE: the citations might not point to the correct paper if the package authors have not
#' supplied correct citation information for their package.
#' The output is written to the current log file, if specified.
#'
#' @examples
#'
#' citations()
#'
#' plot_tsne(merged_sample)
#'
#' # Rtsne added to citations
#' citations()
citations <- function() {
  cites <- getOption("notame.citations")
  for (i in seq_along(cites)) {
    log_text(names(cites)[i])
    log_text(capture.output(show(cites[[i]])))
  }
}


#' Summary statistics of finite elements
#'
#' These functions first remove non-finite and missing values, then compute the summary statistic in question.
#' They are helper functions used for computing quality measurements.
#'
#' @param x a numeric vector.
#' @name finite_helpers
NULL

#' @export
#' @rdname finite_helpers
finite_sd <- function(x) {
  sd(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_mean <- function(x) {
  if (all(is.na(x))) {
    return(NA_real_)
  }
  mean(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_median <- function(x) {
  median(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_min <- function(x) {
  min(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_max <- function(x) {
  max(x[is.finite(x)], na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_mad <- function(x) {
  mad(x[is.finite(x)], center = median(x[is.finite(x)], na.rm = TRUE), na.rm = TRUE)
}

#' @export
#' @rdname finite_helpers
finite_quantile <- function(x, ...) {
  unname(quantile(x[is.finite(x)], na.rm = TRUE, ...))
}



# Defaults for NULL values
`%||%` <- function(a, b) {
  suppressWarnings(if (is.null(a)){
    b
  } else if (is.na(a)){
    b
  } else {
    a
  })
}

#' Proportion of NA values in a vector
#'
#' @param x a numeric vector
#'
#' @export
prop_na <- function(x) {
  sum(is.na(x)) / length(x)
}

#' Proportion of non-missing values in a vector
#'
#' @param x a numeric vector
#'
#' @export
prop_found <- function(x) {
  sum(!is.na(x)) / length(x)
}

best_class <- function(x) {
  x <- type.convert(as.character(x), as.is = TRUE)
  if (class(x) == "numeric") {
    x <- x
  } else if (length(unique(x)) < length(x)/4) {
    x <- as.factor(x)
  } else if (is.integer(x)) {
    x <- as.numeric(x)
  } else {
    x <- as.character(x)
  }
  x
}


best_classes <- function(x) {
  as.data.frame(lapply(x, best_class), stringsAsFactors = FALSE)
}

all_unique <- function(x) {
  !any(duplicated(x))
}
