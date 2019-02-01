

lcms_data <- function(data, response_info, responses, response_name_col,
                      order_col, group_col, sample_id_col = NULL,
                      qc_col = NULL, time_col = NULL, subject_id_col = NULL){

  # All responses should have peak info
  if(!response_name_col %in% colnames(response_info)){
    stop("response_name_col does not match a column name in response_info")
  }
  if (length(responses) != nrow(response_info)) {
    stop("The length of responses and number of rows in peak_info does not match")
  }
  if (!all(responses %in% response_info[, response_name_col])){
    stop("Some responses do not have matching entries in response_info")
  }

  # All responses should be column names of data
  if (!all(responses %in% colnames(data))) {
    stop("Some responses not found in column names of 'data'")
  }

  # All given columns should be in data
  columns <- c(order_col, group_col, sample_id_col, qc_col,
               time_col, subject_id_col)
  columns <- columns[!is.null(columns)]
  missing_columns <- columns[!columns %in% colnames(data)]
  if(length(missing_columns)){
    stop(paste(paste(missing_columns, collapse = ", ")), " columns not found in data")
  }

  # Injection order should be unique
  if (length(unique(data[, order_col])) != nrow(data)) {
    stop("Injection order is not unique")
  }

  if (is.null(qc_col)) {
    if("QC" %in% data[, group_col]){
      qc_col <- group_col
    } else {
      stop("Group column does not have any 'QC' labels, please specify a QC column")
    }
  } else {
    if(!"QC" %in% data[, qc_col]){
      stop("QC column does not have any 'QC' labels")
    }
  }

  # Injection order can be used as sample ID
  if (is.null(sample_id_col)) {
    sample_id_col <- order_col
  }
  # Sample ID should be unique
  if (length(unique(data[, sample_id_col])) != nrow(data)) {
    stop("Sample ID is not unique")
  }

  # Construct the object
  x <- list(data = data, response_info = response_info, responses = responses, response_name = response_name_col,
            order = order_col, group = group_col, sample_id = sample_id_col,
            qc = qc_col, time = time_col, subject_id = subject_id_col)
  class(x) <- c("lcms_data", "list")
  x
}



read_from_excel <- function(file, corner_row, corner_column){
  dada <- openxlsx::read.xlsx(file, colNames = FALSE)

  cc <- corner_column
  cr <- corner_row

  # Extract sample information
  sample_info <- as.data.frame(t(dada[1:cr, (cr+1):ncol(dada)]), stringsAsFactors = FALSE)
  colnames(sample_info) <- gsub(" ", "_", c(dada[1:(cr-1), cc], "Datafile"))

  # Exctract compound information
  response_info <- dada[(cr+1):nrow(dada), 1:cc]
  colnames(response_info) <- dada[cr, 1:cc]

  response_data <- dada[idx + cr, (cr+1):ncol(dada)] %>%
    lapply(as.numeric) %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame()
  colnames(response_data) <- response_info$Compound
  response_data$Datafile <- t(dada[cr, (cr+1):ncol(dada)])[, 1]



}


construct_excel <- function(file, corner_row, corner_column, split_by) {
  # Create Compound column by concatenating rounded mass and retention time
  round_mz <- as.numeric(compound.info$`Average Mz`) %>% round(digits = 4) %>%
    as.character() %>% gsub("[.]", "_", .)
  round_rt <- as.numeric(compound.info$`Average Rt(min)`) %>% round(digits = 4) %>%
    as.character() %>% gsub("[.]", "_", .)
  compound.info$Compound <- paste0("COMPOUND_", round_mz, "a", round_rt)
  compound.info <- compound.info %>%
    select(Compound, Mode, everything())
}
