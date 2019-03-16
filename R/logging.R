
#' Initialize log to a file
#'
#' Initialize a log file with the current data and time.
#' All major operations run after this will be logged to the specified file.
#'
#' @section Warning:
#' This overwrites the current contents of the file
#'
#' @param log_file Path to the log file
#'
#' @examples
#' file_name <- "~/log.txt"
#' init_log(file_name)
#' # Print the contents of the file
#' scan(file_name, sep="\n", what = "chracter")
#'
#' @seealso \code{\link{log_text}}, \code{\link{finish_log}}, \code{\link{log_state}}
#'
#' @export
init_log <- function(log_file) {

  options(amp.logging = TRUE, amp.log_file = log_file)

  cat(paste0("Logging started: ", log_file, ".\n"))
  write(paste(date(), "\n", sep=""), log_file)
}

#' Log text to the current log file
#'
#' The specified text is printed and written to the current log file. Does not overwrite the file.
#' Also used internally by many functions in the package
#'
#' @param text The text to be logged
#'
#' @examples
#' file_name <- "~/log.txt"
#' init_log(file_name)
#' log_text("Hello World!")
#' # Print the contents of the file
#' scan(file_name, sep="\n", what = "chracter")
#'
#' @seealso \code{\link{init_log}}, \code{\link{finish_log}}, \code{\link{log_state}}
#'
#' @export
log_text <- function(text) {

  cat(paste(text, "\n", sep=""))
  if (getOption("amp.logging")) {
    if(is.null(getOption("amp.log_file"))) {
      stop("Log file is not defined.")
    }
    write(text, getOption("amp.log_file"), append=TRUE)
  }
}

#' Finish a log
#'
#' Logs the current date and time and session info, and switches logging off.
#'
#' @seealso \code{\link{init_log}}, \code{\link{log_text}}, \code{\link{log_state}}
#'
#' @export
finish_log <- function() {
  # Log end of session info
  log_text(paste("Finished analysis. ", date(), "\nSession info:\n", sep=""))
  log_text(capture.output(sessionInfo()))
  options(amp.logging = FALSE, amp.log_file = NULL)
}


#' Get the current logging state
#'
#' If a log file is currently in use, prints the log file,
#' else tells you that logging is not enabled
#'
#' @seealso \code{\link{init_log}}, \code{\link{log_text}}, \code{\link{finish_log}}
#'
#' @export
log_state <- function() {
  if (!getOption("amp.logging")) {
    cat("Logging is not enabled")
  } else {
    cat(paste("Current log file:", getOption("amp.log_file")))
  }
}
