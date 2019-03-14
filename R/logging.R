
start_log <- function(log_file) {

  options(amp.logging = TRUE, amp.log_file = log_file)

  cat(paste0("Logging started: ", log_file, ".\n")
  write(paste(date(), "\n", sep=""), log_file)
}

log_text <- function(text) {

  cat(paste(text, "\n", sep=""))
  if (getOption("amp.logging")) {
    if(is.null(getOption("amp.log_file"))) {
      stop("Log file is not defined.")
    }
    write(text, getOption("amp.log_file")), append=TRUE)
  }
}

end_log <- function() {
  # Log end of session info
  log_text(paste("Finished analysis. ", date(), "\nSession info:\n", sep=""))
  log_text(capture.output(sessionInfo()))
  options(amp.logging = FALSE, amp.log_file = NULL)
}
