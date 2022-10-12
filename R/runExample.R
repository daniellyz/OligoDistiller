#' @export

runExample <- function() {
  appDir <- system.file("shiny-examples", "myapp", package = "OligoDistiller")
  if (appDir == "") {
    stop("Could not find example directory. Try re-installing `OligoDistiller`.", call. = FALSE)
  }
  
  shiny::runApp(appDir, display.mode = "normal")
}