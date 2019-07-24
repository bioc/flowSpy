#' running the shiny web app.
#'
#' @name runFlowSpyApp
#' @description running the shiny web app.
#' @noRd
#'
runFlowSpyApp <- function(...){
  if(requireNamespace("shiny", quietly=TRUE)){
    message("Starting the flowSpy shiny web app. Use Ctrl-C to stop.")
    shiny::runApp(appDir=system.file("flowSpyShiny",
                                     package="flowSpyShiny"), ...)
  }else{
    stop("Package shiny not installed!")
  }
}



