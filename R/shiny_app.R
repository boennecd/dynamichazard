#' A Shiny app to illustrates model and method
#'
#' Runs a shiny app where you try different model specification on simulated data
#' @export
#' @examples
#'\dontrun{
#'dynamichazard::ddhazard_app()
#'}
ddhazard_app = function() {
  for(pkg in c("shiny", "formatR")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg, " needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }

  shiny::runApp(system.file('shiny', package = 'dynamichazard'))
}
