#' @title ddhazard Demo
#'
#' @description
#' \code{ddhazard_app} runs a shiny app with demonstration of models.
#'
#' @param quietly \code{TRUE} if no messages should be printed when the app is run.
#' @param ... starting values for the shiny app.
#'
#' @details
#' Runs a shiny app where you try different model specifications on simulated data.
#'
#' @examples
#'\dontrun{
#' dynamichazard::ddhazard_app()
#' dynamichazard::ddhazard_app(seed = 1, more_options = TRUE)
#'}
#' @export
ddhazard_app = function(quietly = F, ...) {
  for(pkg in c("shiny", "formatR")){
    if (!requireNamespace(pkg, quietly = TRUE)) {
      stop(pkg, " needed for this function to work. Please install it.",
           call. = FALSE)
    }
  }

  input_args <- list(...)

  ui <- NULL
  server <- NULL
  source(system.file('shiny/app.R', package = 'dynamichazard'), local = TRUE)

  shiny::shinyApp(ui = ui, server = server)
}
