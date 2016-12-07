#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(dynamichazard)

test_sim_func_exp <- with(environment(ddhazard), test_sim_func_exp)
test_sim_func_logit <- with(environment(ddhazard), test_sim_func_logit)

# Define UI for application that draws a histogram
ui <- fluidPage(

   # Application title
   titlePanel("Old Faithful Geyser Data"),

   sidebarLayout(
      sidebarPanel(
         sliderInput("n_series",
                     h3("Number of series to simulate"),
                     min = 100,
                     max = 1e4,
                     step = 100,
                     value = 100),

         sliderInput("ridge_eps",
                     h3("Ridge regresion like penalty factor"),
                     min = 0.001,
                     max = .05,
                     step = .001,
                     value = .01),

         selectInput("est_with",
                     h3("Choose model to estimate with"),
                     choices = c("logit", "exponential"),
                     selected = "exponential"),

         selectInput("sim_with",
                     h3("Choose model to simulate from"),
                     choices = c("logit", "exponential"),
                     selected = "exponential"),

         numericInput("seed",
                      label = h3("RNG seed"),
                      value = 625409)
      ),

      mainPanel({
        textOutput("n_deaths")
        plotOutput("coef_plot")
      })
   )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  sim_input <- reactive({
    set.seed(input$seed)
    f_choice <- if(input$sim_with == "exponential")
      test_sim_func_logit else test_sim_func_logit
    f_choice(
      n_series = input$n_series, n_vars = 5,
      t_max = 10, re_draw = T, beta_start = runif(5),
      intercept_start = -4, sds = c(.25, rep(1, 5)),
      x_range = 1, x_mean = 0)
  })

  fit_input <- reactive({
    sims <- sim_input()
    ddhazard(
      formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
      data = sims$res,
      by = 1,
      Q_0 = diag(10, 6),
      Q = diag(1e-2, 6),
      control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^2,
                     save_data = F, save_risk_set = F,
                     ridge_eps = input$ridge_eps),
      max_T = 10,
      id = sims$res$id, order = 1,
      verbose = F,
      model = input$est_with)
  })

  output$n_deaths <- renderText({
    sprintf("%d of %d dies", sum(sim_input()$res$event), input$n_series)
  })

  output$coef_plot <- renderPlot({
    sims <- sim_input()
    fit <- fit_input()

    matplot(sims$beta[-1, ], lty = 1, type = "l", ylim = range(sims$beta, fit$state_vecs))
    matplot(fit$state_vecs[-1, ], lty = 2, type = "l", add = T)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

