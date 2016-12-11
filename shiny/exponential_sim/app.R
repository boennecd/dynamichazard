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

# Global params
t_max <- 30
max_rugs <- 1e3
start_fun <- function(t_0 = t_0, t_max = t_max) max(0, runif(1, t_0 - t_max, t_max - 1 - 1e-8))

# params for UI
col_w <- 3

# Thanks to this guy http://stackoverflow.com/a/31066997
n_series_stuff <- list(base = 2, exp_min = 8, exp_max = 12)
JScode <- with(n_series_stuff, paste0(
  "$(function() {
  setTimeout(function(){
  var vals = [", base^exp_min, "];
  var powStart = ", exp_min + 1,";
  var powStop = ", exp_max, ";
  for (i = powStart; i <= powStop; i++) {
  var val = Math.pow(", base,", i);
  val = parseFloat(val.toFixed(8));
  vals.push(val);
  }
  $('#n_series').data('ionRangeSlider').update({'values':vals})
  }, 5)})"))

# Define UI for application that draws a histogram
ui <- fluidPage(

  tags$head(tags$script(HTML(JScode))),

 titlePanel("ddhazard demo"),

 plotOutput("coef_plot"),

 textOutput("n_deaths"),

 br(),

 actionButton("compute", "Compute"),

 hr(),

 fluidRow(
   column(col_w,
          sliderInput("n_series",
                      "Number of series to simulate",
                      min = 0,
                      max = n_series_stuff$exp_max - n_series_stuff$exp_min - 1,
                      value = 1),

          sliderInput("obs_time",
                      "Observed time",
                      min = 1,
                      max = t_max,
                      step = 1,
                      value = ceiling(t_max / 3 * 2)),

          selectInput("sim_with",
                      "Choose model to simulate from",
                      choices = c("logit", "exponential"),
                      selected = "exponential"),

          numericInput("seed",
                       label = "RNG seed",
                       value = 65848)),

   column(col_w,  offset = .5,
          sliderInput("ridge_eps",
                      "Ridge regresion like penalty factor",
                      min = 0.00001,
                      max = .05,
                      step = 0.00001,
                      value = .001),

          selectInput("est_with_model",
                      "Choose model to estimate with",
                      choices = c("logit", "exponential", "exponential_binary_only",
                                  "exponential_trunc_time_only"),
                      selected = "exponential"),

          selectInput("est_with_method",
                      "Choose method to use in the E-step",
                      choices = c("UKF", "EKF"),
                      selected = "EKF")),

   column(col_w, offset = .5,
          h4("EKF settings"),

          checkboxInput("use_extra_correction",
                        "Extra correction steps",
                        value = FALSE)),

   column(col_w, offset = .5,

          h4("UKF settings"),

          sliderInput("beta",
                      "Beta",
                      min = 0,
                      max = 2,
                      step = .5,
                      value = 2),

          sliderInput("alpha",
                      "Alpha",
                      min = 1e-2,
                      max = 1,
                      step = 1e-2,
                      value = 1e-1))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {

  n_series_input <- eventReactive(input$compute, {
    n_series_stuff$base^(n_series_stuff$exp_min + input$n_series)
  })

  sim_input <- eventReactive(input$compute, {
    print(n_series_stuff$base^(n_series_stuff$exp_min + input$n_series))

    set.seed(input$seed)
    f_choice <- if(input$sim_with == "exponential")
      test_sim_func_exp else test_sim_func_logit
    f_choice(
      n_series = n_series_input(),
      n_vars = 5,
      t_max = t_max, re_draw = T, beta_start = runif(5, min = -1.5, max = 1.5),
      intercept_start = -3.5, sds = c(.25, rep(.5, 5)),
      x_range = 1, x_mean = 0, lambda = t_max / 10,
      tstart_sampl_func = start_fun)
  })

  fit_input <- eventReactive(input$compute, {
    sims <- sim_input()
    ddhazard(
      formula = survival::Surv(tstart, tstop, event) ~ . - id - tstart - tstop - event,
      data = sims$res,
      by = 1,
      Q_0 = if(input$est_with_method == "UKF") diag(1, 6) else diag(10, 6),
      Q = diag(1e-1, 6),
      control = list(est_Q_0 = F, eps = 10^-2, n_max = 10^2,
                     save_data = F, save_risk_set = F,
                     ridge_eps = input$ridge_eps,
                     method = input$est_with_method,
                     beta = input$beta, alpha = input$alpha,
                     NR_eps = if(input$use_extra_correction) 1e-1 else NULL),
      max_T = input$obs_time,
      id = sims$res$id, order = 1,
      verbose = F,
      model = input$est_with_model)
  })

  output$n_deaths <- renderText({
    sprintf("%d of %d dies", sum(sim_input()$res$event), n_series_input())
  })

  output$coef_plot <- renderPlot({
    sims <- sim_input()
    fit <- fit_input()

    matplot(seq_len(dim(sims$beta)[1]) - 1, sims$beta, lty = 1, type = "l",
            ylim = range(sims$beta, fit$state_vecs), xaxt='n',
            ylab = expression(beta), xlab = "Time")
    matplot(seq_len(dim(fit$state_vecs)[1]) - 1, fit$state_vecs, lty = 2, type = "l", add = T)

    axis(side = 1,at = seq_len(dim(sims$beta)[1]) - 1,tick = FALSE)

    # Add rug plots to illustrate survivers and deaths
    death_times <- sims$res$tstop[sims$res$event==1]
    if(length(death_times) > max_rugs){
      death_times <- sample(death_times, max_rugs, replace = F)
    }
    rug(sims$res$tstop[sims$res$event==1], line = -.25, col = rgb(0,0,0,.05))

    surv_times <- sims$res$tstop[sims$res$event==0]
    if(length(sims$res$tstop[sims$res$event==0]) > max_rugs){
      surv_times <- sample(surv_times, size = max_rugs, replace = F)
    }
    rug(surv_times, line = .75, col = rgb(0,0,0,.05))

  })
}

# Run the application
shinyApp(ui = ui, server = server)

