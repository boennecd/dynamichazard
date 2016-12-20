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

# Save memory on server
ddhazard <- function(..., control){
  control$save_data = F
  control$save_risk_set = F
  dynamichazard::ddhazard(..., control = control)
}

# Global params
t_max <- 30
max_rugs <- 5e2
start_fun <- function(t_0 = t_0, t_max = t_max) max(0, runif(1, t_0 - t_max, t_max - 1 - 1e-8))

# params for UI
col_w <- 4
col_w_out <- 4

get_em <- function(n_lines) paste0("height:",n_lines * 1.5,"em;")

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

 div(fixedRow(
   column(3,
          h3("Intro"),
          div("Illustrates simulated data and a fit. The true coeffecient are the continous curves and the predicted coeffecients are the dashed curves.", style = get_em(3)),
          div(textOutput("rug_explanation"), style = get_em(5)),
          div("See the ddhazard vignette for further details", style = get_em(1))
   ),
   column(3,
          h3("Output"),
          div(textOutput("n_deaths"), style = get_em(1)),
          div(htmlOutput ("LR_out"), style = get_em(1)),
          div(textOutput("MSE"), style = get_em(1)),
          div(textOutput("model_text"), style = get_em(3)),
          tags$textarea(
            id = 'demo', style = 'display: none;',
            verbatimTextOutput("fit_call_txt")
          )
   ),
   column(5,
          h3("Fit call"),
          div(verbatimTextOutput("fit_call_txt"), style = "height:20em;")
   )
 ), style = "max-width:110em;min-width: 100em;"),

 br(),
 hr(),
 br(),

 div(plotOutput("coef_plot"), style = "max-width:1000px;min-width=200px"),

 br(),
 hr(),
 br(),

 div(style = "width: 60em;",
   fluidRow(
   column(col_w,
          h4("Simulation settings"),

          sliderInput("n_series",
                      "Number of series to simulate",
                      min = 0,
                      max = n_series_stuff$exp_max - n_series_stuff$exp_min - 1,
                      value = 1),

          selectInput("sim_with",
                      "Choose model to simulate from",
                      choices = c("logit", "exponential"),
                      selected = "exponential"),

          radioButtons("sim_fix_options",
                       label = "Number of fixed covariates",
                       choices = list("Zero" = 1,
                                      "Intercept and two coef" = 2,
                                      "All but one coef" = 3),
                       selected = 1),

          conditionalPanel(
            "input.more_options",
            sliderInput("obs_time",
                        "Observed time",
                        min = 1,
                        max = t_max,
                        step = 1,
                        value = t_max),

            numericInput("seed",
                         label = "RNG seed",
                         value = 65848))),

   column(col_w,
          h4("Estimation settings"),

          selectInput("est_with_model",
                      "Choose model to estimate with",
                      choices = c("logit", "exponential", "exponential_binary_only",
                                  "exponential_trunc_time_only"),
                      selected = "exponential"),

          selectInput("est_with_method",
                      "Choose method to use in the E-step",
                      choices = c("UKF", "EKF"),
                      selected = "EKF"),

          radioButtons("est_fix_options",
                       label = "Number of fixed covariates",
                       choices = list("Zero" = 1,
                                      "Intercept and two coef" = 2,
                                      "All but one coef" = 3),
                       selected = 1),

          conditionalPanel(
            "input.more_options",
            sliderInput("order",
                        "Randowm walk order in estimation",
                        min = 1,
                        max = 2,
                        step = 1,
                        value = 1),
            sliderInput("ridge_eps",
                        "Ridge regresion like penalty factor",
                        min = 0.00001,
                        max = .05,
                        step = 0.00001,
                        value = .001),
            selectInput("fixed_terms_method",
                        "Estimate fixed effect in",
                        choices = c("E_step", "M_step"),
                        selected = "M_step"))),

   conditionalPanel(
     "input.est_with_method == 'EKF'",
     column(col_w,
            h4("EKF settings"),

            checkboxInput("use_extra_correction",
                          "Extra correction steps",
                          value = FALSE))),

   conditionalPanel(
     "input.est_with_method == 'UKF'",
     column(col_w,
            h4("UKF settings"),

            sliderInput("beta",
                        "Beta",
                        min = 0,
                        max = 2,
                        step = .5,
                        value = 0),

            sliderInput("alpha",
                        "Alpha",
                        min = 1e-2,
                        max = 1,
                        step = 1e-2,
                        value = 1)))
  ),


  div(style="float:right;",
      flowLayout(
        checkboxInput("more_options", label = "Show more options", value = FALSE, width = "12em"),
        actionButton("compute", "Compute", width = "8em")))
  )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  data <- NA

  n_series_input <- eventReactive(input$compute, {
    n_series_stuff$base^(n_series_stuff$exp_min + input$n_series)
  })

  sim_input <- eventReactive(input$compute, {
    set.seed(input$seed)
    f_choice <- if(input$sim_with == "exponential")
      test_sim_func_exp else test_sim_func_logit

    n_fixed <- if(input$sim_fix_options == 1)
      0 else if(input$sim_fix_options == 2)
        3 else if(input$sim_fix_options == 3)
          5

    n_varying <- 5 - n_fixed + (n_fixed > 0)

    f_choice(
      n_series = n_series_input(),
      n_vars = 5,
      t_max = t_max, re_draw = T, beta_start = runif(5, min = -1.5, max = 1.5),
      intercept_start = -3.5,
      sds = c(.25, rep(.5, 5)),
      x_range = 1, x_mean = 0, lambda = t_max / 10,
      is_fixed = if(n_fixed == 0) c() else 1:n_fixed,

      tstart_sampl_func = start_fun)
  })

  n_fixed_when_est <- eventReactive(input$compute, {
    est_fix_options <- input$est_fix_options

    if(est_fix_options == 1){
      return(0)
    } else if(est_fix_options == 2){
      return(3)
    } else if(est_fix_options == 3){
      return(5)
    } else
      stop("est_fix_options option not implemented")
  })

  fit_quote_input <- eventReactive(input$compute, {
    sims <- sim_input()
    data <<- sims$res

    n_fixed <- n_fixed_when_est()

    if(n_fixed == 0){
      form <- formula(survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3 + x4 + x5)

    } else if(n_fixed == 3){
      form <- formula(survival::Surv(tstart, tstop, event) ~ ddFixed(1) +
        ddFixed(x1) + ddFixed(x2) + x3 + x4 + x5)

    } else if(n_fixed == 5){
      form <- formula(survival::Surv(tstart, tstop, event) ~ ddFixed(1) +
        ddFixed(x1) + ddFixed(x2) + ddFixed(x3) + ddFixed(x4) + x5)

    } else
      stop("n_fixed is not implemented")

    q_0_len <- 6 - n_fixed
    q_0_term <- ifelse(input$est_with_method == "UKF", 1, 10)

    if(input$order == 2){
      if(q_0_len > 1){
        Q_0 <- bquote(diag(c(rep(.(q_0_term), .(q_0_len)),
                             rep(.1, .(q_0_len)))))
      } else
        Q_0 <- bquote(diag(c(.(q_0_term), .1)))

    } else{
      if(q_0_len > 1){
        Q_0 <- bquote(diag(.(q_0_term), .(q_0_len)))
      } else
        Q_0 <- q_0_term

    }

    control_list <- list(eps = 10^-2,
                         ridge_eps = input$ridge_eps,
                         method = input$est_with_method,
                         fixed_terms_method  = input$fixed_terms_method)

    if(input$est_with_method == "UKF"){
      control_list <- c(control_list,
                        list(beta = input$beta, alpha = input$alpha))
    } else if (input$est_with_method == "EKF"){
      if(input$use_extra_correction)
        control_list <- c(control_list,
                          list(NR_eps = .1))
    }

    Q <- if(6 - n_fixed > 1)
      bquote(diag(.1, .(6 - n_fixed))) else 0.1

    list(quote = bquote(ddhazard(
        formula = .(form),
        data = data,
        by = 1,
        Q_0 = .(Q_0),
        Q = .(Q),
        max_T = .(input$obs_time),
        id = data$id,
        order = .(input$order),
        model = .(input$est_with_model),
        control = .(control_list))),
      data = data)
  })

  output$fit_call_txt <- renderText({
    eval_quote <- fit_quote_input()$quote
    out <- formatR::tidy_source(text = capture.output(eval_quote), width.cutoff = 40)
    out <- out$text.tidy
    out <- paste0(out, "\n\n# data is the simulated data set")
  })

  fit_input <- eventReactive(input$compute, {
    tmp <- fit_quote_input()
    eval_quote <- tmp$quote
    data <- tmp$data
    eval(eval_quote)
  })

  output$n_deaths <- renderText({
    sprintf("%d of %d dies", sum(sim_input()$res$event), n_series_input())
  })

  output$LR_out <- renderText({
    LR <- fit_input()$LR
    out <- sprintf(paste0("The used learning rate is ", ifelse(LR == 1, "%d", "%.3f")), LR)

    if(LR < 1)
      out <- paste0("<font color=\"#FF0000\"><b>", out, "</b></font>")

    out
  })

  output$MSE <- renderText({
    n_fixed <- n_fixed_when_est()
    fit <- fit_input()

    state_vecs <- if(n_fixed == 0)
      fit$state_vecs[1:input$obs_time + 1, 1:6] else
        cbind(sapply(fit$fixed_effects, rep, times = input$obs_time),
              fit$state_vecs[1:input$obs_time + 1, 1:(6 - n_fixed)])

    sprintf("MSE for coeffecients is %.3f",
      mean.default((sim_input()$beta[1:input$obs_time + 1, ] - state_vecs)^2))
  })

  output$rug_explanation <- renderText({
    text_out <- if(input$sim_with == "exponential")
      "The lines above x-axis indicates when a death is observed. " else
        "The lines above x-axis indicates when a death is observed. A small jitter is added to distinguish them (this is only when we simulate from the logit model). "

    paste0(text_out, "The lines below the x-axis is when a stop time is observed that is not a death.")
  })

  output$coef_plot <- renderPlot({
    sims <- sim_input()
    fit <- fit_input()
    n_fixed <- n_fixed_when_est()

    matplot(seq_len(dim(sims$beta)[1]) - 1, sims$beta, lty = 1, type = "l",
            ylim = range(sims$beta, fit$state_vecs, fit$fixed_effects), xaxt='n',
            ylab = expression(beta), xlab = "Time")
    matplot(seq_len(dim(fit$state_vecs)[1]) - 1, fit$state_vecs[, 1:(6 - n_fixed)],
            lty = 2, type = "l", add = T,
            col = (1+n_fixed):6)

    if(n_fixed > 0)
      abline(h = fit$fixed_effects, col = 1:n_fixed, lty = 2)

    axis(side = 1, at = seq_len(dim(sims$beta)[1]) - 1,tick = FALSE)

    # Add rug plots to illustrate survivers and deaths
    death_times <- sims$res$tstop[sims$res$event==1]
    if(length(death_times) > max_rugs){
      death_times <- sample(death_times, max_rugs, replace = F)
    }

    if(input$sim_with == "logit")
      death_times <- jitter(death_times, amount = .2)

    rug(death_times, line = -.25, col = rgb(0,0,0,.1), lwd = 1)

    surv_times <- sims$res$tstop[sims$res$event==0]
    if(length(sims$res$tstop[sims$res$event==0]) > max_rugs){
      surv_times <- sample(surv_times, size = max_rugs, replace = F)
    }
    rug(surv_times, line = .75, col = rgb(0,0,0,.1), lwd = 1)

  })

  output$model_text <- renderText({
    result <- fit_input()

    out <- paste0("Model is estimated with order ", result$order, " and with the ", result$method, " in the E-step.",
                  " ", result$n_iter, " iteration was used in the EM algorithm.")

    if(result$method == "EKF"){
      out <- paste0(out, ifelse(is.null(result$control$NR_eps), " No extra", " Extra"),
                    " iterations are used in correction step")
    } else{
      out <- paste0(out, " Alpha and beta are ", result$control$alpha, " and ",
                    result$control$beta)
    }

    return(out)
  })
}

# Run the application
shinyApp(ui = ui, server = server)

