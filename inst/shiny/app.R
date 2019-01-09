library(dynamichazard)
if(!require(shiny))
  stop("Requires 'shiny' to run the demo")

test_sim_func_exp <- dynamichazard:::test_sim_func_exp
exp_model_names <- dynamichazard:::exp_model_names

# Global params
quietly <- if(exists("quietly")) quietly else FALSE
cat_l <- if(quietly) function(...) { invisible() } else cat
print_l <- if(quietly) function(...) { invisible() } else print

t_max <- 30
start_fun <- function(t_0 = t_0, t_max = t_max)
  max(0, runif(1, t_0 - t_max, t_max - 1 - 1e-8))

# Starting arguments
start_args <- list(
  n_series = 2, sim_with = "exponential", sim_fix_options = 1,
  obs_time = t_max, seed = 65848, covar_range = c(-.5, .5),
  sd_intercept = .2, sd_coef = .5, est_with_model = "exponential",
  est_with_method = "EKF", est_fix_options = 1, LR = 1,
  order = 1, denom_term = 1, fixed_terms_method = "M_step",
  use_extra_correction = F, beta = 0, alpha = 1,
  SMA_version = "woodbury", GMA_max_rep = 25,
  GMA_NR_eps = 2, more_options = F, debug = FALSE,
  n_threads = max(1, parallel::detectCores(logical = FALSE)))

if(exists("input_args")){
  if(any(is.na(arg_match <- match(names(input_args), names(start_args)))))
    stop("These input arguments are not recognized: ",
         paste0(names(input_args)[is.na(arg_match)], collapse = "\t"))

  start_args[arg_match] <- input_args
}

cat_l("Starting app with these arguments:\n")
print_l(start_args)

# Plot function to show density
breaks <-  function(model){
  if(model != "logit")
    return(seq(0, t_max, length.out = 101))

  return(0:t_max)
}

rug_dens <- function(case_times, control_count, model){
  breaks <- breaks(model)
  n_breaks <- length(breaks)
  case_count <- hist(case_times, plot = F, breaks = breaks)$counts

  b <- 2
  cases_cols_scale <- log((case_count + 1) / max(control_count), b)
  controls_cols_scale <- log((control_count + 1) / max(control_count), b)

  ma <- max(cases_cols_scale, controls_cols_scale)
  mi <- min(cases_cols_scale, controls_cols_scale)

  a <- .8
  cases_cols_scale <- a * (cases_cols_scale - mi) / (ma - mi)
  controls_cols_scale <- a * (controls_cols_scale - mi) / (ma - mi)

  # Add rect for controls
  y_consts <- c(.01, .03)
  y_dif <- diff(par("usr")[3:4])
  y0 <- par("usr")[3] - y_consts[1] * y_dif
  y1 <- y0 -  y_consts[2] * y_dif
  cols <- rgb(0, 0, 0, controls_cols_scale)

  rect(head(breaks, -1), rep(y0, length(breaks) - 1),
       breaks[-1], rep(y1, length(breaks) - 1),
       col = cols, border = NA, xpd = TRUE)

  # Add rect for cases
  y0 <- par("usr")[3] + y_consts[1] * y_dif
  y1 <- y0 +  y_consts[2] * y_dif
  cols <- rgb(0, 0, 0, cases_cols_scale)
  rect(head(breaks, -1), rep(y0, length(breaks) - 1),
       breaks[-1], rep(y1, length(breaks) - 1),
       col = cols, border = NA, xpd = TRUE)
}

# params for UI
col_w <- 4
col_w_out <- 4

get_em <- function(n_lines) paste0("height:",n_lines * 1.5,"em;")

# Thanks to this guy http://stackoverflow.com/a/31066997
get_JS_code_for_log_slider <- function(input_name, base, exp_min, exp_max){
  paste0("$(function() {
  setTimeout(function(){
  var vals = [", base^exp_min, "];
  var powStart = ", exp_min + 1,";
  var powStop = ", exp_max, ";
  for (i = powStart; i <= powStop; i++) {
  var val = Math.pow(", base,", i);
  val = parseFloat(val.toFixed(8));
  vals.push(val);
  }
  $('#", input_name,"').data('ionRangeSlider').update({'values':vals})
  }, 5)})")
}

n_series_stuff <- list(base = 2, exp_min = 6, exp_max = 19)

JScode <- get_JS_code_for_log_slider(
  "n_series", n_series_stuff$base, n_series_stuff$exp_min, n_series_stuff$exp_max)

denom_term_stuff <- list(base = 10, exp_min = -6, exp_max = -1)
GMA_NR_eps_stuff <- list(base = 10, exp_min = -6, exp_max = 1)

JScode <- paste(
  JScode,

  get_JS_code_for_log_slider(
  "denom_term", denom_term_stuff$base, denom_term_stuff$exp_min, denom_term_stuff$exp_max),

  get_JS_code_for_log_slider(
    "GMA_NR_eps", GMA_NR_eps_stuff$base, GMA_NR_eps_stuff$exp_min, GMA_NR_eps_stuff$exp_max),

  sep = "\n")

# Define UI for application that draws a histogram
ui <- fluidPage(

 tags$head(tags$script(HTML(JScode))),

 fluidRow(
   style = "width: 1280px;",
   column(
     9, offset = 3,
     titlePanel("ddhazard demo"))),

 fluidRow(
   style = "width: 1280px;",
   column(
     3,

     wellPanel(
       h4("Simulation settings"),

       sliderInput("n_series",
                   "Number of series to simulate",
                   min = 0,
                   max = n_series_stuff$exp_max - n_series_stuff$exp_min - 1,
                   value = start_args$n_series),

       selectInput("sim_with",
                   "Choose model to simulate from",
                   choices = c("logit", "cloglog", "exponential"),
                   selected = start_args$sim_with),

       radioButtons("sim_fix_options",
                    label = "Number of fixed coefficients",
                    choices = list("Zero" = 1,
                                   "Intercept and two coef" = 2,
                                   "All but one coef" = 3),
                    selected = start_args$sim_fix_options),

       conditionalPanel(
         "input.more_options",
         sliderInput("obs_time",
                     "Observed time",
                     min = 1,
                     max = t_max,
                     step = 1,
                     value = start_args$obs_time),

         numericInput("seed",
                      label = "RNG seed",
                      value = start_args$seed),

         sliderInput("covar_range", "Range to draw covariates from uniformly",
                     min = -2, max = 2, value = start_args$covar_range, step = .5),

         sliderInput("sd_intercept", "Standard deviations for intercept",
                     min = .1, max = 1, value = start_args$sd_intercept, step = .1),

         sliderInput("sd_coef", "Standard deviations for coefficients",
                     min = .1, max = 1, value = start_args$sd_coef, step = .1))
     ),

     wellPanel(
       h4("Estimation settings"),

       selectInput("est_with_model",
                   "Choose model to estimate with",
                   choices = c("logit", "cloglog", "exponential"),
                   selected = start_args$est_with_model),

       selectInput("est_with_method",
                   "Choose method to use in the E-step",
                   choices = c("UKF", "EKF", "SMA", "GMA"),
                   selected = start_args$est_with_method),

       radioButtons("est_fix_options",
                    label = "Number of fixed coefficients",
                    choices = list("Zero" = 1,
                                   "Intercept and two coef" = 2,
                                   "All but one coef" = 3),
                    selected = start_args$est_fix_options),

       conditionalPanel(
         "input.more_options",
         sliderInput("LR",
                     "Learning rate",
                     min = .1,
                     max = 2,
                     step = .1,
                     value = start_args$LR),

         sliderInput("order",
                     "Randowm walk order in estimation",
                     min = 1,
                     max = 2,
                     step = 1,
                     value = start_args$order),

         conditionalPanel(
           "input.est_with_method == 'EKF' || input.est_with_method == 'UKF'",
           sliderInput("denom_term",
                       "denom_term for extra term in denominators",
                       min = 0,
                       max = denom_term_stuff$exp_max - denom_term_stuff$exp_min - 1,
                       value = start_args$denom_term)),

         selectInput("fixed_terms_method",
                     "Estimate fixed effect in",
                     choices = c("E_step", "M_step"),
                     selected = start_args$fixed_terms_method),

         checkboxInput("debug",
                       "Print debug information to file",
                       value = start_args$debug))
     ),


     conditionalPanel(
       "input.more_options",
         wellPanel(
           conditionalPanel(
             "input.est_with_method == 'EKF'",
             h4("EKF settings"),

             checkboxInput("use_extra_correction",
                           "Extra correction steps",
                           value = start_args$use_extra_correction)),

           conditionalPanel(
             "input.est_with_method == 'UKF'",
             h4("UKF settings"),

             sliderInput("beta",
                         "Beta",
                         min = 0,
                         max = 2,
                         step = .5,
                         value = start_args$beta),

             sliderInput("alpha",
                         "Alpha",
                         min = 1e-2,
                         max = 1,
                         step = 1e-2,
                         value = start_args$alpha)),

           conditionalPanel(
             "input.est_with_method == 'SMA'",
             h4("SMA method settings"),

             selectInput("SMA_version",
                         "Computation version",
                         choices = c("woodbury", "cholesky"),
                         selected = start_args$SMA_version)),

           conditionalPanel(
             "input.est_with_method == 'GMA'",
             h4("GMA method settings"),

             sliderInput("GMA_max_rep",
                         "GMA_max_rep parameter",
                         min = 1,
                         max = 25,
                         step = 1,
                         value = start_args$GMA_max_rep),

             sliderInput("GMA_NR_eps",
                         "GMA_NR_eps parameter",
                         min = 0,
                         max = GMA_NR_eps_stuff$exp_max - GMA_NR_eps_stuff$exp_min - 1,
                         value = start_args$GMA_NR_eps))
         )),

     wellPanel(
       checkboxInput("more_options", label = "Show more options",
                     value = start_args$more_options, width = "12em"))

     ),

   column(
     12 - 3,
     plotOutput("coef_plot"),


     fluidRow(
       column(
         6,
         h3("Intro"),
         div("Illustrates simulated data and a fit. The true coefficients are the continous curves and the predicted coefficients are the dashed curves. Shaded areas are 95% confidence bounds from smoothed covariance matrices", style = get_em(4)),
         div(textOutput("rug_explanation"), style = get_em(4)),
         div("See the ddhazard vignette for further details", style = get_em(1))
       ),

       column(
         6,
         h3("Output"),
         div(textOutput("n_deaths"), style = get_em(1)),
         div(htmlOutput ("LR_out"), style = get_em(1)),
         div(textOutput("MSE"), style = get_em(1)),
         div(textOutput("model_text"), style = get_em(3))
       )
     ),

     fluidRow(
       column(
         6,
         div()
       ),

       column(
         6,
         h3("Fit call"),
         div(verbatimTextOutput("fit_call_txt"), style = "height:40em;")
       )
     ),

     h3("First 250 rows of simulated data"),
     dataTableOutput("sim_dat")
   )
 )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
  n_series_input <- reactive({
    n_series_stuff$base^(n_series_stuff$exp_min + input$n_series)
  })

  denom_term <- reactive({
    denom_term_stuff$base^(denom_term_stuff$exp_min + input$denom_term)
  })

  GMA_NR_eps <- reactive({
    GMA_NR_eps_stuff$base^(GMA_NR_eps_stuff$exp_min + input$GMA_NR_eps)
  })

  sim_input <- reactive({
    f_choice <- if(input$sim_with == "exponential")
      quote(test_sim_func_exp) else quote(test_sim_func_discrete)

    n_fixed <- if(input$sim_fix_options == 1)
      0 else if(input$sim_fix_options == 2)
        3 else if(input$sim_fix_options == 3)
          5

    n_varying <- 5 - n_fixed + (n_fixed > 0)

    x_range <- diff(input$covar_range)
    x_mean <- mean(input$covar_range)

    sim_quote <- bquote(
      .(f_choice)(
        n_series = .(n_series_input()),
        n_vars = 5,
        t_max = .(t_max), re_draw = T, beta_start = runif(5, min = -1.5, max = 1.5),
        intercept_start = -3.5,
        sds = c(.(input$sd_intercept), rep(.(input$sd_coef), 5)),
        x_range = .(x_range), x_mean = .(x_mean), lambda = .(5 / t_max),
        is_fixed = .(if(n_fixed == 0) c() else 1:n_fixed),

        tstart_sampl_func = .(start_fun)))
    if(input$sim_with == "cloglog")
      sim_quote[["linkfunc"]] <- "cloglog"
    else if(input$sim_with == "logit")
      sim_quote[["linkfunc"]] <- "logit"

    sim_exp <- bquote({
      set.seed(.(input$seed))
      dat <- .(sim_quote)

      dat$res[dat$res$event == 0 & dat$res$tstop > t_max, "tstop"] <- t_max
    })

    out <- formatR::tidy_source(text = capture.output(sim_exp), width.cutoff = 50,
                                output = FALSE)
    out <- out$text.tidy
    cat_l(paste0("Data is simulated with this call: \n", out, "\n"))

    eval(sim_exp)
    dat
  })

  output$sim_dat <- renderDataTable({
    dat <- sim_input()$res
    dat <- dat[seq_len(min(nrow(dat), 250)), ]
    dat[, c("x1", "x2", "x3", "x4", "x5")] <-
      round(dat[, c("x1", "x2", "x3", "x4", "x5")], 2)

    dat
  },
  options = list(pageLength = 10, searching = FALSE))

  n_fixed_when_est <- reactive({
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

  fit_quote_input <- reactive({
    sims <- sim_input()
    data <- sims$res

    debug <- input$debug

    n_fixed <- n_fixed_when_est()

    if(n_fixed == 0){
      form <- survival::Surv(tstart, tstop, event) ~ x1 + x2 + x3 + x4 + x5

    } else if(n_fixed == 3){
      form <- survival::Surv(tstart, tstop, event) ~ ddFixed_intercept() +
        ddFixed(x1) + ddFixed(x2) + x3 + x4 + x5

    } else if(n_fixed == 5){
      form <- survival::Surv(tstart, tstop, event) ~ ddFixed_intercept() +
        ddFixed(x1) + ddFixed(x2) + ddFixed(x3) + ddFixed(x4) + x5

    } else
      stop("n_fixed is not implemented")

    q_0_len <- 6 - n_fixed
    q_0_term <- ifelse(
      input$est_with_method %in% c("UKF", "GMA"),
      1, 1e5)

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

    control_list <- list(eps = 1e-3, n_max = 100,
                         method = input$est_with_method,
                         n_threads = start_args$n_threads)

    if(input$est_with_method == "UKF"){
      control_list <- c(control_list,
                        list(beta = input$beta, alpha = input$alpha,
                             denom_term = denom_term()))

    } else if (input$est_with_method == "EKF"){
      control_list$denom_term <- denom_term()
      if(input$use_extra_correction)
        control_list <- c(control_list,
                          list(NR_eps = 1e-5))

    } else if(input$est_with_method == "SMA"){
      control_list$posterior_version = input$SMA_version
    } else if(input$est_with_method == "GMA"){
      control_list <- c(control_list,
                        list(GMA_max_rep = input$GMA_max_rep,
                             GMA_NR_eps = GMA_NR_eps()))
    }

    if(n_fixed > 0){
      control_list <- c(control_list, list(fixed_terms_method  = input$fixed_terms_method))
    }

    if(input$LR != 1){
      control_list$LR <- input$LR
    }

    if(debug > 0){
      control_list$debug <- debug
    }

    control_list <- as.call(c(quote(ddhazard_control), control_list))

    Q <- if(6 - n_fixed > 1)
      bquote(diag(.1, .(6 - n_fixed))) else 0.1

    quote <- bquote(ddhazard(
      formula = .(form),
      data = data,
      by = 1,
      Q_0 = .(Q_0),
      Q = .(Q),
      max_T = .(input$obs_time),
      id = data$id,
      order = .(input$order),
      model = .(input$est_with_model),
      control = .(control_list)))

    if(debug > 0){
      tmp_file <- tempfile()
      tmp_file <- paste0(gsub("\\\\", "/", tmp_file), ".txt")
      quote <- bquote({
        # Open this file to get information about the estimation
        f <- file(.(tmp_file))
        sink(f)
        tryCatch(
          .(quote),
          finally = {
            sink()
            close(f)
          })
      })
    }

    list(quote = quote,
         data = data,
         debug_file = if(debug) tmp_file else NULL)
  })

  output$fit_call_txt <- renderText({
    eval_quote <- fit_quote_input()$quote
    out <- formatR::tidy_source(text = capture.output(eval_quote), width.cutoff = 30,
                                output = FALSE)
    out <- out$text.tidy

    cat_l(paste0("Model is fitted with this call: \n", out, "\n"))

    if(input$debug)
      out <- paste0(
        out, "\n\nCall file.edit('", fit_quote_input()$debug_file, "') to see debug info")

    paste0(out, "\n\n# data is the simulated data set")
  })

  fit_input <- reactive({
    tmp <- fit_quote_input()
    eval_quote <- tmp$quote
    data <- tmp$data
    t <- system.time(fit <- eval(eval_quote))

    list(fit = fit, time = t)
  })

  output$n_deaths <- renderText({
    sprintf("%d of %d dies", sum(sim_input()$res$event), n_series_input())
  })

  output$LR_out <- renderText({
    LR <- fit_input()$fit$LR
    out <- sprintf(paste0("The used learning rate is ", ifelse(LR == 1, "%d", "%.3f")), LR)

    if(LR < 1)
      out <- paste0("<font color=\"#FF0000\"><b>", out, "</b></font>")

    out
  })

  output$MSE <- renderText({
    n_fixed <- n_fixed_when_est()
    fit <- fit_input()$fit

    state_vecs <- if(n_fixed == 0)
      fit$state_vecs[1:input$obs_time + 1, 1:6] else
        cbind(sapply(fit$fixed_effects, rep, times = input$obs_time),
              fit$state_vecs[1:input$obs_time + 1, 1:(6 - n_fixed)])

    sprintf("MSE for coefficients is %.3f",
      mean.default((sim_input()$beta[1:input$obs_time + 1, ] - state_vecs)^2))
  })

  output$rug_explanation <- renderText({
    "The density above x-axis indicates when a death is observed. The density below the x-axis indicates the size of the control set. The densities are log scaled"
  })

  output$coef_plot <- renderPlot({
    sims <- sim_input()
    fit <- fit_input()$fit
    n_fixed <- n_fixed_when_est()

    cols <- palette()[1:ncol(sims$beta)]
    cols_conf <- apply(sapply(cols, col2rgb)/255, 2,
                       function(x)
                         rgb(x[1], x[2], x[3], alpha=.1))

    par(mai = rep(1, 4))
    x <- seq_len(dim(sims$beta)[1]) - 1
    matplot(x, sims$beta, lty = 1, type = "l",
            ylim = range(sims$beta, fit$state_vecs, fit$fixed_effects),
            ylab = "Coefficients", xlab = "Time", cex.lab = 1.4,
            frame = FALSE, axes=F, xlim = c(0 - .2, t_max + .2), xaxs = "i",
            col = cols)
    matplot(x, fit$state_vecs[, 1:(6 - n_fixed)],
            lty = 2, type = "l", add = T,
            col = cols[(1+n_fixed):6])

    for(i in (1+n_fixed):6){
      j <- i - n_fixed
      fac <- qnorm(.5 - .95 / 2)
      lb = fit$state_vecs[, j] + fac * sqrt(fit$state_vars[j, j, ])
      ub = fit$state_vecs[, j] - fac * sqrt(fit$state_vars[j, j, ])

      polygon(c(x, rev(x)), c(ub, rev(lb)),
              col = cols_conf[i], border = NA)
    }

    axis(1, lwd.ticks = 0, at = c(-10, seq(0, 30, 3), 100))
    axis(2, at = c(-10, axTicks(2), 10))

    if(n_fixed > 0)
      abline(h = fit$fixed_effects, col = cols[1:n_fixed], lty = 2)

    n_breaks <- length(breaks(input$sim_with))
    risk_obj <- with(
      sims$res,
      get_risk_obj(Surv(tstart, tstop, event), t_max / (n_breaks - 1), max_T = t_max,
                   id = id, is_for_discrete_model = fit$model == "logit"))
    count_at_risk <- unlist(lapply(risk_obj$risk_sets, length))

    # Add rug plots to illustrate survivers and deaths
    death_times <- sims$res$tstop[sims$res$event==1]

    rug_dens(death_times, count_at_risk,  input$sim_with)
  })

  output$model_text <- renderText({
    tmp <- fit_input()
    result <- tmp$fit

    out <- paste0("Model is estimated with order ", result$order, " and with the ", result$method, " in the E-step.",
                  " ", result$n_iter, " iterations was used in the EM algorithm.")

    if(result$method == "EKF"){
      out <- paste0(out, ifelse(is.null(result$control$NR_eps), " No extra", " Extra"),
                    " iterations are used in correction step.")
    } else if (result$method == "UKF"){
      out <- paste0(out, " Alpha and beta are ", result$control$alpha, " and ",
                    result$control$beta, ".")
    }

    t <- tmp$t
    fmt <- "%.3f"
    out <-
      paste0(out, "\n The computation time of estimation was",
             " user: ", sprintf(fmt, t["user.self"]),
             ", system: ", sprintf(fmt, t["sys.self"]),
             ", elapsed: ", sprintf(fmt, t["elapsed"]))

    return(out)
  })
}

# Run the application
shinyApp(ui = ui, server = server)
