#' @title Static GLM fit for survival models
#' @description Function used to get design matrix and weights for a static fit for survivals models where observations are binned into intervals
#'
#' @param formula \code{\link[survival]{coxph}} like formula with \code{\link[survival]{Surv}(tstart, tstop, event)} on the left hand site of \code{~}
#' @param data Data frame or environment containing the outcome and co-variates
#' @param by Length of each intervals that cases are binned into
#' @param max_T The end time of the last bin
#' @param id The id for each row in \code{data}. This is important when variables are time varying
#' @param init_weights Weights for the rows \code{data}. Useful with skewed sampling and will be used when computing the final weights
#' @param risk_obj A pre-computed result from a \code{\link{get_risk_obj}}. Will be used to skip some computations
#' @param use_weights \code{TRUE} if weights should be used. See details
#' @param is_for_discrete_model \code{TRUE} if the model is for a discrete hazard model like the logistic model. Affects how deaths are included when individuals have time varying coefficients
#' @param c_outcome,c_weights,c_end_t Alternative names to use for the added columns described in the return section. Useful if you already have a column named \code{Y}, \code{t} or \code{weights}
#' @details
#' This function is used to get the data frame for e.g. a \code{glm} fit that is comparable to a \code{\link{ddhazard}} fit in the sense that it is a static version. For example, say that we bin our time periods into \code{(0,1]}, \code{(1,2]} and \code{(2,3]}. Next, consider an individual who dies at time 2.5. He should be a control in the the first two bins and should be a case in the last bin. Thus the rows in the final data frame for this individual is \code{c(Y = 1, ..., weights = 1)} and \code{c(Y = 0, ..., weights = 2)} where \code{Y} is the outcome, \code{...} is the co-variates and \code{weights} is the weights for the regression. Consider another individual who does not die and we observe him for all three periods. Thus, he will yield one row with \code{c(Y = 0, ..., weights = 3)}
#'
#' This function use similar logic as the \code{ddhazard} for individuals with time varying co-variates (see the vignette "ddhazard" for details)
#'
#' If \code{use_weights = FALSE} then the two individuals will yield three rows each. The first individual will have \code{c(Y = 0, t = 1, ..., weights = 1)}, \code{c(Y = 0, t = 2, ..., weights = 1)}, \code{c(Y = 1, t = 3, ..., weights = 1)} while the latter will have three rows \code{c(Y = 0, t = 1, ..., weights = 1)}, \code{c(Y = 0, t = 2, ..., weights = 1)}, \code{c(Y = 0, t = 3, ..., weights = 1)}. This kind of data frame is useful if you want to make a fit with e.g. \code{\link[mgcv]{gam}} function in the \code{mgcv} package as described en Tutz et. al (2016) (see reference)
#'
#' @return
#' Returns a data frame with the design matrix from the formula where the following is added (column names will differ if you specified them): column \code{Y} for the binary outcome, column \code{weights} for weights of each row and additional rows if applicable. A column \code{t} is added for the stop time of the bin if \code{use_weights = FALSE}
#'
#' @seealso
#' \code{\link{ddhazard}}, \code{\link{static_glm}}
#'
#' @references
#' Tutz, Gerhard, and Matthias Schmid. \emph{Nonparametric Modeling and Smooth Effects}. Modeling Discrete Time-to-Event Data. Springer International Publishing, 2016. 105-127.
#'
#' @export
get_survival_case_weights_and_data = function(
  formula, data, by, max_T, id, init_weights, risk_obj,
  use_weights = T, is_for_discrete_model = T,
  c_outcome = "Y",
  c_weights = "weights",
  c_end_t = "t"){
  X_Y = get_design_matrix(formula, data, predictors = F)
  formula_used <- X_Y$formula_used
  attr(data, "class") <- c("data.table", attr(data, "class"))

  if(missing(init_weights))
    init_weights = rep(1, nrow(data))

  if(any(colnames(data) == c_outcome)){
    warning("Column called '", c_outcome, "' will be replaced")
    data = data[, colnames(data) != c_outcome]
  }

  if(any(colnames(data) == c_weights)){
    warning("Column called '", c_weights, "' will be replaced")
    data = data[, colnames(data) != c_weights]
  }

  if(!use_weights && any(colnames(data) == c_end_t)){
    warning("Column called '", c_end_t, "' will be replaced")
    data = data[, colnames(data) != c_end_t]
  }

  compute_risk_obj <- missing(risk_obj) || is.null(risk_obj)

  if(!compute_risk_obj && nrow(data) <
     max(unlist(lapply(risk_obj$risk_sets, max))))
    stop("risk_obj has indicies out site of data. Likely the risk_set comes from a different data set")

  if(compute_risk_obj){
    if(missing(id)){
      warning("You did not parse and ID argument. This can affact result for discrete models with time-coefficients effects")
      id = 1:nrow(data)
    }

    risk_obj <- get_risk_obj(
      Y = X_Y$Y, by = by,
      max_T = ifelse(missing(max_T), max(X_Y$Y[X_Y$Y[, 3] == 1, 2]), max_T),
      id = id, is_for_discrete_model = is_for_discrete_model)

    risk_obj$event_times <- tail(risk_obj$event_times, -1)
  }

  if(use_weights){
    new_weights = rep(0, nrow(data))
    new_case_rows = list()

    for(i in seq_along(risk_obj$risk_sets)){
      time_ = risk_obj$event_times[i]
      r_set = risk_obj$risk_sets[[i]]

      is_case = risk_obj$is_event_in[r_set] == i - 1
      r_set_is_case <- r_set[is_case]
      r_set_not_case <- r_set[!is_case]

      new_case_rows <- c(
        new_case_rows, list(c(
          list(Y = rep(1, sum(is_case))),
          data[r_set_is_case],
          list(weights = init_weights[r_set_is_case]))))

      new_weights[r_set_not_case] = new_weights[r_set_not_case] + init_weights[r_set_not_case]
    }

    do_keep <- new_weights > 0
    n_keep <- sum(do_keep)
    do_keep <- which(do_keep)
    X = c(
      list(Y = rep(0, n_keep)),
      data[do_keep],
      list(weights = new_weights[do_keep]))
    X <- rbindlist(c(list(X), new_case_rows))
    attr(X, "colnames")[c(1, ncol(X))] <- c(c_outcome, c_weights)
    class(X) <- "data.frame"

  } else {
    n_rows_final <- sum(unlist(lapply(risk_obj$risk_sets, length)))
    X <- data.frame(matrix(NA, nrow = n_rows_final, ncol = 3 + ncol(data),
                           dimnames = list(NULL, c(colnames(data),
                                                   c_outcome, c_end_t, c_weights))))

    j <- 1
    for(i in seq_along(risk_obj$risk_sets)){
      time_ = risk_obj$event_times[i]
      r_set = risk_obj$risk_sets[[i]]

      n_risk <- length(r_set)
      is_case = risk_obj$is_event_in[r_set] == i - 1

      X[j:(j + n_risk -1), ] <- cbind(
        data[r_set, ],
        is_case, rep(time_, n_risk), rep(1, n_risk))

      j <- j + n_risk
    }

    for(i in which(sapply(data, is.factor)))
      X[[i]] <- factor(levels(data[[i]])[X[[i]]], levels = levels(data[[i]]))
  }

  list(X = X, formula_used = formula_used)
}


#' Function to make a static glm fit
#' @inheritParams get_survival_case_weights_and_data
#' @param ... arguments passed to \code{\link{glm}} or \code{\link[speedglm]{speedglm}}. If \code{only_coef = TRUE} then the arguments are passed to \code{\link{glm.control}} if \code{\link{glm}} is used
#' @param family \code{"logit"} or \code{"exponential"} for the static equivalent model of \code{\link{ddhazard}}
#' @param model \code{TRUE} if you want to save the design matrix used in \code{\link{glm}}
#' @param weights weights if a skewed sample or similar is used
#' @param speedglm Depreciated.
#' @param only_coef \code{TRUE} if only coefficients should be returned. This will only call the \code{\link[speedglm]{speedglm.wfit}} or \code{\link{glm.fit}} which will be faster.
#' @param mf model matrix for regression. Needed when \code{only_coef = TRUE}
#' @param method_use method to use for estimation. \code{\link{glm}} uses \code{\link{glm.fit}}, \code{\link[speedglm]{speedglm}} uses \code{\link[speedglm]{speedglm.wfit}} and \code{parallelglm} uses a parallel \code{C++} version \code{\link{glm.fit}} which only gives the coefficients.
#' @param n_threads number of threads to use when \code{method_use} is \code{"parallelglm"}.
#'
#' @details
#' Method to fit a static model corresponding to a \code{\link{ddhazard}} fit. The method uses weights to ease the memory requirements. See \code{\link{get_survival_case_weights_and_data}} for details on weights
#'
#' @return
#' The returned list from the \code{\link{glm}} call or just coefficients depending on the value of \code{only_coef}
#'
#' @export
static_glm = function(
  formula, data, by, max_T, ..., id, family = "logit", model = F, weights, risk_obj = NULL,
  speedglm = F, only_coef = FALSE, mf, method_use = c("glm", "speedglm", "parallelglm"),
  n_threads = getOption("ddhazard_max_threads")){
  if(only_coef && missing(mf))
    stop("mf must be supplied when only_coef = TRUE")

  if(!missing(mf) && nrow(mf) != nrow(data))
    stop("data and mf must have the same number of rows")

  if(speedglm)
    warning(sQuote("speedglm"), " have been depreciated. Use ", sQuote("method_use"))

  if(length(method_use) > 1)
    method_use <- method_use[1]

  if(only_coef){
    # We mark the row numbers as some may be removed and the order may be
    # changed
    col_for_row_n <- ncol(data) + 1
    data[[col_for_row_n]] <- 1:nrow(data)
    col_for_row_n <- colnames(data)[col_for_row_n]
  }

  if(family %in% c("binomial", "logit")){
    family <- binomial()

    tmp = get_survival_case_weights_and_data(
      formula = formula, data = data, by = by, max_T = max_T, id = id,
      init_weights = weights, risk_obj = risk_obj)

    formula <- tmp$formula_used
    X <- tmp$X
    rm(tmp)

    data <- X
    formula <- update(formula, Y ~ ., data = data)

  } else if(family == "exponential"){
    family <- poisson()
    # TODO: can be quicker when we just want the outcome
    X_Y = get_design_matrix(formula, data)
    X_Y$X <- X_Y$X[, -1] # remove the intercept

    formula <- X_Y$formula_used

    is_before_max_T <- X_Y$Y[, 1] < max_T
    data <- data[is_before_max_T, ]
    X_Y$Y <- X_Y$Y[is_before_max_T, ]

    X <- cbind(Y = X_Y$Y[, 3] & (X_Y$Y[, 2] <= max_T),
               data,
               log_delta_time = log(pmin(X_Y$Y[, 2], max_T) - X_Y$Y[, 1]),
               weights = rep(1, nrow(data)))

    data <- X
    formula <- update(formula, Y ~ . + offset(log_delta_time), data = data)

  } else
    stop("family '", family, "' not implemented in static_glm")

  if(only_coef){
    new_order <- data[[col_for_row_n]]
    mf <- mf[new_order, , drop = FALSE]
  }

  offset <- if(family$family == "poisson")
    data$log_delta_time else rep(0, nrow(data))

  if(method_use == "speedglm" && requireNamespace("speedglm", quietly = T)){
    if(only_coef){


      fit <- speedglm::speedglm.wfit(
        X = mf, y = data$Y, weights = data$weights,
        family = family, offset = offset, ...)

      return(fit$coefficients)
    }

    return(drop(speedglm::speedglm(
      formula = formula, data = data, family = family, model = model,
      weights = data$weights, ...)))

  } else if(method_use == "glm"){
    if(only_coef){
      ctrl <- do.call(glm.control, list(...))

      fit <- eval(bquote(
        glm.fit(x = mf, y = data$Y, weights = data$weights,
                family = .(family), control = .(ctrl), offset = offset)))

      return(fit$coefficients)
    }

    return(eval(bquote(
      glm(formula = .(formula), data = data,
          family = .(family), model = model,
          weights = weights, ...))))

  } else if(method_use == "parallelglm" && only_coef){
    epsilon <- list(...)$epsilon
    if(is.null(epsilon))
      epsilon <- glm.control()$epsilon

    out <- drop(
      parallelglm(X = t(mf), Ys = data$Y,
                  weights = data$weights, offsets = offset, beta0 = numeric(),
                  family = family$family,
                  tol = epsilon, nthreads = n_threads))

    return(structure(out, names = dimnames(mf)[[2]]))

  } else
    stop(sQuote("method_use"), " not implemented with ", sQuote("only_coef"),
         " = ", only_coef)
}
