#' @title Get data.frame for Discrete Time Survival Models
#' @description Function used to get \code{data.frame} with weights for a static fit for survivals.
#'
#' @inheritParams ddhazard
#' @param init_weights weights for the rows in \code{data}. Useful e.g., with skewed sampling.
#' @param risk_obj a pre-computed result from a \code{\link{get_risk_obj}}. Will be used to skip some computations.
#' @param use_weights \code{TRUE} if weights should be used. See details.
#' @param is_for_discrete_model \code{TRUE} if the model is for a discrete hazard model is used like the logistic model.
#' @param c_outcome,c_weights,c_end_t alternative names to use for the added columns described in the return section. Useful if you already have a column named \code{Y}, \code{t} or \code{weights}.
#'
#' @details
#' This function is used to get the \code{data.frame} for e.g. a \code{glm} fit that is comparable to a \code{\link{ddhazard}} fit in the sense that it is a static version. For example, say that we bin our time periods into \code{(0,1]}, \code{(1,2]} and \code{(2,3]}. Next, consider an individual who dies at time 2.5. He should be a control in the the first two bins and should be a case in the last bin. Thus the rows in the final data frame for this individual is \code{c(Y = 1, ..., weights = 1)} and \code{c(Y = 0, ..., weights = 2)} where \code{Y} is the outcome, \code{...} is the covariates and \code{weights} is the weights for the regression. Consider another individual who does not die and we observe him for all three periods. Thus, he will yield one row with \code{c(Y = 0, ..., weights = 3)}.
#'
#' This function use similar logic as the \code{ddhazard} for individuals with time varying covariates (see the vignette \code{vignette("ddhazard", "dynamichazard")} for details).
#'
#' If \code{use_weights = FALSE} then the two previously mentioned individuals will yield three rows each. The first individual will have \code{c(Y = 0, t = 1, ..., weights = 1)}, \code{c(Y = 0, t = 2, ..., weights = 1)}, \code{c(Y = 1, t = 3, ..., weights = 1)} while the latter will have three rows \code{c(Y = 0, t = 1, ..., weights = 1)}, \code{c(Y = 0, t = 2, ..., weights = 1)}, \code{c(Y = 0, t = 3, ..., weights = 1)}. This kind of data frame is useful if you want to make a fit with e.g. \code{\link[mgcv]{gam}} function in the \code{mgcv} package as described en Tutz et. al (2016).
#'
#' @return
#' Returns a \code{data.frame} where the following is added (column names will differ if you specified them): column \code{Y} for the binary outcome, column \code{weights} for weights of each row and additional rows if applicable. A column \code{t} is added for the stop time of the bin if \code{use_weights = FALSE}. An element \code{Y} with the used \code{Surv} object is added if \code{is_for_discrete_model = FALSE}.
#'
#' @seealso
#' \code{\link{ddhazard}}, \code{\link{static_glm}}
#'
#' @references
#' Tutz, Gerhard, and Matthias Schmid. \emph{Nonparametric Modeling and Smooth Effects}. Modeling Discrete Time-to-Event Data. Springer International Publishing, 2016. 105-127.
#'
#' @examples
#'library(dynamichazard)
#'# small toy example with time-varying covariates
#'dat <- data.frame(
#'  id     = c(   1,    1, 2,     2),
#'  tstart = c(   0,    4, 0,     2),
#'  tstop  = c(   4,    6, 2,     6),
#'  event  = c(   0,    1, 0,     0),
#'  x1     = c(1.09, 1.29, 0, -1.16))
#'
#'get_survival_case_weights_and_data(
#'  Surv(tstart, tstop, event) ~ x1, dat, by = 1, id = dat$id)$X
#'get_survival_case_weights_and_data(
#'  Surv(tstart, tstop, event) ~ x1, dat, by = 1, id = dat$id,
#'  use_weights = FALSE)$X
#'
#' @export
get_survival_case_weights_and_data <- function(
  formula, data, by, max_T, id, init_weights, risk_obj,
  use_weights = T, is_for_discrete_model = T,
  c_outcome = "Y",
  c_weights = "weights",
  c_end_t = "t"){
  #####
  # checks

  X_Y <- get_design_matrix(formula, data, predictors = F)
  formula_used <- X_Y$formula_used

  if(missing(init_weights))
    init_weights <- rep(1, nrow(data))

  c_outcome <- change_new_var_name(c_outcome, data = data)
  c_weights <- change_new_var_name(c_weights, data = data)
  c_end_t   <- change_new_var_name(c_end_t, data = data)

  #####
  # find data frame and return

  # TODO: an expression is saved here as other code have passed a missing max_T
  #       argument along where we do not need to evaluate max_T. Re-solve the
  #       the issue higher up and make the evalutation here regardless of
  #       whether we need max_T to simplify the code.
  set_max_T_expr <- expression(
    max_T <- ifelse(
      missing(max_T), min(
        max(X_Y$Y[X_Y$Y[, 3] == 1, 2]),
        max(X_Y$Y[X_Y$Y[, 3] == 0, 2])),
      max_T))

  out <- if(is_for_discrete_model){
    compute_risk_obj <- missing(risk_obj) || is.null(risk_obj)

    if(missing(id) && compute_risk_obj){
      warning("You did not parse and ID argument. This can affact result for discrete models with time-coefficients effects")
      id = 1:nrow(data)
    }

    if(compute_risk_obj)
      eval(set_max_T_expr)

    .get_survival_case_weights_and_data_discrete(
      risk_obj, compute_risk_obj, data, X_Y, by, max_T, use_weights, c_outcome,
      c_end_t, c_weights, id, formula_used, init_weights)
  } else{
    eval(set_max_T_expr)

    .get_survival_case_weights_and_data_continous(
      X_Y, max_T, data, formula_used, c_outcome)
  }

  return(out)
}

#' @importFrom utils tail
change_new_var_name <- function(c_x, data){
  if(any(names(data) == c_x)){
    c_x_new <- tail(make.unique(c(names(data), c_x)), 1)
    warning("Column called ", sQuote(c_x), " is already in the data.frame. Will use", sQuote(c_x_new), "instead.")

    c_x_new
  } else
    c_x
}

#' @importFrom utils tail
.get_survival_case_weights_and_data_discrete <- function(
  risk_obj, compute_risk_obj, data, X_Y, by, max_T, use_weights, c_outcome,
  c_end_t, c_weights, id, formula_used, init_weights){
  if(!compute_risk_obj && nrow(data) <
     max(unlist(lapply(risk_obj$risk_sets, max))))
    stop("risk_obj has indicies out site of data. Likely the risk_set comes from a different data set")

  if(compute_risk_obj){
    risk_obj <- get_risk_obj(Y = X_Y$Y, by = by, max_T = max_T, id = id,
                             is_for_discrete_model = TRUE)

    risk_obj$event_times <- tail(risk_obj$event_times, -1)
  }

  if(use_weights){
    new_weights = rep(0, nrow(data))
    new_case_rows <- NULL

    for(i in seq_along(risk_obj$risk_sets)){
      time_ = risk_obj$event_times[i]
      r_set = risk_obj$risk_sets[[i]]

      is_case = risk_obj$is_event_in[r_set] == i - 1
      r_set_is_case <- r_set[is_case]
      r_set_not_case <- r_set[!is_case]

      new_case_rows <- c(
        new_case_rows, list(c(
          list(Y = rep(1, sum(is_case))),
          data[r_set_is_case, ],
          list(weights = init_weights[r_set_is_case]))))

      new_weights[r_set_not_case] <-
        new_weights[r_set_not_case] + init_weights[r_set_not_case]
    }

    do_keep <- new_weights > 0
    n_keep <- sum(do_keep)
    do_keep <- which(do_keep)
    X = c(
      list(Y = rep(0, n_keep)),
      data[do_keep, ],
      list(weights = new_weights[do_keep]))
    X <- .rbind_list(c(list(X), new_case_rows))
    names(X)[c(1, ncol(X))] <- c(c_outcome, c_weights)

  } else {
    # TODO: change to more general setup that does not convert all columns to
    #       double values first
    n_rows_final <- sum(unlist(lapply(risk_obj$risk_sets, length)))
    X <- data.frame(matrix(NA, nrow = n_rows_final, ncol = 3 + ncol(data),
                           dimnames = list(
                             NULL, c(colnames(data),
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

.get_survival_case_weights_and_data_continous <- function(
  X_Y, max_T, data, formula_used, c_outcome){
  # round times if needed
  Y <- X_Y$Y
  rm(X_Y)
  for(i in 1:2){
    ord <- order(Y[, i]) - 1L
    Y[, i] <- round_if_almost_eq(Y[, i], ord, max_T)
  }

  # keep only observations before max_T
  is_before_max_T <- Y[, 1] < max_T
  data <- data[is_before_max_T, ]
  Y <- Y[is_before_max_T, ]

  X <- eval(parse(text = paste0(
    "cbind(", c_outcome, " = Y[, 3] & Y[, 2] <= max_T, data)")))

  list(X = X, formula_used = formula_used, Y = Y)
}


#' @title  Static glm Fit
#' @inheritParams ddhazard
#' @inheritParams get_survival_case_weights_and_data
#' @param ... arguments passed to \code{\link{glm}} or \code{\link[speedglm]{speedglm}}. If \code{only_coef = TRUE} then the arguments are passed to \code{\link{glm.control}} if \code{\link{glm}} is used.
#' @param family \code{"logit"}, \code{"cloglog"}, or \code{"exponential"} for a static equivalent model of \code{\link{ddhazard}}.
#' @param model \code{TRUE} if you want to save the design matrix used in \code{\link{glm}}.
#' @param speedglm depreciated.
#' @param only_coef \code{TRUE} if only coefficients should be returned. This will only call the \code{speedglm::speedglm.wfit} or \code{\link{glm.fit}} which will be faster.
#' @param mf model matrix for regression. Needed when \code{only_coef = TRUE}
#' @param method_use method to use for estimation. \code{\link{glm}} uses \code{\link{glm.fit}}, \code{speedglm::speedglm} uses \code{speedglm::speedglm.wfit} and \code{parallelglm_quick} and \code{parallelglm_QR} uses a parallel \code{C++} estimation method.
#' @param n_threads number of threads to use when \code{method_use} is \code{"parallelglm"}.
#'
#' @description
#' Method to fit a static model corresponding to a \code{\link{ddhazard}} fit.
#' The method uses weights to ease the memory requirements. See
#' \code{\link{get_survival_case_weights_and_data}} for details on weights.
#'
#' The \code{parallelglm_quick} and \code{parallelglm_QR} methods are similar
#' to two methods used in \code{bam} function in the \code{mgcv} package (see
#' the \code{`use.chol`} argument or Wood et al. 2015). \code{parallelglm_QR}
#' is more stable but slower. See Golub (2013) section 5.3 for a comparison of
#' the Cholesky decomposition method and the QR method.
#'
#' @return
#' The returned list from the \code{\link{glm}} call or just coefficients depending on the value of \code{only_coef}.
#'
#' @references
#' Wood, S.N., Goude, Y. & Shaw S. (2015) Generalized additive models for large datasets. Journal of the Royal Statistical Society, Series C 64(1): 139-155.
#'
#' Golub, G. H., & Van Loan, C. F. (2013). Matrix computations (4th ed.). JHU Press.
#'
#' @examples
#'library(dynamichazard)
#'fit <- static_glm(
#'  Surv(time, status == 2) ~ log(bili), pbc, id = pbc$id, max_T = 3600,
#'  by = 50)
#'fit$coefficients
#'
#'
#' @export
static_glm <- function(
  formula, data, by, max_T, ..., id, family = "logit", model = F, weights,
  risk_obj = NULL, speedglm = F, only_coef = FALSE, mf,
  method_use = c("glm", "speedglm", "parallelglm_quick", "parallelglm_QR"),
  n_threads = getOption("ddhazard_max_threads")){
  if(only_coef && missing(mf))
    stop("mf must be supplied when ", sQuote("only_coef = TRUE"))

  if(!missing(mf) && nrow(mf) != nrow(data))
    stop("data and mf must have the same number of rows")

  if(speedglm)
    warning(sQuote("speedglm"), " have been depreciated. Use ",
            sQuote("method_use"))

  method_use <- method_use[1]

  formula_org <- formula

  if(only_coef){
    # we mark the row numbers as some may be removed and the order may be
    # changed
    col_for_row_n <- ncol(data) + 1
    data[[col_for_row_n]] <- 1:nrow(data)
    col_for_row_n <- colnames(data)[col_for_row_n]
  }

  # find new column variable names if needed
  c_outcome <- change_new_var_name("Y", data = data)
  c_weights <- change_new_var_name("weights", data = data)
  c_end_t   <- change_new_var_name("t", data = data)

  if(!missing(mf)){
    # we only need the outcome variable and weights from the next part
    formula <- update(formula, . ~ 1)
    data <- data[, c(all.vars(formula), col_for_row_n)]
  }

  if(family %in% c("binomial", "logit", "cloglog")){
    family <- binomial(switch(
      family, binomial = "logit", logit = "logit", cloglog = "cloglog"))

    tmp <- get_survival_case_weights_and_data(
      formula = formula, data = data, by = by, max_T = max_T, id = id,
      init_weights = weights, risk_obj = risk_obj,
      c_outcome = c_outcome,
      c_weights = c_weights,
      c_end_t = c_end_t)

    formula <- tmp$formula_used
    X <- tmp$X
    rm(tmp)

    data <- X
    formula <- eval(substitute(
      update(formula, Y ~ ., data = data),
      list(Y = as.name(c_outcome))))

  } else if(family == "exponential"){
    family <- poisson()

    tmp <- get_survival_case_weights_and_data(
      formula = formula, data = data, by = by, max_T = max_T, id = id,
      init_weights = weights, risk_obj = risk_obj,
      is_for_discrete_model = FALSE,
      c_outcome = c_outcome,
      c_weights = c_weights,
      c_end_t = c_end_t)

    formula <- tmp$formula_used
    data <- tmp$X
    Y <- tmp$Y
    rm(tmp)

    data[[c_weights]] <- if(missing(weights))
      rep(1, nrow(data)) else weights
    data[["log_delta_time"]] <- log(pmin(Y[, 2], max_T) - Y[, 1])

    eval(substitute(
      formula <- update(formula, Y ~ . + offset(log_delta_time), data = data),
      list(Y = as.name(c_outcome))))

  } else
    stop("family ", sQuote("family"), " not implemented in static_glm")

  if(only_coef){
    new_order <- data[[col_for_row_n]]
    mf <- mf[new_order, , drop = FALSE]
  }

  offset <- if(family$family == "poisson")
    data$log_delta_time else rep(0, nrow(data))

  cl <- match.call()
  cl <- cl[!names(cl) %in% names(formals(static_glm))]
  cl[[1L]] <- quote(.static_glm_fit)

  cl[c("method_use", "only_coef", "c_outcome", "c_weights", "family",
       "offset", "data", "mf", "n_threads", "formula", "model")] <- list(
         quote(method_use), quote(only_coef), quote(c_outcome),
         quote(c_weights), quote(family), quote(offset), quote(data),
         quote(mf), quote(n_threads), quote(formula), quote(model))

  eval(cl, environment())
}

.static_glm_fit <- function(
  method_use, only_coef, c_outcome, c_weights, family, offset, data, mf,
  n_threads, formula, model, ...){
  if(method_use == "speedglm" && requireNamespace("speedglm", quietly = T)){
    if(only_coef){
      fit <- speedglm::speedglm.wfit(
        X = mf, y = data[[c_outcome]], weights = data[[c_weights]],
        family = family, offset = offset, ...)

      return(fit$coefficients)
    }

    return(drop(speedglm::speedglm(
      formula = formula, data = data, family = family, model = model,
      weights = data[[c_weights]], ...)))

  } else if(method_use == "glm"){
    if(only_coef){
      ctrl <- do.call(glm.control, list(...))

      fit <- eval(bquote(
        glm.fit(
          x = mf, y = data[[.(c_outcome)]], weights = data[[.(c_weights)]],
          family = family, control = .(ctrl), offset = offset)))

      return(fit$coefficients)
    }

    return(eval(bquote(
      glm(formula = .(formula), data = data,
          family = family, model = model,
          weights = .(data[[c_weights]]), ...))))

  } else if(grepl("^parallelglm", method_use) && only_coef){
    epsilon <- list(...)$epsilon
    if(is.null(epsilon))
      epsilon <- glm.control()$epsilon

    method. <- gsub("(^parallelglm_)([a-zA-Z]+$)", "\\2", method_use)

    out <- drop(
      parallelglm(
        X = t(mf), Ys = data[[c_outcome]],
        weights = data[[c_weights]], offsets = offset, beta0 = numeric(),
        family = family$family, method = method.,
        tol = epsilon, nthreads = n_threads))

    return(structure(out, names = dimnames(mf)[[2]]))

  }

  stop(sQuote("method_use"), " not implemented with ", sQuote("only_coef"),
       " = ", only_coef)
}
