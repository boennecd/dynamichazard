# Find the design matrix and returns it left hand site and right hand site of
# regression equation
get_design_matrix = function(formula, data, response = T, predictors = T){
  if(!response && !predictors)
    stop("get_design_matrix is called where neither respone or predictors is requested")

  Call <- match.call()
  indx <- match(c("formula", "data"),  names(Call), nomatch = 0)

  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")

  # We may call this expression twice. Thus, we save it unevaluted
  get_form_exp <-
      quote(terms(formula, data = data, specials = ddfixed_specials))

  temp$formula <- eval(get_form_exp)

  # Check if we have a fixed intercept
  # We want the following to give the same:
  # ... ~ ddFixed(1) + ...
  # ... ~ -1 + ddFixed(1) + ...
  # ... ~ ddFixed(rep(1, nrow(data))) + ...
  # ... ~ -1 + ddFixed(rep(1, nrow(data))) + ...
  # ... ~ ddFixed_intercept(nrow(data)) + ...
  # ... ~ -1 + ddFixed_intercept(nrow(data)) + ...

  is_fixed <- unlist(attr(temp$formula, "specials")[ddfixed_specials])
  if(length(is_fixed) > 0){
    fixed_label <- dimnames(attr(temp$formula, "factors"))[[1]][is_fixed]
    inter_match_regexp <- "(ddFixed\\(\\s*1\\s*\\))|(ddFixed\\(\\s*rep\\(\\s*1,)|(ddFixed_intercept\\()"
    is_fixed_intercept <- which(grepl(inter_match_regexp, fixed_label, perl = T))

    if(length(is_fixed_intercept) > 1)
      stop("Found more than one term that match these regexp patterns '",
           is_fixed_intercept, "'. These are used to detect fixed intercept. There should at most be one match")

    if(length(is_fixed_intercept) == 1){
      # We assume that future that methods that is going to extract the data
      # the formula object is going to use a data.frame or list and call
      # model.frame
      # Then, ls()[1] will give the name of the first column (if a data.frame)
      # is used in the line:
      #   variables <- eval(predvars, data, env)
      formula <- eval(parse(text = paste0(
        "update(formula, .~ ",
        " + ddFixed_intercept(length(eval(parse(text = ls()[1]))))",
        " + . -1 - ", fixed_label[is_fixed_intercept], ")")))

      # have to do this again
      temp$formula <- eval(get_form_exp)
    }}

  if(!response)
    temp$formula <- delete.response(temp$formula)

  class(temp$formula) <- c("ddformula", class(temp$formula))

  if(!predictors){
    old_form <- temp$formula
    temp$formula <- update(old_form, . ~ -1)
    mf <- eval(temp, parent.frame())

    Y <- model.response(mf)
    temp$formula <- old_form

    X <- NULL
    fixed_terms <- NULL

  } else{
    mf <- eval(temp, parent.frame())

    Y <- if(response) model.response(mf) else NULL

    if(!is.null(Y) && attr(Y, "type") == "right"){ # Change outcome if formula was Surv(stop, event) ~
      Y <- cbind(rep(0, nrow(Y)), Y)
      dimnames(Y)[[2]] <- c("start", "stop", "status")
      attr(Y, "type") <- "counting"
    }

    Terms <- terms(mf)
    fixed_terms_indicies <- unlist(attr(Terms, "specials")[ddfixed_specials])

    X <- model.matrix(Terms, mf)
    if(is.null(fixed_terms_indicies)){
      fixed_terms_indicies <- c()
    } else {
      fixed_terms_indicies <- which(
        attr(X, "assign") %in% (fixed_terms_indicies - attr(Terms, "response")))
    }

    fixed_terms <- X[, fixed_terms_indicies, drop = F]
    if(length(fixed_terms_indicies) > 0)
      X <- X[, -fixed_terms_indicies, drop = F]
  }

  # Change fixed_terms name if ddFixed(1) or similar was used
  is_ddFixed_intercept <- grepl("^ddFixed_intercept\\(.+\\)1$", colnames(fixed_terms))
  if(any(is_ddFixed_intercept)){
    old_name <- colnames(fixed_terms)[is_ddFixed_intercept]
    new_name <- "ddFixed((Intercept))"

    attr(temp$formula, "term.labels")[
      attr(temp$formula, "term.labels") == old_name] <- new_name
    colnames(attr(temp$formula, "factors"))[
      colnames(attr(temp$formula, "factors")) == old_name] <- new_name
    rownames(attr(temp$formula, "factors"))[
      rownames(attr(temp$formula, "factors")) == old_name] <- new_name

    if(predictors)
      colnames(fixed_terms)[is_ddFixed_intercept] <- new_name
  }

  list(X = X, fixed_terms = fixed_terms, Y = Y, formula = temp$formula)
}

#' model.frame and model.matrix for ddformula
#' @param formula Same as \code{\link{model.frame.default}}
#' @param subset Same as \code{\link{model.frame.default}}
#' @param na.action Same as \code{\link{model.frame.default}}
#' @param xlev Same as \code{\link{model.frame.default}}
#' @param object Same as \code{\link{model.matrix.default}}
#' @param data Same as \code{\link{model.matrix.default}}
#' @param contrasts.arg Same as \code{\link{model.matrix.default}}
#' @param xlev Same as \code{\link{model.matrix.default}}
#'
#' @description
#' Functions added to handle fixed (time-invariant) intercept coefficient for \code{\link{ddhazard}}. \code{model.frame.ddformula} always has \code{drop.unused.levels = FALSE} regardless of the input.
#'
#' @export
model.frame.ddformula <- function (
  formula, data = NULL, subset = NULL, na.action = na.fail,
  xlev = NULL, ...){
  call <- match.call()
  call$drop.unused.levels = FALSE
  call[[1]] <- stats::model.frame.default
  eval(call, parent.frame())
}

#' @rdname model.frame.ddformula
#' @export
model.matrix.ddformula  <- function(
  object, data = environment(object), contrasts.arg = NULL,
  xlev = NULL, ...){
  ans <- stats::model.matrix.default(object, data, contrasts.arg, xlevl)

  fix_inter <- attr(object, "specials")$ddFixed_intercept
  if(is.null(fix_inter) || length(fix_inter) == 0)
    return(ans)

  if(length(fix_inter) > 1)
    stop("Cannot specify more than one intercept term with ddFixed_intercept")

  if(attr(object, "intercept") > 0)
    stop("Cannot have intercept with formula when ddFixed_intercept is used. Repeat the call with ~ -1 + ... pr ~ 0 + ...")

  # Remove the zero level for the intercept
  is_inter_term <- which(attr(ans, "assign") == (fix_inter - attr(object, "response")))

  keep <- 1:ncol(ans) != is_inter_term[1]
  attri <- attributes(ans)
  attri$assign <- attri$assign[keep]
  attri$dim[2] <- attri$dim[2] - 1
  attri$dimnames[[2]] <- attri$dimnames[[2]][keep]
  ans <- ans[, keep]
  attributes(ans) <- attri

  ans
}

#' Auxiliary functions for fixed effects
#' @description
#' Functions used in formula of \code{\link{ddhazard}} for time-invariant effects. \code{ddFixed_intercept} is only used for the intercept.
#' @param object Expression that would be used in formula. E.g. \code{x} or \code{poly(x, degree = 3)}
#' @param n Number of rows in the data frame the data for estimation
#'@examples
#'# All these call with give the same result where
#'# 'data' is a hypothetical data frame. We can get a
#'# time-invariant estimate for x1 by:
#'\dontrun{
#' ddhazard(Surv(stop, event) ~ ddFixed(x1), data)
#'}
#'
#'# All of the calls below will yield the same result
#'# with a time-invariant intercept:
#'\dontrun{
#' ddhazard(Surv(stop, event)
#'     ~ ddFixed(1) + x1, data)
#' ddhazard(Surv(stop, event)
#'     ~ -1 + ddFixed(1) + x1, data)
#' ddhazard(Surv(stop, event)
#'     ~ ddFixed(rep(1, nrow(data))) + x1, data)
#' ddhazard(Surv(stop, event)
#'     ~ -1 + ddFixed(rep(1, nrow(data))) + x1, data)
#' ddhazard(Surv(stop, event)
#'     ~ ddFixed_intercept(nrow(data)) + x1, data)
#' ddhazard(Surv(stop, event)
#'     ~ -1 + ddFixed_intercept(nrow(data)) + x1, data)
#'}
#' @export
ddFixed <- function(object){
  object
}

#' @rdname ddFixed
#' @export
ddFixed_intercept <- function(n){
  out <- factor(rep(1, n), levels = 0:1)
  out
}

ddfixed_specials <- c("ddFixed_intercept", "ddFixed")
