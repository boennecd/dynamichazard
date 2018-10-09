get_design_matrix <- function(
  formula, data, response = T, predictors = T,
  Terms = NULL, has_fixed_intercept = NULL, xlev = NULL, fixed = NULL,
  random = NULL){
  if(is.null(fixed) && is.null(random))
    return(.get_design_matrix_one_frm (
      formula = formula, data = data, response = response,
      predictors = predictors, Terms = Terms,
      has_fixed_intercept = has_fixed_intercept, xlev = xlev))

  .get_design_matrix_two_frms(
    data = data, response = response, predictors = predictors, Terms = Terms,
    xlev = xlev, fixed = fixed, random = random)
}

# find the design matrix and returns the left hand side and right hand side of
# objects in the formula
.get_design_matrix_one_frm <- function(
  formula, data, response = T, predictors = T,
  Terms = NULL, has_fixed_intercept = NULL, xlev = NULL){
  if(!response && !predictors)
    stop(sQuote("get_design_matrix"),
         " is called where neither respone or predictors is requested")

  if(!all(
    is.null(Terms) == c(is.null(has_fixed_intercept), is.null(xlev))))
    stop(sQuote("Terms"), ", ", sQuote("has_fixed_intercept"), " and ",
         sQuote("xlev"), " must either all be NULL or not")

  if(is_first_call <- is.null(Terms)){
    #####
    # See if there is a fixed intercept. If so, update the formula
    tt <- terms(formula, data = data, specials = "ddFixed_intercept")
    fixed_inter_terms <- attr(tt, "specials")[["ddFixed_intercept"]]
    if(length(fixed_inter_terms) > 1)
      stop("more than one ", sQuote("ddFixed_intercept"), " term")
    has_fixed_intercept <- length(fixed_inter_terms) > 0
    formula_used <- if(has_fixed_intercept){
      attr(tt, "intercept") <- 1
      tt <- drop.terms(tt, fixed_inter_terms - attr(tt, "response"),
                       keep.response = TRUE)
      attributes(tt) <- attributes(formula)
      tt
    } else
      formula
  }

  do_keep_random_intercept <- .keep_random_intercept(formula, data = data)

  #####
  # Get model.frame, model.matrix and outcome. The two latter if needed
  if(is_first_call){
    Call <- match.call()
    indx <- match(c("formula", "data"),  names(Call), nomatch = 0)

    temp <- Call[c(1, indx)]
    temp[[1]] <- as.name("model.frame")

    formula_used <- terms(formula_used, data = data, specials = "ddFixed")
    temp[[2]] <- formula_used

    if(!response)
      temp$formula <- delete.response(temp$formula)
    if(!predictors)
      temp$formula <- update(temp$formula, . ~ -1)

    mf <- eval(temp, parent.frame())

  } else {
    if(!response)
      Terms <- delete.response(Terms)
    if(!predictors)
      Terms <- update(Terms, . ~ -1)

    mf <- model.frame(Terms, data, xlev = xlev)

    formula_used <- Terms
  }

  Y <- if(response) model.response(mf) else NULL

  if(!predictors){
    Terms <- NULL
    X <- NULL
    fixed_terms <- NULL
    xlev <- NULL

  } else {
    Terms <- terms(mf)
    X <- model.matrix(Terms, mf)
    xlev <- .getXlevels(Terms, mf)

    fixed_terms_indicies <- attr(Terms, "specials")[["ddFixed"]]
    if(is.null(fixed_terms_indicies)){
      fixed_terms_indicies <- c()
    } else {
      fixed_terms_indicies <- which(
        attr(X, "assign") %in%
          (fixed_terms_indicies - attr(Terms, "response")))
    }

    fixed_terms <- X[, fixed_terms_indicies, drop = F]
    if(length(fixed_terms_indicies) > 0)
      X <- X[, -fixed_terms_indicies, drop = F]
    if(has_fixed_intercept){
      fixed_terms <- cbind(`(Intercept)` = 1, fixed_terms)
      if(!do_keep_random_intercept)
        X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  }

  list(X = X, fixed_terms = fixed_terms, Y = .set_Y_to_counting(Y),
       formula_used = formula_used,
       terms = Terms, has_fixed_intercept = has_fixed_intercept, xlev = xlev)
}

.keep_random_intercept <- function(formula, data){
  tt <- terms(formula, specials = "ddFixed_intercept", data = data)
  fixed_inter_terms <- attr(tt, "specials")[["ddFixed_intercept"]]
  if(length(fixed_inter_terms) == 0)
    return(FALSE)

  if(length(fixed_inter_terms) > 1)
    stop("more than one ", sQuote("ddFixed_intercept"), " term")

  fixed_expr <- attr(tt, "variables")[[fixed_inter_terms + 1L]]
  attr(eval(fixed_expr), "random_intercept")
}

.get_design_matrix_two_frms <- function(
  data, response = T, predictors = T, Terms = NULL, xlev = NULL, fixed = NULL,
  random = NULL){
  if(!is.null(Terms) || !is.null(xlev))
    stop("method not implemented")
  stopifnot(inherits(fixed, "formula"), inherits(random, "formula"))

  . <- function(formula, data, ...){
    tt <- terms(formula, data, specials = c("ddFixed", "ddFixed_intercept"))
    if(any(sapply(attr(tt, "specials"), length) > 0))
      stop("Do not use ", sQuote("ddFixed"), " or ",
           sQuote("ddFixed_intercept"),
           " when you use ", sQuote("random"), " and ", sQuote("fixed"))

    cl <- match.call()
    cl[[1L]] <- quote(model.frame)
    eval(cl, parent.frame())
  }

  if(!predictors)
    fixed <- update(fixed, . ~ 1)

  mf <- .(fixed, data)
  Y <- if(response)
    model.response(mf) else NULL

  if(predictors){
    Terms_fixed <- terms(mf)
    xlev_fixed <- .getXlevels(Terms_fixed, mf)
    fixed_terms <- model.matrix(Terms_fixed, mf)

    mf <- model.frame(random, data)
    Terms_random <- terms(mf)
    xlev_random <- .getXlevels(Terms_random, mf)
    X <- model.matrix(Terms_random, mf)
  } else
    X <- fixed_terms <- Terms_fixed <- Terms_random <- xlev_fixed <-
    xlev_random <- NULL

  list(X = X, fixed_terms = fixed_terms, Y = .set_Y_to_counting(Y),
       terms = list(fixed = Terms_fixed, random = Terms_random),
       xlev  = list(fixed = xlev_fixed , random = xlev_random))
}

.set_Y_to_counting <- function(Y){
  # Change outcome if formula it was Surv(stop, event) ~
  if(!is.null(Y) && attr(Y, "type") == "right"){
    Y <- cbind(0, Y)
    dimnames(Y)[[2]] <- c("start", "stop", "status")
    attr(Y, "type") <- "counting"
  }
  Y
}

#' Auxiliary Functions for Fixed Effects
#' @description
#' Functions used in formula of \code{\link{ddhazard}} for time-invariant effects. \code{ddFixed_intercept} is only used for the intercept.
#'
#' @param object expression that would be used in formula. E.g. \code{x} or \code{poly(x, degree = 3)}.
#'
#'@examples
#'# we can get a time-invariant effect of `x1` by
#'\dontrun{
#' ddhazard(Surv(stop, event) ~ ddFixed(x1), data)
#'}
#'
#'# all of the calls below will yield the same result with a time-invariant
#'# intercept:
#'\dontrun{
#' ddhazard(Surv(stop, event) ~ ddFixed_intercept() + x1, data)
#' ddhazard(Surv(stop, event) ~ -1 + ddFixed_intercept() + x1, data)
#'}
#' @export
ddFixed <- function(object){
  if(all(object == 1)){
    # check that the function is not called as part of predict where we may
    # just have ones

    funcs <- lapply(sys.calls(), "[[", 1)
    if(!any(sapply(funcs, function(x)
      any(sapply(x, function(z) is.name(z) && z == as.name("predict"))))))
      stop("All elements in call to ", sQuote("ddFixed"), " is one.",
           " This is potentially a depreciated way of fixing the intercept.",
           " See ", sQuote("?ddFixed"))
  }

  object
}

#' @rdname ddFixed
#'
#' @param random_intercept \code{TRUE} if a zero mean time-varying process
#' should be included at as an additional term. Only relevant in stationary
#' models. See the \code{type} argument in \code{\link{PF_EM}}.
#'
#' @export
ddFixed_intercept <- function(random_intercept = FALSE){
  if(!is.logical(random_intercept) || length(random_intercept) != 1L)
    stop(sQuote("random_intercept"), " needs to be a scalar logical")

  structure(list(), random_intercept = random_intercept)
}
