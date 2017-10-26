# find the design matrix and returns the left hand side and right hand side of
# objects in the formula
get_design_matrix = function(
  formula, data, response = T, predictors = T,
  Terms = NULL, has_fixed_intercept = NULL, xlev = NULL){
  if(!response && !predictors)
    stop(sQuote("get_design_matrix"),
         " is called where neither respone or predictors is requested")

  if(!all(is.null(Terms) == c(is.null(has_fixed_intercept), is.null(xlev)))) # TODO: test error message
    stop(sQuote("Terms"), ", ", sQuote("has_fixed_intercept"), " and ",
         sQuote("xlev"), " must either all be NULL or not")

  is_first_call <- is.null(Terms)

  if(is_first_call){
    #####
    # See if there is a fixed intercept. If so, update the formula
    tt <- terms(formula, data = data, specials = "ddFixed_intercept")
    fixed_inter_terms <- attr(tt, "specials")[["ddFixed_intercept"]]
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
        attr(X, "assign") %in% (fixed_terms_indicies - attr(Terms, "response")))
    }

    fixed_terms <- X[, fixed_terms_indicies, drop = F]
    if(length(fixed_terms_indicies) > 0)
      X <- X[, -fixed_terms_indicies, drop = F]
    if(has_fixed_intercept){
      fixed_terms <- cbind(`(Intercept)` = rep(1, nrow(fixed_terms)),
                           fixed_terms)
      X <- X[, colnames(X) != "(Intercept)", drop = FALSE]
    }
  }

  # Change outcome if formula it was Surv(stop, event) ~
  if(!is.null(Y) && attr(Y, "type") == "right"){
    Y <- cbind(rep(0, nrow(Y)), Y)
    dimnames(Y)[[2]] <- c("start", "stop", "status")
    attr(Y, "type") <- "counting"
  }

  list(X = X, fixed_terms = fixed_terms, Y = Y, formula_used = formula_used,
       terms = Terms, has_fixed_intercept = has_fixed_intercept, xlev = xlev)
}

#' Auxiliary functions for fixed effects
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
#' ddhazard(Surv(stop, event) ~ ddFixed_intercept(what_ever) + x1, data)
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
#' @export
ddFixed_intercept <- function(object) NULL
