# Find the design matrix and returns it left hand site and right hand site of
# regression equation
get_design_matrix = function(formula, data, response = T){
  Call <- match.call()
  indx <- match(c("formula", "data"),  names(Call), nomatch = 0)

  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")

  specials <- "ddFixed"

  # We may call this expression twice. Thus, we save it unevaluted
  get_form_exp <- quote(
    if(response)
      terms(formula, data = data, specials = specials) else
        eval(bquote(terms(update(formula, rep(1, .(nrow(data))) ~ ., data = data), data = data,
                          specials = specials))) # remove right hand site of formula
  )

  temp$formula <- eval(get_form_exp)

  # Check if we have a fixed intercept
  is_fixed <- attr(temp$formula, "specials")$ddFixed
  if(length(is_fixed) > 0){
    fixed_label <- dimnames(attr(temp$formula, "factors"))[[1]][is_fixed]
    is_fixed_intercept <- which(grepl("^ddFixed\\(\\ *1\\ *\\)$", fixed_label, perl = T))

    if(length(is_fixed_intercept) == 1){
      formula <- eval(parse(text = paste0(
        "update(formula, .~ ",
        " + ddFixed(rep(1,  nrow(",  deparse(substitute(data)),
        "))) + . -1 - ", fixed_label[is_fixed_intercept], ")")))

      # have to do this again
      temp$formula <- eval(get_form_exp)
    }}

  environment(temp$formula) <- parent.frame() # Needed if we use fixed intercept
  mf <- eval(temp, parent.frame())

  Y <- if(response) model.extract(mf, "response") else NULL

  if(!is.null(Y) && attr(Y, "type") == "right"){ # Change outcome if formula was Surv(stop, event) ~
    Y <- cbind(rep(0, nrow(Y)), Y)
    dimnames(Y)[[2]] <- c("start", "stop", "status")
    attr(Y, "type") <- "counting"
  }

  Terms <- terms(mf)
  fixed_terms_indicies <- attr(Terms, "specials")$ddFixed

  # From coxph
  if(length(fixed_terms_indicies) > 0){
    # First deal with fixed effects
    temppred <- attr(terms, "predvars")
    Terms1 <- Terms[fixed_terms_indicies - 1]
    fixed_terms <- model.matrix(Terms1, mf)
    fixed_terms <- fixed_terms[, colnames(fixed_terms) != "(Intercept)", drop = F] # remove intercept

    # Then deal with dynamic effects
    Terms2 <- Terms[-fixed_terms_indicies + 1]
    if (!is.null(temppred)) {
      attr(Terms2, "predvars") <- temppred[-(1 + fixed_terms_indicies)]
    }
    X <- model.matrix(Terms2, mf)
    renumber <- match(colnames(attr(Terms2, "factors")),
                      colnames(attr(Terms, "factors")))
    attr(X, "assign") <- c(0, renumber)[1 + attr(X, "assign")]
  } else{
    X <- model.matrix(Terms, mf)
    fixed_terms <- matrix(nrow = nrow(X), ncol = 0)

  }

  # Change fixed_terms name if ddFixed(1) or similar was used
  is_ddFixed_intercept <- grepl("^ddFixed\\(rep\\(\\s*1,\\s*nrow\\(.*\\)\\)\\)$", colnames(fixed_terms))
  if(any(is_ddFixed_intercept)){
    old_name <- colnames(fixed_terms)[is_ddFixed_intercept]
    new_name <- "ddFixed((Intercept))"

    attr(temp$formula, "term.labels")[
      attr(temp$formula, "term.labels") == old_name] <- new_name
    colnames(attr(temp$formula, "factors"))[
      colnames(attr(temp$formula, "factors")) == old_name] <- new_name
    rownames(attr(temp$formula, "factors"))[
      rownames(attr(temp$formula, "factors")) == old_name] <- new_name

    colnames(fixed_terms)[is_ddFixed_intercept] <- new_name
  }

  list(X = X, fixed_terms = fixed_terms, Y = Y, formula = temp$formula)
}

#' Function used in formula of \code{\link{ddhazard}} for time-invariant effects
#' @param object Expression that would be used in formula. E.g. \code{x} or \code{poly(x, degree = 3)}
#' @export
ddFixed <- function(object){
  object
}
