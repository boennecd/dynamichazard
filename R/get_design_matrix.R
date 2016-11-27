# Find the design matrix and returns it left hand site and right hand site of
# regression equation
get_design_matrix = function(formula, data, response = T){
  Call <- match.call()
  indx <- match(c("formula", "data"),  names(Call), nomatch = 0)
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")

  specials <- "ddFixed"
  temp$formula <- if(response)
    temp$formula = terms(formula, data = data, specials = specials) else
      temp$formula <- eval(bquote(terms(update(formula, rep(1, .(nrow(data))) ~ .), data = data,
                                  specials = specials))) # remove right hand site of formula

  mf <- eval(temp, parent.frame())

  Y <- if(response) model.extract(mf, "response") else NULL

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

  list(X = X, fixed_terms = fixed_terms, Y = Y, formula = temp$formula)
}

#' Function used in formula of \code{\link{ddhazard}} for time-invariant effects
#' @param object Expression that would be used in formula. E.g. \code{x} or \code{poly(x, degree = 3)}
#' @export
ddFixed <- function(object){
  object
}
