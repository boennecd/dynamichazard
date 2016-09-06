#' Find the design matrix and returns it LHS and RHS
#' @export
get_design_matrix = function(formula, data, response = T){
  Call <- match.call()
  indx <- match(c("formula", "data"),  names(Call), nomatch = 0)
  temp <- Call[c(1, indx)]
  temp[[1]] <- as.name("model.frame")

  temp$formula <- if(response)
    temp$formula = terms(formula, data = data) else
      temp$formula <- eval(bquote(terms(update(formula, rep(1, .(nrow(data))) ~ .), data = data))) # remove right hand site of formula

  mf <- eval(temp, parent.frame())

  Y <- if(response) model.extract(mf, "response") else NULL

  Terms <- terms(mf)
  X <- model.matrix(Terms, mf)

  list(X = X, Y = Y, formula = temp$formula)
}
