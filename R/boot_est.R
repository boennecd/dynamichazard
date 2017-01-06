#' @export
ddhazard_boot <- function(ddhazard_fit,  strata, unique_id, R = 100,
                          do_stratify_with_event = T, do_sample_weights = F,
                          print_errors = F){
  if(is.null(ddhazard_fit$risk_set) || is.null(ddhazard_fit$data))
    stop("Cannot bootstrap estimates when ddhazard has been called with control = list(save_risk_set = F, save_data = F, ...)")

  id <- ddhazard_fit$id
  X_Y <- get_design_matrix(formula = ddhazard_fit$formula, data = ddhazard_fit$data)

  # Find unique id and whether the individuals do die
  if(missing(unique_id))
    unique_id <- unique(id)

  if(missing(strata)){
    strata <- as.factor(rep(1, length(unique_id)))
  } else if(!is.factor(strata))
    strata <- as.factor(strata)

  if(length(strata) != length(unique_id)){
    stop("The strata argument should have same length as the unique ids")
  }
  if(do_stratify_with_event){
    is_event <- as.factor(tapply(ddhazard_fit$risk_set$is_event_in, id,
                                 function(x) sum(x > -1) == 1))
    strata <- interaction(strata, is_event, drop = T)
  }

  # Define function for boot
  statistic <- function(data, ran.gen){
    if(do_sample_weights){
      ws <- data
    } else{
      ws <- xtabs(~ran.gen)
    }

    # find which ids are included
    included_ids <- names(ws)
    if(is.numeric(id))
      included_ids <- as.numeric(included_ids)
    which_rows_are_included <- which(id %in% included_ids)

    # Adjust weights vector
    ws <- ws[match(id[which_rows_are_included], included_ids)]

    # Adjust X_Y
    boot_X_Y <- X_Y
    for(el_name in c("X", "fixed_terms", "Y")){
      boot_X_Y[[el_name]] <- boot_X_Y[[el_name]][which_rows_are_included, , drop = F]
    }
    # Transpose due to column-major ordering in c++
    boot_X_Y$X <- t(boot_X_Y$X)
    boot_X_Y$fixed_terms <- t(boot_X_Y$fixed_terms)

    # Adjust risk set
    boot_risk_set <- ddhazard_fit$risk_set
    boot_risk_set$is_event_in <-
      boot_risk_set$is_event_in[which_rows_are_included]

    index_map <- cbind(old_index = which_rows_are_included,
                       new_index = 1:length(which_rows_are_included))
    for(i in seq_along(boot_risk_set$risk_sets)){
      tmp <-
        boot_risk_set$risk_sets[[i]][boot_risk_set$risk_sets[[i]] %in%
                                       which_rows_are_included]

      boot_risk_set$risk_sets[[i]] <-
        index_map[match(tmp, index_map[,'old_index']), 'new_index']
    }

    # Estimate and return state vectors
    est <- NULL
    tryCatch(
      suppressWarnings(est <- ddhazard_no_validation(
        a_0 = ddhazard_fit$state_vecs[1,], Q_0 = ddhazard_fit$Q_0,
        F_ = ddhazard_fit$F_, verbose = F, Q = ddhazard_fit$Q,
        risk_set= boot_risk_set, X_Y = boot_X_Y,
        order = ddhazard_fit$order, model = ddhazard_fit$model,
        LR = ddhazard_fit$LR,
        n_fixed_terms_in_state_vec =
          ifelse(ddhazard_fit$control$fixed_terms_method == "E_step", ncol(X_Y$fixed_terms), 0),
        weights = ws,
        control = ddhazard_fit$control)$a_t_d_s),
      error = function(e){
        if(print_errors)
          print(e)
        est <<- rep(NA_real_, length(ddhazard_fit$state_vecs))
      })

    est
  }


  # Bootstrap estimate
  data <- rep(1, length(unique_id))
  names(data) <- unique_id
  boot_est <-
    boot::boot(data = data,
               statistic = statistic,
               R = R,
               sim = ifelse(do_sample_weights, "parametric", "ordinary"),
               ran.gen = function(data, ...){
                 f <- function(data){
                   tmp <- runif(length(data))
                   tmp <- tmp * (length(data)/ sum(tmp))
                   names(tmp) <- names(data)
                   tmp
                 }

                 out <- vector()
                 for(l in levels(strata)){
                   out <- c(out, f(data[strata == l]))
                 }

                 return(out)
               },
               strata = strata)

  n_fails <- sum(apply(boot_est$t, 1, function(x) any(is.na(x))))
  if(n_fails > 0)
    warning("Failed to estimate ", n_fails, " times")

  warning("Implement fixed effects in E-step")

  class(boot_est) <- c("ddhazard_boot", class(boot_est))
  boot_est$t_names <- c(sapply(colnames(ddhazard_fit$state_vecs), function(x) paste(x, 1:nrow(ddhazard_fit$state_vecs), sep = ":")))
  colnames(boot_est$t) <- boot_est$t_names
  colnames(boot_est$t0) <- colnames(ddhazard_fit$state_vecs)

  boot_est
}

# from boot.print
#' @export
print.ddhazard_boot <-
  function (x, digits = getOption("digits"), index = 1L:ncol(boot.out$t), ...)
{
  boot.out <- x
  sim <- boot.out$sim
  cl <- boot.out$call
  t <- matrix(boot.out$t[, index], nrow = nrow(boot.out$t))
  allNA <- apply(t, 2L, function(t) all(is.na(t)))
  ind1 <- index[allNA]
  index <- index[!allNA]
  t <- matrix(t[, !allNA], nrow = nrow(t))
  rn <- boot.out$t_names
  if (length(index) == 0L)
    op <- NULL
  else {
    t0 <- boot.out$t0[index]
    op <- cbind(
      t0,
      apply(t, 2L, mean, na.rm = TRUE) - t0,
      apply(t, 2L, mean, na.rm = TRUE, trim = .025) - t0,
      sqrt(apply(t, 2L, function(t.st) var(t.st[!is.na(t.st)]))))
    dimnames(op) <-
      list(rn, c("original", " bias  "," bias (truncated)"," std. error"))
  }

  cat("Bootstrap Statistics :\n")
  if (!is.null(op))
    print(op, digits = digits)
  invisible(boot.out)
}
