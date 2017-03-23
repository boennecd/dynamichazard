#' @title Bootstrap for \code{\link{ddhazard}}
#' @param ddhazard_fit Returned object from a \code{\link{ddhazard}} call
#' @param strata Strata to sample within. These need to be on an individual by individual basis and not rows in the design matrix
#' @param unique_id Unique ids where entries match entries of \code{strata}
#' @param R Number of bootstrap estimates
#' @param do_stratify_with_event \code{TRUE} if sampling should be by strata of whether the individual has an event. An interaction factor will be made if \code{strata} is provided
#' @param do_sample_weights \code{TRUE} if weights should be sample instead of individuals
#' @param LRs Learning rates in decreasing order which will be used to estimate the model.
#' @param print_errors \code{TRUE} if errors should be printed when estimations fails
#'
#' @description
#' See the vignette 'Bootstrap illustration'. The \code{do_stratify_with_event} may be useful when either cases or non-cases are very rare to ensure that the model estimation succeds.
#'
#' @return
#' An object like returned from the \code{\link[boot]{boot}} function
#'
#' @seealso
#' \code{\link{ddhazard}}, \code{\link[=plot.fahrmeier_94]{plot}}
#'
#' @export
ddhazard_boot <- function(ddhazard_fit,  strata, unique_id, R = 100,
                          do_stratify_with_event = F, do_sample_weights = F,
                          LRs = ddhazard_fit$control$LR * 2^(0:(-4)), print_errors = F){
  if(is.unsorted(-LRs, strictly = T) && any(LRs <=0))
    stop("LRs need to be stricly decreasing postive numbers")

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

  # Change starting value for fixed effects
  ddhazard_fit$control$fixed_parems_start <- ddhazard_fit$fixed_effects

  # Define indicies for fixed effects in the state space vector
  n_fixed <- ncol(X_Y$fixed_terms)
  n_varying <- ncol(ddhazard_fit$state_vecs)/ddhazard_fit$order
  n_periods <- nrow(ddhazard_fit$state_vecs)
  n_out <- n_fixed +                           # number of parameters returned
    n_varying * n_periods * ddhazard_fit$order # by statistic function

  # We need to make adjustments if fixed effects are estimated in the E-step
  a_0 <- ddhazard_fit$state_vecs[1,]
  if(ddhazard_fit$control$fixed_terms_method == "E_step" &&
     n_fixed > 0){
    indicies_fixed_coef_in_state <- 1:n_fixed +  n_varying

    # We need to add indicies to Q
    Q_new <- ddhazard_fit$Q_0 # This has the right dimension
    Q_new[,] <- 0
    Q_new[1:n_varying, 1:n_varying] <- ddhazard_fit$Q
    ddhazard_fit$Q <- Q_new

    # We need the right F matrix
    ddhazard_fit$F_ <-  get_F(
      ddhazard_fit$order, n_varying, n_fixed, T)

    # and to a_0
    if(ddhazard_fit$order == 1){
      a_0 <- c(a_0, ddhazard_fit$control$fixed_parems_start)
    } else if(ddhazard_fit$order == 2){
      first_half <- 1:(length(a_0)/2)
      a_0 <-
        c(a_0[first_half], ddhazard_fit$control$fixed_parems_start, a_0[-first_half])
    } else
      stop("Method not implemented for order equal to ", ddhazard_fit$order)
  } else
    indicies_fixed_coef_in_state <- vector()

  # Define function for boot
  statistic <- function(data, ran.gen){
    if(do_sample_weights){
      ws <- data
    } else{
      ws <- xtabs(~ ran.gen)
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

    # Adjust risk set
    boot_risk_set <- ddhazard_fit$risk_set
    boot_risk_set$is_event_in <-
      boot_risk_set$is_event_in[which_rows_are_included]

    # Permutate if needed
    if(ddhazard_fit$control$permu)
      eval(get_permu_data_exp(boot_X_Y[1:3], boot_risk_set, ws))

    # Add fixed terms to design matrix if fixed terms are estimated in the
    # E-step
    if(ddhazard_fit$control$fixed_terms_method == "E_step" &&
       n_fixed > 0){
      boot_X_Y$X <- cbind(boot_X_Y$X, boot_X_Y$fixed_terms)
      boot_X_Y$fixed_terms <- matrix(nrow = nrow(boot_X_Y$X), ncol = 0)
    }

    # Transpose due to column-major ordering in c++
    boot_X_Y$X <- t(boot_X_Y$X)
    boot_X_Y$fixed_terms <- t(boot_X_Y$fixed_terms)

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
    did_succed <- FALSE
    for(l in LRs){
      tryCatch({
        suppressWarnings(est <- ddhazard_no_validation(
          a_0 = a_0, Q_0 = ddhazard_fit$Q_0,
          F_ = ddhazard_fit$F_, verbose = F, Q = ddhazard_fit$Q,
          risk_set= boot_risk_set, X_Y = boot_X_Y,
          order = ddhazard_fit$order, model = ddhazard_fit$model,
          LR = l,
          n_fixed_terms_in_state_vec =
            ifelse(ddhazard_fit$control$fixed_terms_method == "E_step", n_fixed, 0),
          weights = ws,
          control = ddhazard_fit$control))

        out <- c(est$a_t_d_s[, !seq_len(ncol(est$a_t_d_s)) %in% indicies_fixed_coef_in_state])
        if(n_fixed > 0){
          if(ddhazard_fit$control$fixed_terms_method == "E_step")
            out <- c(out, est$a_t_d_s[1, indicies_fixed_coef_in_state]) else
              out <- c(out, est$fixed_effects)
        }
        did_succed <- TRUE
        est <- out
        }, error = function(e){
          if(l > min(LRs))
            return(invisible())

          if(print_errors)
            print(e)
          est <<- rep(NA_real_, n_out)
        })

      if(did_succed)
        break
    }

    est
  }


  # Bootstrap estimates
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
  if(n_fails == R){
    stop("Failed to estimate with all samples")
  } else if(n_fails > 0)
    warning("Failed to estimate ", n_fails, " times")

  class(boot_est) <- c("ddhazard_boot", class(boot_est))
  boot_est$t_names <-
    c(c(sapply(colnames(ddhazard_fit$state_vecs),
               function(x) paste(x, 1:nrow(ddhazard_fit$state_vecs) - 1, sep = ":t"))),
      names(ddhazard_fit$fixed_effects))

  colnames(boot_est$t) <- boot_est$t_names
  names(boot_est$t0) <- boot_est$t_names
  boot_est$sim <- "ordinary"

  boot_est
}

get_frac_n_weights <- function(R, a){
  #####
  # Function to make linear interpolation with weights from normal density

  k <- floor((R + 1) * a)
  if(k == 0 || k == R)
    stop("Sample of ", R, " is too small to give ", a, " confidence bounds")

  w <- (qnorm(a) - qnorm(k/(R + 1))) / (qnorm((k + 1) / (R + 1)) - qnorm(k/(R + 1)))
  w_k <- 1 - w
  w_k_p_1 <- w

  list(k = k, w_k = w_k, w_k_p_1 = w_k_p_1)
}

