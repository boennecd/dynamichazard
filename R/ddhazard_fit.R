#' Function to perform EM algorithm for dynamic discrete hazard models. This
#' version perform no validation
#' @export
ddhazard_fit = function(X, tstart, tstop, is_event_in_bin,
                        a_0 = rep(0, n_parems), # ignorant prior
                        Q_0 = diag(.1, n_parems), # something large
                        F_ = diag(1, n_parems), # assume first order AR
                        Q = Q_0, risk_obj, n_parems,
                        n_max = 10^2, eps = 10^-3,
                        verbose = F, save_all_output = F,
                        order_ = 1,
                        est_Q_0 = T){
  if(order_ > 1)
    warning("No test made with order of ", order_)

  # Initalize
  d = risk_obj$d

  n_parems = length(a_0) / order_

  i = 0
  a_prev = a_0
  conv_criteria = function(prev, new_)
    norm(prev - new_, type = "2") / (norm(prev, type = "2") + 10^-10)
  conv_values = rep(NA_real_, n_max)

  if(save_all_output)
    all_output = list()

  # Run EM
  repeat{
    i = i + 1

    # E-step
    inner_output <-
      gen_kalman_filter_cpp(a_0 = a_0, Q = Q, Q_0 = Q_0, F_ = F_,
                            risk_sets = risk_obj$risk_sets,
                            I_len =  risk_obj$I_len,
                            d = risk_obj$d,
                            X = X, start = tstart,
                            stop = tstop, is_event_in_bin = risk_obj$is_event_in,
                            order_ = order_)

    # M-step
    with(inner_output, {
      a_0 <<- a_t_d_s[1, ]
      if(est_Q_0)
        Q_0 <<- V_t_d_s[, , 1]

      Q_ =  Q
      Q_[, ] = 0
      for(t in 2:(d + 1)){
        delta_t = risk_obj$I_len[t - 1]
        B = B_s[, , t - 1]

        V_less = V_t_d_s[, , t - 1]
        V = V_t_d_s[, , t]

        a_less = a_t_d_s[t - 1, ]
        a = a_t_d_s[t, ]

        Q_ = Q_ + ((a - F_ %*% a_less) %*% t(a - F_ %*% a_less) +
                     V -
                     F_ %*% B %*% V -
                     t(F_ %*% B %*% V) +
                     F_ %*% V_less %*% t(F_)) / delta_t
      }

      Q <<- Q_ / d
    })

    if(max((Q - t(Q))[1:n_parems, 1:n_parems]) >  sqrt(.Machine$double.eps))
      warning("Max diff in Q - Q^T is: ", max(Q - t(Q)))
    if(max(Q_0 - t(Q_0)) >  sqrt(.Machine$double.eps))
      warning("Max diff in Q_0 - Q^T_0 is: ", max(Q_0 - t(Q_0)))

    Q = (Q + t(Q)) / 2.0
    Q_0 = (Q_0 + t(Q_0)) / 2.0

    if(order_ > 1){ # CHANGED # TODO: I figure I should set the primaery element to zero, right?
      tmp_Q = Q[1:n_parems, 1:n_parems]
      Q[,] = 0
      Q[1:n_parems, 1:n_parems] = tmp_Q
    }

    conv_values[i] = conv_criteria(a_prev, a_0)
    if(save_all_output)
      all_output[[i]] = c(inner_output, list(conv_values = conv_values[i]))

    if(verbose && i %% 5 == 0)
      message("Finished iteration ", i, " with convergence criteria ", conv_values[i])
    if(conv_values[i] < eps){
      break
    } else{
      is_large_change = abs((a_prev - a_0) / (a_prev + 10^-2)) > 5
      if(any(is_large_change) && verbose){
        c_names = sapply(colnames(X), function(s_){ # TODO>this could done once and not in each iteration
          n_char = nchar(s_)
          if(n_char <= 15 + 3 + 3) s_ else paste0(substring(s_, 1, 15), "...", substring(s_, n_char - 3, n_char),
                                                  collapse = "")
        })

        message("These variable have big relative changes ([name], [new value], [previous value]):\n",
                paste0("(", apply(cbind(c_names, format(a_0, digits = 2), format(a_prev, digits = 2))[
                  is_large_change, , drop = F],
                  1, paste0, collapse = ", "), ")", collapse = "\t"))
      }
    }

    if(i == n_max){
      message("Function did not converge in ", n_max, " iterations")
      break
    }

    a_prev = a_0
  }

  if(save_all_output)
    all_output
  else
    c(inner_output, list(n_iter = i, conv_values = conv_values[seq_len(i)],
                         Q = Q, Q_0 = Q_0))
}
