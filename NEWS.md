# dynamichazard 0.6.7
* solve issue because of changes in `all.equal`.
* fix other minor issues with test on CRAN.

# dynamichazard 0.6.6
* two errors in `PF_get_score_n_hess` are fixed. One is that a off diagonal 
  block in the observed information matrix was not computed. The other is that
  parts of the score and observed information matrix was only correct if 
  parts of them were multiplied by the duplication matrix.
* `nlopt` is no longer used in mode optimization. A Newton–Raphson method 
  is used instead. This seems a bit faster in some cases and does not fail in 
  some cases where `nlopt` did.
* A `fix_seed` argument is added to `PF_control`. `fix_seed = FALSE` combined
  with averaging and a low number of particles seems to yield better results.
* fixed a bug in `PF_EM` when some periods do not have any observations.

# dynamichazard 0.6.5
* Minor bug fixes.

# dynamichazard 0.6.4
* Particle filtering implementation is changed. One may get slightly different 
  results.
* Estimation with a moderate amount of fixed parameters is much faster now with 
  the particle filters.
* Averaging is possible with `PF_EM`.
* `PF_get_score_n_hess` is added to compute the approximate negative 
  observation matrix and score vector.
* `nu` in `PF_control` scales the scale matrix to get an identical covariance 
  matrix.
* Fixed bug in `get_Q_0`.
* Fixed bug in M-step in EM algorithm with particle smoothers.

# dynamichazard 0.6.2
* The cloglog link function is added.
* `predict.ddhazard` has been re-written. The output with `type = "term"` 
has changed. It now yields a list of lists. Each list contains a list for each
`new_data` row. The time zero index is no longer included if `tstart` and 
`tstop` is not matched in `new_data`. Parallel computing is no longer supported. 
It likely did not yield any reduction in computation time with the previous
implementation. Calls with `type = "term"` now uses the `tstart` and `tstop` 
argument and supports predictions in the future. A covariance matrix is added 
to the terms in the predictions. 
* A method has been added to plot the survival curve for a given observation
which may have time-varying covariates.

# dynamichazard 0.6.1
* Fixed a bug in `type == "VAR"` in particle filters in the smoothing proposal distribution. 
This has a major impact for most calls. 
* Fixed a bug in `type == "VAR"` in particle filters where the transition from the 
time zero state to time one was not used in the M-step estimation. This only has 
a larger impact for short series. 
* Fixed a bug in `method == "bootstrap_filter"` where a wrong covariance matrix was 
used for the proposal distribution. 
* Fixed a bug in `method == "AUX_normal_approx_w_particles"` where a wrong covariance matrix was 
used for the proposal distribution. 
* Fixed a bug in the `logLik.PF_clouds`. The log-likelihood approximation was 
too high especially for the auxiliary particle filters. 
* An option is added to use a (multivariate) t-distribution as the proposal 
distribution in particle filters. 
* A few miscellaneous functions have been added for particle filter methods.

# dynamichazard 0.6.0
* One can fit first order Vector vector autoregression models with the particle
filter. 
* The averages used in the Taylor approximation in the backward smoothing have 
changed so the results differ for `PF_EM`.
* A bug is fixed with the `..._w_particles` methods so results have changed.
* Mode estimation is now done in the proposal distribution with more than one
iteration.
* A `random` and `fixed` argument is added to `PF_EM` as an alternative way 
to specify the random and fixed effect parts. 

# dynamichazard 0.5.2
* Fixed effects is estimated faster with `PF_EM` and can be estimated 
with `model = "exponential"`.
* The covariance matrix in `ddhazard` objects are no longer degenerate (e.g., 
in the case where a second order random walk is used). Instead the dimension 
is equal to the dimension of the error term.
* The seed argument in `PF_EM` has been moved from the `control` list. Further, 
there is a `PF_control` which should preferably be used to construct the object 
for the `control` argument of `PF`.
* A more stable and parallel estimation method have been added to `static_glm`.
* `PF_EM` uses $Q_0$ instead of $Q$ for the artificial prior and a bug have 
been fixed for sampling in the initial state in the backward filter. This have 
changed the output.
* Fixed bug in `PF_EM` with seed argument. The new way to get reproducible is 
to call `f1 <- PF_EM(...); .GlobalEnv$.Random.seed" <- f1$seed; f2 <- eval(f1$call)`
kinda as in `simulate.lm`.

# dynamichazard 0.5.0
* Fixed issues with close to equal non-integer stop and start times.
* Re-factored code to use more generic setup for distribution families.
* Changed EKF and UKF to use Poisson counts instead of the three different types used before for exponential distributed arrival times. This means that the `model` argument to `ddhazard` should be changed from `"exp_bin"`, `"exp_clip_time"` or `"exp_clip_time_w_jump"` to `"exponential"`.

# dynamichazard 0.4.0
* Parallel version of `glm` is used to find the first state vector.
* Default values of `ddhazard` control argument is changed.
* Particle filter and smoothers alternative are available.

# dynamichazard 0.3.0
The following has been changed or added:

* Refactored C++, R and test codes.
* A bug have been fixed in the `get_risk_obj` when `is_for_discrete_model = TRUE` in the call. The issue was that individuals who were right censored in the middle of an interval were included despite that we do not know that they survive the entire interval. This will potentially affect the output for logit fits with  `ddhazard`.
* The `ddhazard_boot` now provides the option of different learning rates to be used rather than one if the first fit succeeds.
* Error in computation of likelihood have been fixed.
* Added the option to use relative change in the likelihood (not including the prior) as a convergence criteria instead of relative change in coefficients. The new option is selected by by calling `ddhazard` with `control = list(criteria = "delta_likeli", ...)`. The relative change in coefficient seems "preferable" as a default since it tends to not converge when the fit is has large "odd" deviation due to a few observations. The likelihood method though stops earlier for model does not have such deviation.
* The default for kappa in the UKF have been changed to yield a weight on the first sigma point of 0.1 rather than 0.
* Fixed memory leak in `ddhazard`.
* New example have been added to the bootstrap vignette and other minor changes have been made.
* Added a `hatvalues` method for `ddhazard`. These described "ddhazard" vignette and examples of usage are shown the vignette "Diagnostics".
* Added information about the `residuals` method and a vignette "Diagnostics" with examples of usage of the `residuals` function.
* Fixed bug with default starting value with fixed effects.
* Re-wrote the ddhazard vignette to have more a consistent notation.
* Implemented new filter which use an series of rank-one approximation of the posterior in the correction step.
* Global mode approximation have been added.
* Added description of the posterior approximation methods to the "ddhazard" vignette.
* Added section about weights to the "ddhazard" vignette.
* Removed the `rug` call in the shiny app demo, fixed a bug with the simulation function for the logit model and added the computation time of the estimation to the output.
* Inverses has been replaced by pseudo inverses. This will only have implications the matrices are rank deficient. The old option can be used by calling `ddhazard` with `control = list(use_pinv = FALSE, ...)`.
* Refactored the UKF code. This version no longer supports the `exp_combined` method.

# dynamichazard 0.2.0
The following have been added:

* Weights can be used in estimation. This is done by setting the `weights` argument when calling `ddhazard`.
* Bootstrap of confidence bounds have been added. The function to bootstrap the estimates is `ddhazard_boot`. See the new vignette 'Bootstrap_illustration' for details.
* S3 generic methods for `print` is added.
* The method name for the variables have changed. What was previously denoted as truncated in now denoted as clipped.
