# dynamichazard 0.3.0
The following has been changed or added:

* Refactored C++, R and test codes
* A bug have been fixed in the `get_risk_obj` when `is_for_discrete_model = TRUE` in the call. The issue was that individuals who were right censored in the middle of an interval was included despite that we do not know that they survive the entire interval. This will potentially affect the output for logit fits with  `ddhazard`
* The `ddhazard_boot` now provides the option of different learning rates to be used rather than one
* Error in computation of likelihood have been fixed
* Added the option to use relative change in the likelihood (not including the prior) as a convergence criteria instead of relative change in coefficients. The new option is selected by by calling `ddhazard` with `control = list(criteria = "delta_likeli", ...)`. The relative change in coefficient seems "preferable" as a default since it tends to not converge when the fit is has large "odd" deviation due to a few observations. The likelihood method though stops earlier for model does not have such deviation
* The default for kappa in the UKF have been changed to yield a weight on the first sigma point of 0.1 rather than 0
* Fixed memory leak in `ddhazard`
* New example have been added to the bootstrap vignette and other minor changes have been made
* Added information about the `residuals` method and a vignette "Diagnostics" with examples of usage of the `residuals` function
* Added a `hatvalues` method for `ddhazard`. These described "ddhazard" vignette and examples of usage are shown the vignette "Diagnostics"
* Fixed bug with default starting value with fixed effects
* Re-wrote the ddhazard vignette to have more a consistent notation
* Implemented new filter which use an series of rank-one approximation of the posterior in the correction step
* Added description of the posterior approximation method to the "ddhazard" vignette
* Added section about weights to the "ddhazard" vignette
* Removed the `rug` call in the shiny app demo, fixed a bug with the simulation function for the logit model and added the computation time of the estimation
* Inverses has been replaced by pseudo inverses. This will only have implications the matrices are rank deficient. The old option can be used by calling `ddhazard` with `control = list(use_pinv = FALSE, ...)`

# dynamichazard 0.2.0
The following have been added:

* Weights can be used in estimation. This is done by setting the `weights` argument when calling `ddhazard`
* Bootstrap of confidence bounds have been added. The function to bootstrap the estimates is `ddhazard_boot`. See the new vignette 'Bootstrap_illustration' for details
* S3 generic methods for `print` is added
* The method name for the variables have changed. What was previously denoted as truncated in now denoted as clipped
