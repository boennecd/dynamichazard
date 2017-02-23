# dynamichazard 0.3.0
The following has been changed or added:

* Inverses has been replaced by pseudo inverses. This will only have implications the matrices are rank deficient. The old option can be used by calling `ddhazard` with `control = list(use_pinv = FALSE, ...)`
* A bug have been fixed in the `get_risk_obj` when `is_for_discrete_model = TRUE` in the call. The issue was that individuals who were right censored in the middle of an interval was included despite that we do not know that they survive the entire interval. This will potentially affect the output for logit fits with  `ddhazard`
* The `ddhazard_boot` now provides the option of different learning rates to be used rather than one
* Error in computation of likelihood have been fixed
* Added the option to use relative change in the likelihood (not including the prior) as a convergence criteria instead of relative change in coefficients. The new option is selected by by calling `ddhazard` with `control = list(criteria = "delta_likeli", ...)`. The relative change in coefficient seems "preferable" as a default since as it tends to not converge when the fit is has large "odd" deviation due to a few observations. The likelihood method though stops earlier for model does not have such deviation
* The default for kappa in the UKF have been changed to yield a weight on the first sigma point of 0.1 rather than 0
* Fixed memory leak in `ddhazard`
* New example have been added to the bootstrap vignette and other minor changes have been made

# dynamichazard 0.2.0
The following have been added:

* Weights can be used in estimation. This is done by setting the `weights` argument when calling `ddhazard`
* Bootstrap of confidence bounds have been added. The function to bootstrap the estimates is `ddhazard_boot`. See the new vignette 'Bootstrap_illustration' for details
* S3 generic methods for `print` is added
* The method name for the variables have changed. What was previously denoted as truncated in now denoted as clipped
