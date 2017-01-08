[![Build Status](https://travis-ci.org/boennecd/dynamichazard.svg?branch=master,osx)](https://travis-ci.org/boennecd/dynamichazard)

dynamichazard
=============

The goal of dynamichazard is to estimate binary regression models with time-varying effects. The time-varying effects are estimated with state space models where the coefficients follow a given order random walk. The advantageous of using state space models is that you can extrapolate beyond the last observed time period

The estimation methods are implemented such:

-   They have linear time complexity in both time and the number of observations
-   Computation is done in `c++` with the `Armadillo` library. Consequently, they are fast and computation time can be reduced further with an optimized BLAS and LAPACK library
-   Use a `coxph` like formula from the survival package. Therefore, the methods are easily applied to panel datasets
-   Allows for time-varying covariates

For more details, see the ddhazard vignette at <https://cran.r-project.org/web/packages/dynamichazard/vignettes/ddhazard.pdf>

Installation
------------

You can install dynamichazard from github with:

``` r
install.packages("devtools")
devtools::install_github("dynamichazard/boennecd")
```

You can also download the package from CRAN by calling:

``` r
installed.packages("dynamichazard")
```

Example
-------

We will use the `aids` data set from the `JMbayes` package. The data set is from a a randomized clinical trial between two drugs for HIV or aids patients. The event is when the patient die. Patient are right censored at the end of the study. The data set is longitudinal/panel format with rows for patient. Though, the only longitudinal variable is the `CD4` count (T-cells count) which is presumably affected by the drug. Thus, we will not use it in the model. The other the other columns of interest are:

-   `AZT` is one of the two enrolment criteria. It indicates whether the patient was enrolled due to intolerance to the drug zidovudine or whether the drug failed prior to the study start
-   `prevOI` is the other enrolment criteria. Patients are enrolled either due AIDS diagnosis or two CD4 counts of 300 or fewer. The variable indicates which is the case
-   `drug` is either `ddC` or `ddI` depending on which of the two drugs the patient is randomly assigned to
-   `gender`

The analysis is given below with comments:

``` r
library(dynamichazard)
#> Loading required package: survival
library(survival)
library(JMbayes) # Contain the aids data set
#> Loading required package: MASS
#> Loading required package: nlme
#> Loading required package: splines

# We remove the data we dont neeed
aids <- aids[aids$Time == aids$stop, ]
aids <- aids[, !colnames(aids) %in% c("Time", "death", "obstime", "CD4")]

# A look at head of data
head(aids)
#>    patient drug gender prevOI         AZT start  stop event
#> 3        1  ddC   male   AIDS intolerance    12 16.97     0
#> 7        2  ddI   male noAIDS intolerance    18 19.00     0
#> 10       3  ddI female   AIDS intolerance     6 18.53     1
#> 14       4  ddC   male   AIDS     failure    12 12.70     0
#> 18       5  ddI   male   AIDS     failure    12 15.13     0
#> 19       6  ddC female   AIDS     failure     0  1.90     1
max(aids$stop)                  # Last observation time
#> [1] 21.4
max(aids$stop[aids$event == 1]) # Last person with event
#> [1] 19.07

# Fit model
fit <- ddhazard(
  Surv(stop, event) ~ AZT + gender + drug + prevOI,
  aids,
  model = "exp_trunc_time_w_jump", # The model we use from ddhazard. In short, 
                                   # this model assumes that event times are 
                                   # exponentially distributed
  by = .5,                         # Length of time intervals in state space 
                                   # model
  max_T = 19,                      # Last period we observe when modeling
  Q = diag(.1, 5),                 # Covariance matrix for state equation in 
                                   # first iteration
  Q_0 = diag(10000, 5),            # Covariance matrix for the prior
  control = list(eps = .1))
#> a_0 not supplied. One iteration IWLS of static glm model is used

# Plot the estimates. Dashed lines are 95% confidence bounds
plot(fit)
```

![](README-unnamed-chunk-2-1.png)

``` r

# Bootstrap the estimates
boot_out <- ddhazard_boot(fit, R = 10000) # R is number of bootstrap samples
#> Warning in ddhazard_boot(fit, R = 10000): Failed to estimate 322 times

# Plot bootstrap estimates. Dashed lines are 2.5% and 97.5% quantiles of the 
# bootstrap estimates. Transparent lines are bootstrap estimates
plot(fit, ddhazard_boot = boot_out)
#> Only plotting 500 of the boot sample estimates
```

![](README-unnamed-chunk-2-2.png)

Bootstrapping only slightly changes the confidence bounds. It seems that:

-   It is hard to tell the difference between the two drugs. The `ddi` may be more effective in the latter period (the estimates is negative) though the point-wise confidence bounds still contains 0. Further, this comment neglect that the confidence bounds are point-wise
-   Having aids rather than being enrolled (only) due to two CD4 counts of 300 or fewer increase your risk of dying
-   Males seems to be at lower risk in the first period

An example of a paper analyzing the CD4 count can be found in Guo & Carlin (2004). They also fit a static model (time-invariant coefficients) of the survival times with an exponential model. The estimates are comparable with those above as expected

References
==========

Guo, X., & Carlin, B. P. (2004). Separate and joint modeling of longitudinal and event time data using standard computer packages. *The American Statistician*, *58*(1), 16–24.
