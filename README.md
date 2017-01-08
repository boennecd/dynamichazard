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

-   `prevOI` is a factor for whether previous opportunistic infection of AIDS is present at study entry
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

# A look at head of data
head(aids, 25)
#>    patient  Time death       CD4 obstime drug gender prevOI         AZT
#> 1        1 16.97     0 10.677078       0  ddC   male   AIDS intolerance
#> 2        1 16.97     0  8.426150       6  ddC   male   AIDS intolerance
#> 3        1 16.97     0  9.433981      12  ddC   male   AIDS intolerance
#> 4        2 19.00     0  6.324555       0  ddI   male noAIDS intolerance
#> 5        2 19.00     0  8.124038       6  ddI   male noAIDS intolerance
#> 6        2 19.00     0  4.582576      12  ddI   male noAIDS intolerance
#> 7        2 19.00     0  5.000000      18  ddI   male noAIDS intolerance
#> 8        3 18.53     1  3.464102       0  ddI female   AIDS intolerance
#> 9        3 18.53     1  3.605551       2  ddI female   AIDS intolerance
#> 10       3 18.53     1  6.164414       6  ddI female   AIDS intolerance
#> 11       4 12.70     0  3.872983       0  ddC   male   AIDS     failure
#> 12       4 12.70     0  4.582576       2  ddC   male   AIDS     failure
#> 13       4 12.70     0  2.645751       6  ddC   male   AIDS     failure
#> 14       4 12.70     0  1.732051      12  ddC   male   AIDS     failure
#> 15       5 15.13     0  7.280110       0  ddI   male   AIDS     failure
#> 16       5 15.13     0  8.602325       2  ddI   male   AIDS     failure
#> 17       5 15.13     0  8.602325       6  ddI   male   AIDS     failure
#> 18       5 15.13     0  6.708204      12  ddI   male   AIDS     failure
#> 19       6  1.90     1  4.582576       0  ddC female   AIDS     failure
#> 20       7 14.33     0  6.782330       0  ddC   male   AIDS intolerance
#> 21       7 14.33     0  5.385165       2  ddC   male   AIDS intolerance
#> 22       7 14.33     0  4.472136       6  ddC   male   AIDS intolerance
#> 23       7 14.33     0  3.162278      12  ddC   male   AIDS intolerance
#> 24       8  9.57     1  3.464102       0  ddI female noAIDS intolerance
#> 25       8  9.57     1  1.000000       2  ddI female noAIDS intolerance
#>    start  stop event
#> 1      0  6.00     0
#> 2      6 12.00     0
#> 3     12 16.97     0
#> 4      0  6.00     0
#> 5      6 12.00     0
#> 6     12 18.00     0
#> 7     18 19.00     0
#> 8      0  2.00     0
#> 9      2  6.00     0
#> 10     6 18.53     1
#> 11     0  2.00     0
#> 12     2  6.00     0
#> 13     6 12.00     0
#> 14    12 12.70     0
#> 15     0  2.00     0
#> 16     2  6.00     0
#> 17     6 12.00     0
#> 18    12 15.13     0
#> 19     0  1.90     1
#> 20     0  2.00     0
#> 21     2  6.00     0
#> 22     6 12.00     0
#> 23    12 14.33     0
#> 24     0  2.00     0
#> 25     2  6.00     0

# We remove the data we dont neeed
aids <- aids[aids$Time == aids$stop, ]
aids <- aids[, !colnames(aids) %in% c("Time", "death", "obstime", "CD4")]

# A look at head of data
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
boot_out <- ddhazard_boot(fit, R = 1000) # R is number of bootstrap samples
#> Warning in ddhazard_boot(fit, R = 1000): Failed to estimate 30 times

# Plot bootstrap estimates. Dashed lines are 2.5% and 97.5% quantiles of the 
# bootstrap estimates. Transparent lines are bootstrap estimates
plot(fit, ddhazard_boot = boot_out)
#> Only plotting 500 of the boot sample estimates
```

![](README-unnamed-chunk-2-2.png)

Bootstrapping only slightly changes the confidence bounds. It seems that: \* It is hard to tell the difference between the two drugs. The `ddi` may be more effective in the latter period (the estimates is negative) though the point-wise confidence bounds still contains 0. Further, this comment neglect that the confidence bounds are point-wise \* Having aids rather than two CD4 counts of 300 or fewer increase your risk of dying \* Males seems to be at lower risk in the first period
