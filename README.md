[![Build Status](https://travis-ci.org/boennecd/dynamichazard.svg?branch=master,osx)](https://travis-ci.org/boennecd/dynamichazard) \# dynamichazard

The goal of dynamichazard is to \[TO BE WRITTEN\]

Installation
------------

You can install dynamichazard from github with:

``` r
# install.packages("devtools")
devtools::install_github("dynamichazard/boennecd")
```

Example
-------

This is a basic example which shows you how to solve a common problem: \[TO BE WRITTEN\]

``` r
# See https://archive.ics.uci.edu/ml/machine-learning-databases/pima-indians-diabetes

# We download the data
diabetes <-
  read.table("https://archive.ics.uci.edu/ml/machine-learning-databases/pima-indians-diabetes/pima-indians-diabetes.data",
             sep = ",")
colnames(diabetes) <-
  c("times_pregnant", "glucose_cocen", "blood_pres", "skin_fold",
    "insulin", "BMI", "diabetes_pedigree", "age", "has_diabetes")

head(diabetes)
#>   times_pregnant glucose_cocen blood_pres skin_fold insulin  BMI
#> 1              6           148         72        35       0 33.6
#> 2              1            85         66        29       0 26.6
#> 3              8           183         64         0       0 23.3
#> 4              1            89         66        23      94 28.1
#> 5              0           137         40        35     168 43.1
#> 6              5           116         74         0       0 25.6
#>   diabetes_pedigree age has_diabetes
#> 1             0.627  50            1
#> 2             0.351  31            0
#> 3             0.672  32            1
#> 4             0.167  21            0
#> 5             2.288  33            1
#> 6             0.201  30            0
nrow(diabetes) # Number of samples before cleaning
#> [1] 768

# Remove observations with invalid values
diabetes <- diabetes[diabetes$BMI > 0, ]
diabetes <- diabetes[diabetes$blood_pres > 0, ]
diabetes <- diabetes[diabetes$skin_fold > 0, ]
diabetes$has_diabetes[diabetes$insulin >= 200]
#>  [1] 1 1 1 0 1 1 1 1 0 1 0 0 1 1 1 0 0 0 0 0 1 1 1 1 1 1 0 1 0 0 1 0 0 0 1
#> [36] 1 0 1 0 0 1 1 0 1 1 1 0 0 1 0 1 1 1 0 1 1 1 0 0 0 0 0 1 1 1 0 1 0 1 0
#> [71] 1 0 1 0 1 0 1 0 0 0 1 0 0 0 0 1 0 1 1

nrow(diabetes) # Number of samples after cleaning
#> [1] 537



# We load the two libraries we will need
library(survival); library(dynamichazard)
min(diabetes$age)                             # Youngest person in study
#> [1] 21
max(diabetes$age)                             # Oldest person in study
#> [1] 81
max(diabetes$age[diabetes$has_diabetes == 1]) # Last person who get diagnosed
#> [1] 70
                                              # with diabetes

# Notice that we have few cases and control in the end as illustrated below
diabetes$t <- diabetes$age - min(diabetes$age) + 1 # Time scale in study
xtabs(~ has_diabetes + t, diabetes)
#>             t
#> has_diabetes  1  2  3  4  5  6  7  8  9 10 11 12 13 14 15 16 17 18 19 20
#>            0 45 53 28 33 25 20 15 21  8 10  6  5  5  7  3  4  8  4  8  4
#>            1  4  7  4  7 12  6  5  7 10  4  8  5  7  2  4  7  4  6  1  5
#>             t
#> has_diabetes 21 22 23 24 25 26 27 28 29 30 31 32 33 34 35 36 37 38 39 40
#>            0  5  8  2  0  4  4  1  3  2  1  1  0  0  0  1  1  1  1  1  3
#>            1  7  3  9  1  5  6  2  1  2  3  5  4  3  2  1  2  1  3  1  1
#>             t
#> has_diabetes 41 42 43 45 50 61
#>            0  1  1  3  1  0  1
#>            1  0  1  0  0  1  0

# Thus, we set the last time (max_T) slighly lower later
# Fit the model
fit <- ddhazard(
  Surv(t, has_diabetes) ~
    ddFixed(1) +                        # These effects are time invariant
    ddFixed(diabetes_pedigree) +
    glucose_cocen + blood_pres +        # The others are time-varying
    skin_fold + BMI + times_pregnant
    ,
  data = diabetes,
  by = 1,                               # time interval lengths
  max_T = 42,                           # last period we observe
  Q_0 = diag(1e5, 5), Q = diag(1e-4, 5) # Covariances mats for state equation
  )
#> a_0 not supplied. One iteration IWLS of static glm model is used

# Plot the effects estimates with 95% point-wise confidence bounds
plot(fit)
```

![](README-unnamed-chunk-2-1.png)

``` r

# Print time invariant estimates
fit$fixed_effects
#>       (Intercept) diabetes_pedigree 
#>        -4.6035970         0.2231969

# Bootstrap the estimates
boot_out <- ddhazard_boot(fit, R = 1e2) # R is number of bootstrap samples
#> Warning in ddhazard_boot(fit, R = 100): Failed to estimate 27 times

# Plot bootstrapped estimates (transparent lines) and 2.5% and 97.5% quantiles
# (dashed-lines)
plot(fit, ddhazard_boot = boot_out)
```

![](README-unnamed-chunk-2-2.png)

``` r

# Make boxplot for the bootstrap estimate of the fixed effect
boxplot(t(tail(t(boot_out$t), 2)))
```

![](README-unnamed-chunk-2-3.png)
