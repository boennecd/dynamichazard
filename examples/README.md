# Examples

This directory contains the following examples

### First order vector auto-regression
with the particle filter and smoothers. 
See [firstOrderVAR.Rmd](firstOrderVAR.Rmd) for the source code 
and [firstOrderVAR.html](https://htmlpreview.github.io/?https://github.com/boennecd/dynamichazard/blob/master/examples/firstOrderVAR.html)
for the returned `.html` file. 

### First order vector auto-regression with restricted parameters
with the particle filter and smoothers. Here we estimate models with fewer
parameters by restricting the model.
See [restrictedVAR.Rmd](restrictedVAR.Rmd) for the source code 
and [restrictedVAR.html](https://htmlpreview.github.io/?https://github.com/boennecd/dynamichazard/blob/master/examples/restrictedVAR.html)
for the returned `.html` file. 

### First order random walk model
with the particle filter and smoothers and the extended Kalman filter. 
See [RW.Rmd](RW.Rmd) for the source code 
and [RW.html](https://htmlpreview.github.io/?https://github.com/boennecd/dynamichazard/blob/master/examples/RW.html)
for the returned `.html` file.

### IID model
with the particle filter and smoothers with a comparison with `lme4::glmer`.
See [iid.Rmd](iid.Rmd) for the source code 
and [iid.html](https://htmlpreview.github.io/?https://github.com/boennecd/dynamichazard/blob/master/examples/iid.html)
for the returned `.html` file.

### Overhead with Multithreading
The file shows how much is gained by adding more threads on one hardware, 
operating system, and compiler. The results will likely differ on another 
setup. See [overhead.Rmd](overhead.Rmd) for the source code 
and [overhead.html](https://htmlpreview.github.io/?https://github.com/boennecd/dynamichazard/blob/master/examples/overhead.html)
for the returned `.html` file.
