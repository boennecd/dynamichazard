## ----setup_knitr, echo=FALSE, cache=FALSE--------------------------------
knitr::render_sweave()

#####
# Hook to set par
with(new.env(), {
  par_default <- function(cex_mult = 1, ...){
    cex <- .75 * cex_mult

    list(
      mar = c(5, 5, 2, 2),
      bty = "L",
      xaxs = "i",
      pch=16,
      cex= cex,
      cex.axis = 1.25,
      cex.lab = 1.4,
      lwd= 1)
  }

  knitr::knit_hooks$set(
    par_1x1 =
      function(before, options, envir) {
        if(!options$par_1x1)
          return()

        if (before){
          par(mfcol = c(1, 1))
          par(par_default(.8))
        }
      },

    par_3x3 =
      function(before, options, envir) {
        if(!options$par_3x3)
          return()

        if (before){
          par(mfcol = c(3, 3))
          tmp <- par_default(.8)
          tmp$mar <- tmp$mar + c(0.5, 0, 0, 0)
          par(tmp)
        }
    })
})

######
# Chunk options
knitr::opts_chunk$set(
  echo = TRUE, warning = F, message = F, dpi = 144,
  cache = T,
  fig.align = "center",
  abbreviate_after = 60,
  stop = FALSE,

  # See opts_hooks definition
  fig.height = -1, fig.width = -1)

#####
# Alter width and height of figures
knitr::opts_hooks$set(
  fig.height = function(options) {
    if(options$fig.height > 0)
      return(options)

    if (!is.null(options$par_3x3) && options$par_3x3){
      options$fig.height <- 8
    } else
      options$fig.height <- 3

    options
  },

  fig.width = function(options) {
    if(options$fig.width > 0)
      return(options)

    if (!is.null(options$par_3x3) && options$par_3x3) {
      options$fig.width <- 9
    } else
      options$fig.width <- 5

    options
  })

## ----setup_other, echo=FALSE, cache=FALSE--------------------------------
#####
# R options
options(digits = 3, scipen=7, width = 60)

#####
# Remove comments from xtable
print.xtable <- with(new.env(), {
  org_fun <- print.xtable
  function(..., comment = FALSE)
    org_fun(..., comment = comment)
})

#####
# Define colors that are good for red-green colorblinds
cols <-cbind(
  r=c(91,0,23,255,8,255,4,0,0,255,255),
  g=c(0,255,169,232,0,208,255,0,79,21,0),
  b=c(12,233,255,0,91,198,4,255,0,205,0))
cols <- apply(cols, 1, function(x)
  rgb(x[1], x[2], x[3], maxColorValue = 255))

# barplot(rep(1, length(cols)), col = cols)
palette(c("black", cols, "darkgray"))

#####
# Load and attach libaries
library(splines)
library(stringr)
library(splines)
library(plotrix)
library(grDevices)
library(mvtnorm)
library(tcltk)
library(xtable)
library(zoo, quietly = T, warn.conflicts = FALSE)

## ----load_hd_dat, echo = FALSE-------------------------------------------
#####
# Load hd data from sub folder
hd_dat <- readRDS("HDS/HDs.RDS")

# Few have data from time zero so we set a few days in as time zero
new_start <- 24 * 4
hd_dat$tstart <- pmax(new_start, hd_dat$tstart)
hd_dat$tstart <- hd_dat$tstart - new_start
hd_dat$tstop <- hd_dat$tstop - new_start

# We need to remove the records that ends before or at the starting time
# sum(hd_dat$tstart >= hd_dat$tstop) # Number of rows thrown away
hd_dat <- hd_dat[hd_dat$tstart < hd_dat$tstop, ]
hd_dat$serial_number <- droplevels(hd_dat$serial_number)

# Re-scale time to months
tmp <- 24 * 30
hd_dat$tstart  <- hd_dat$tstart / tmp
hd_dat$tstop <- hd_dat$tstop / tmp

# Make sure that data is sorted
hd_dat <- hd_dat[order(hd_dat$serial_number, hd_dat$tstart), ]

#####
# Fill in blanks with carry the last observation forward
# Define function to fill in the blanks
library(zoo, quietly = T)
func <- function(x)
  na.locf0(c(0, x))[-1]
func <- compiler::cmpfun(func)

# Use the function
for(n in colnames(hd_dat)["smart_12" == colnames(hd_dat)]){
  hd_dat[[n]] <- unlist(
    tapply(hd_dat[[n]], as.integer(hd_dat$serial_number), func),
    use.names = F)
}

## ----remove_versions_w_few_and_winsorize, echo = FALSE-------------------
#####
# Remove version with few unique disks
n_per_model <-
  xtabs(~ model, hd_dat, subset = !duplicated(serial_number))

# We take those larger than a given size
factor_cut <- 400
models_to_keep <- names(n_per_model)[n_per_model >= factor_cut]
hd_dat <- hd_dat[hd_dat$model %in% models_to_keep, ]
hd_dat$model <- droplevels(hd_dat$model)

#####
# Winsorize
win_lvl <- .99
hd_dat$smart_12 <- pmin(hd_dat$smart_12, quantile(
  hd_dat$smart_12, win_lvl))

## ----define_get_pretty_model_factors, echo = FALSE-----------------------
#####
# Define function to make factor levels shorter
library(stringr)
get_pretty_model_factors <- function(x){
  f <- function(lvls){
    lvls <- str_replace(lvls, "^model", "")
    lvls <- str_replace(lvls, "^[A-z]+\\ ", "")

    lvls
  }

  if(class(x) == "fahrmeier_94"){
    colnames(x$state_vecs) <- paste("Param.", f(colnames(x$state_vecs)))
    return(x)
  }

  f(x)
}

## ----model_stats, echo = FALSE, results='asis', cache=FALSE--------------
#####
# Make data frame to find stats
library(dynamichazard)
tmp_dat <- get_survival_case_weights_and_data(
  Surv(tstart, tstop, fails) ~ model,
  data = hd_dat, by = 1, max_T = 60, use_weights = F,
  id = hd_dat$serial_number)

#####
# Make time cut variable and find # disk and # failures
tmp_dat$X$t_cut <- cut(tmp_dat$X$t, breaks = seq(0, 60, 20),
                       right = FALSE)
stats <- by(tmp_dat$X, list(tmp_dat$X$model, tmp_dat$X$t_cut), function(x){
  c("#D" = length(unique(x$serial_number)),
    "#F" = sum(tapply(x$Y, x$serial_number, any, default = 0)))
})

.names <- dimnames(stats)
stats <- sapply(stats, function(x) if(is.null(x)) c(0, 0) else x)

#####
# Format final tabel
n_models <- length(.names[[1]])
rnames <- .names[[1]]
tbl_dat <- lapply(1:length(.names[[2]]), function(i){
  x_vals <- stats[, 1:n_models + n_models * (i - 1L)]
  cnames <- paste(rownames(x_vals), .names[[2]][i])
  structure(t(x_vals), dimnames = list(rnames, cnames))
})
tbl_dat <- do.call(cbind, tbl_dat)
tbl_dat <- tbl_dat[order(tbl_dat[, 1], decreasing = TRUE), ]

# Convert to string
tbl_dat <- structure(
  apply(tbl_dat, 2, sprintf, fmt = "%d"),
  dimnames = dimnames(tbl_dat))

# Format 2 row header
# See https://stackoverflow.com/a/32490565/5861244
tmp_lvls <- paste0("$t\\in", levels(tmp_dat$X$t_cut), "$")
headers_lvl1 <- paste(
  "\\multicolumn{2}{c}{", tmp_lvls, "}", collapse = " & ")
headers_lvl2 <- paste(substr(colnames(tbl_dat), 0, 2), collapse = " & ")

#####
# Find first observation and first failure and add to table
tbl_xtra <- by(hd_dat, hd_dat$model, function(x) min(x$tstop[x$fails == 1]), simplify = FALSE)
tbl_xtra <- do.call(rbind, tbl_xtra)
tbl_xtra[] <- sprintf(tbl_xtra, fmt = "%.3f")

rmatch <- match(row.names(tbl_dat), row.names(tbl_xtra))

tbl_dat <- cbind(tbl_xtra[rmatch], tbl_dat)

# Set header of final tabel
headers_lvl1 <- paste(
  paste(rep(" & ", ncol(tbl_xtra)), # Not -1 as we have the row.names
        collapse = ""),
  "&", headers_lvl1)
headers_lvl2 <- paste(
  "HD version", paste(
    colnames(tbl_xtra), collapse = " & "),
  headers_lvl2, sep = " & ")

headers_lvl2 <- gsub("#", "\\\\#", headers_lvl2)

# Make xtable output
rownames(tbl_dat) <- get_pretty_model_factors(rownames(tbl_dat))

xtbl <- xtable(tbl_dat,
       caption = "Summary statistics for each hard disk versions. The hard disk version is indicated by the first column. The number of disks is abbreviated as '\\#D' and total failures is abbreviated as '\\#F'. The $t\\in[x, y)$ indicates which time interval that the figures applies to.",
       label = "tab:modelDat")

align(xtbl)[-1] <- "r"
align(xtbl)[1] <- "l"
print_out <- capture.output(print(
  xtbl,
  include.colnames = FALSE,
  sanitize.text.function = force,
  hline.after = 0,
  add.to.row = list(
    pos = list(-1),
    command = paste(headers_lvl1, " \\\\[0.33em] \n", headers_lvl2, " \\\\\ \n"))))

print_out <- sapply(
  print_out, gsub,
  pattern = "\\begin{table}[ht]",
  replacement = "\\begin{sidewaystable}",
  fixed = TRUE,
  USE.NAMES = FALSE)

print_out <- sapply(
  print_out, gsub,
  pattern = "\\end{table}",
  replacement = "\\end{sidewaystable}",
  fixed = TRUE,
  USE.NAMES = FALSE)

print_out <- paste0(print_out, collapse = "\n")
cat(print_out)

# Cleanup
rm(tmp_dat, stats, tbl_dat, tbl_xtra)

#####
# Note: You can check the figures against https://www.backblaze.com/blog/hard-drive-failure-rates-q2-2016/

## ----first_hd_fit--------------------------------------------------------
library(splines)
library(dynamichazard)
# Assign model formula
frm <- Surv(tstart, tstop, fails) ~
      -1 +              # I remove the intercept to not get a reference level
                        # for the hard disk model factor
      model +           # model is the hard disk version
      ns(smart_12,                   # Use a natural cubic spline for power
         knots = c(3, (1:5)*10),     # cycle count
         Boundary.knots = c(0, 60))

# Fit model
system.time(           # Used to get the computation time
  ddfit <- ddhazard(
    formula = frm,
    hd_dat,
    id = hd_dat$serial_number,
    Q_0 = diag(1, 24), # Covariance matrix for first state
    Q = diag(.01, 24), # Covariance matrix for transition
    by = 1,            # Length of intervals
    max_T = 60,        # Last time we observe when estimating
    control = list(
      method = "EKF",
      eps = .001)))

## ----ddfit_change_fac_lvls, echo = FALSE---------------------------------
ddfit <- get_pretty_model_factors(ddfit)

## ----ST3TB,echo=FALSE,par_1x1 = TRUE, fig.cap="Parameter for the hard disk version ST3000DM001. It is shown in a single plot as it differs from the parameters for the other factor levels shown in figure~\\ref{fig:otherfaclvl1} by being a lot larger around month 30.",cache=FALSE----
plot(ddfit, cov_index = 7)

## ----otherfaclvl1, echo = FALSE, par_3x3 = TRUE, fig.cap="Parameters for factor levels for the hard disk version with EKF with a single iteration in the correction step.",cache=FALSE----
plot(ddfit, cov_index = c(1:6, 8:10))

## ----smart_12_plot, echo = FALSE, par_1x1 = TRUE, fig.cap="Plot of predicted terms on the linear predictor scale for different values of number of power cycles.",cache=FALSE----
tmp <- data.frame(
  model = rep(hd_dat$model[1], 3),
  smart_12 = c(1, 10, 40))

preds <- predict(ddfit, tmp, type = "term", sds = T)

is_ns <- grepl(
  "^ns\\(smart_12", dimnames(preds$terms)[[3]])

# Find lower and upper bounds
preds$lbs <- preds$terms[,, is_ns]  - 1.96 * preds$sds[,, is_ns]
preds$ubs <- preds$terms[,, is_ns] + 1.96 * preds$sds[,, is_ns]

# Plot
cols <- rev(1:3)
xs <- ddfit$times
plot(range(xs), range(preds$lbs, preds$ubs), type = "n",
     xlab = "Time", ylab = "Power cycle term", col = cols)
abline(h = 0, lty = 2)
matplot(xs, preds$terms[,, is_ns], type = "l", add = T, lty = 1,
        col = cols)

# Add confidence bounds
for(i in 1:dim(preds$terms)[2]){
  icol <- adjustcolor(cols[i], alpha.f = 0.1)
  polygon(c(xs, rev(xs)), c(preds$ubs[,i], rev(preds$lbs[,i])),
          col = icol, border = NA)
  lines(xs, preds$ubs[,i], lty = 2, col = cols[i])
  lines(xs, preds$lbs[,i], lty = 2, col = cols[i])
}

# Add legend
legend(
  "bottomright", bty = "n",
  legend = paste0(tmp$smart_12, " Cycles"),
  lty = rep(1, nrow(tmp)),
  col = cols,
  cex = par()$cex * 2.5)

## ----setup_for_smart_12_illu, echo=FALSE---------------------------------
min_obs <- 50
quant <- .1
smart_12_illu_cap <- paste0("Plot showing the mean and ", quant * 100, "\\% quantile of the SMART 12 attribute for each version with those hard disk at risk at the start of each interval. The dashed lines are the quantile curves. Values for a given version in a given interval is excluded if there are less than ", min_obs, " hard disk at risk at the start of the interval. The circle radiuses are proportional to the fraction of hard disk that fail in the month. The transparency of the cirles are inversely log proportional to the number of hard disks at risk.")

## ----smart_12_illu,echo = FALSE, fig.cap=paste(smart_12_illu_cap), par_1x1=TRUE,cache=FALSE----
#####
# Compute figures for the plot
library(plotrix)
library(grDevices)
tmp_dat <- get_survival_case_weights_and_data(
  Surv(tstart, tstop, fails) ~ smart_12 + model,
  data = hd_dat, by = 1, max_T = 60, use_weights = F,
  id = hd_dat$serial_number)

smart_12_mean <- by(
  tmp_dat$X, tmp_dat$X$model, function(X)
    tapply(X$smart_12, X$t, mean))

smart_12_quant <- by(
  tmp_dat$X, tmp_dat$X$model, function(X)
    tapply(X$smart_12, X$t, quantile, probs = .1))

n_obs <- by(
  tmp_dat$X, tmp_dat$X$model, function(X)
    tapply(X$smart_12, X$t, length))
max_n_obs <- max(unlist(n_obs))

n_fails <- by(
  tmp_dat$X, tmp_dat$X$model, function(X)
    tapply(X$Y, X$t, sum))

fail_ratios <- mapply("/", n_fails, n_obs)

#####
# Make plot
plot(c(1, 60), c(0, max(unlist(smart_12_mean))), type = "n",
     xlab = "Month", ylab = "Number of power cycles")

for(i in seq_along(smart_12_mean)){
  col <- if(names(smart_12_mean)[i] == "ST3000DM001") "Red" else "Black"

  # We remove the points with less than x observations
  is_in <- n_obs[[i]] > min_obs
  y_mean <- smart_12_mean[[i]][is_in]
  y_quant <- smart_12_quant[[i]][is_in]
  x <- as.integer(names(smart_12_mean[[i]]))[is_in]

  fail_ratio <- fail_ratios[[i]]
  radius <- fail_ratio[is_in] * 10 + 1e-3

  lines(x, y_mean, col = col)
  lines(x, y_quant, col = col, lty = 2)

  n_obs_use <- n_obs[[i]][is_in]
  for(j in seq_along(x)){
    col_circle <- adjustcolor(
      col, alpha.f = 1 - .75 * (log(max_n_obs) - log(n_obs_use[j])) /
        (log(max_n_obs) - log(min_obs)))
    draw.circle(x[j], y_mean[j], radius[j],
                col = col_circle, border = col_circle)
  }
}

#####
# Cleanup
rm(tmp_dat)

## ----refit_ddfit_on_subset, echo = FALSE---------------------------------
#####
# Take subset
hd_dat_sub <- hd_dat
rm(hd_dat)
hd_dat_sub <- hd_dat_sub[hd_dat_sub$model != "ST3000DM001", ]
hd_dat_sub$model <- droplevels(hd_dat_sub$model)

#####
# Re-fit model
new_call <- ddfit$call
new_call$data <- as.name(quote(hd_dat_sub))
new_call$id <- as.call(quote(hd_dat_sub$serial_number))
new_call$Q_0 <- as.call(quote(diag(1, 23)))
new_call$Q <- as.call(quote(diag(.01, 23)))

ddfit <- eval(new_call)
ddfit <- get_pretty_model_factors(ddfit)

## ----subset_smart_12_plot,echo=FALSE,ref.label='smart_12_plot', par_1x1 = TRUE,fig.cap="Similar plot to figure~\\ref{fig:smart_12_plot} for the model without the ST3000DM001 hard disk version.",cache=FALSE----

## ----subset_EKF_xtra-----------------------------------------------------
system.time(
  ddfit_xtr <- ddhazard(
    formula = frm,
    data = hd_dat_sub,     # Data set without ST3000DM001
    id = hd_dat_sub$serial_number,
    Q_0 = diag(1, 24 - 1), # -1 due to removal of a factor level
    Q = diag(.01, 24 - 1),
    by = 1,
    max_T = 60,
    control = list(
      method = "EKF",
      eps = .001,
      NR_eps = .001))) # Tolerance for extra iterations in correction step

## ----define_add_hist, echo = FALSE---------------------------------------
add_hist <- with(new.env(), {
  tmp_dat <- get_survival_case_weights_and_data(
    Surv(tstart, tstop, fails) ~ smart_12,
    data = hd_dat_sub, by = 1, max_T = 60, use_weights = F,
    id = hd_dat_sub$serial_number)

  # Find number of observations through time for each model
  n_obs <- by(
    tmp_dat$X, tmp_dat$X$model, function(X)
      tapply(X$smart_12, X$t, length))
  n_obs <- xtabs(~ model + t, tmp_dat$X)

  xright <- as.numeric(colnames(n_obs))
  xleft <- xright - 1

  # Cleanup
  rm(tmp_dat)

  # Define function to plot histogram in background
  function(i){
    y_lim <- par("usr")[3:4]
    obs_dat <- n_obs[i, ]

    rect(
      xleft = xleft, xright = xright,
      ybottom = rep(y_lim[1], length(xright)),
      ytop <- (.25 * obs_dat) / max(obs_dat) * diff(y_lim) + y_lim[1],
      col = rgb(0, 0, 0, .2),
      border = rgb(1, 1, 1, 0.25),
      lwd = par()$lwd * 2)
  }
})

## ----def_fig.cap_bar_expla, echo=FALSE-----------------------------------
fig.cap_bar_expla <- " Grey transparent bars indicates the number of individuals at risk for the specific model version. Heights are only comparable within the model versions."

## ----subset_EKF_xtra_vs_wo, par_3x3 = TRUE, fig.cap=paste0("Comparison of results with and without extra iterations in the correction step with the EKF. The red curve is the estimate with extra iterations.", fig.cap_bar_expla),cache=FALSE----
for(i in 1:9){
  # Shown in figure \ref&fig:subset_EKF_xtra_vs_wo;
  plot(ddfit, cov_index = i)
  plot(ddfit_xtr, cov_index = i, add = TRUE, col = "red")
  add_hist(i) # Function defined in this paper
}

## ----show_n_iter---------------------------------------------------------
c(one_it = ddfit$n_iter, more_its = ddfit_xtr$n_iter)

## ----EKF_cleanup, echo = FALSE-------------------------------------------
rm(ddfit_xtr, hd_dat)

## ----sigma_pts,echo=FALSE, par_1x1=TRUE, fig.cap = "Illustration of sigma points in the example from equation~\\eqref{eqn:UKFEx}. The dashed lines are the contours of the density given by $\\emNotee{\\vec{a}}{t}{t - 1}$ and $\\emNotee{\\mat{V}}{t}{t - 1}$. The full lines are the direction given by the columns of the Cholesky decomposition. The filled circles are sigma points with $(\\alpha,\\beta,\\kappa) = (1,0,1)$ and the open circles are the sigma points with $(\\alpha,\\beta,\\kappa) = (1/3,0,1)$. The point at $(0,0)$ is a sigma point for both sets for hyperparameters.",cache=FALSE----
#####
# Draw points
library(mvtnorm)
set.seed(7912351)
x.points <- seq(-3, 3,length.out=100)
y.points <- x.points
z <- matrix(0,nrow=100,ncol=100)
mu <- c(0,0)
sigma <- matrix(c(2,1,1,1),nrow=2)
for (i in 1:100) {for (j in 1:100) {
  z[i,j] <- dmvnorm(c(x.points[i],y.points[j]),
                    mean=mu,sigma=sigma)
}}

#####
# Plot contours
plot(c(-3, 3), c(-3, 3), xlab = "", ylab = "", type = "n",
     xlim = c(-3, 3), ylim = c(-3, 3))
contour(x.points, y.points, z, nlevels = 10,
        drawlabels = FALSE, axes = FALSE,
        frame.plot = FALSE, add = TRUE,
        lty = 3)

#####
# Compute Cholesky decomposition
decomp <- chol(sigma)
abline(h = 0)
abline(a = 0, b = decomp[1, 2] / decomp[2, 2])

#####
# Add two sets of sigma points
q <- 2
l1 <- 1
l2 <- -1

pts <- rbind(
  c(0, 0),
  sqrt(q + l1) * t(decomp),
  - sqrt(q + l1) * t(decomp),
  sqrt(q + l2) * t(decomp),
  - sqrt(q + l2) * t(decomp))

points(pts[, 1], pts[, 2],
       pch = c(rep(16, 5), rep(1, 4)), cex = par()$cex * 3)

## ----define_ekf_low_Q_0, echo = FALSE------------------------------------
ukf_ex_Q_0 <- 0.1

## ----ukf_large_Q_0_comp--------------------------------------------------
ddfit_ukf <- ddhazard(
  formula = frm,
  data = hd_dat_sub,
  id = hd_dat_sub$serial_number,
  Q_0 = diag(10, 24 - 1), # Larger value
  Q = diag(.01, 24 - 1),
  by = 1,
  max_T = 60,
  control = list(
    method = "UKF",       # Use the UKF
    eps = .01))           # Decreased to get a fit

## ----ukf_large_Q_0_change_fac_names, echo = FALSE------------------------
ddfit_ukf <- get_pretty_model_factors(ddfit_ukf)

## ----ukf_large_Q_0, par_3x3=TRUE, fig.cap = "Predicted parameters with the UKF used on the hard disk failure dataset where $\\mat{Q}_0$ has large entries in the diagonal.",cache=FALSE----
# Shown in figure \ref&fig:ukf_large_Q_0;
plot(ddfit_ukf, cov_index = 1:9)

## ----ukf_small_Q_0_est, echo = FALSE-------------------------------------
#####
# Re-fit with lower Q_0 entries
new_call <- ddfit_ukf$call
new_call$Q_0 <- diag(ukf_ex_Q_0, 23)
ctrl <- eval(new_call$control)
ctrl$eps <- .001
new_call$control <- ctrl
ddfit_ukf <- eval(new_call)

## ----ukf_small_Q_0, par_3x3=TRUE, echo = FALSE,fig.cap = paste0("Similar plot to figure~\\ref{fig:ukf_large_Q_0} but where the diagonal entries of $\\mat{Q}_0$ are $", ukf_ex_Q_0, "$. The black curve is the estimates from the the EKF with one iteration in the correction step.", fig.cap_bar_expla),cache=FALSE----
for(i in 1:9){
  plot(ddfit, cov_index = i)
  plot(ddfit_ukf, cov_index = i, add = TRUE, col = "red")
  add_hist(i)
}

## ----show_EKF_can_have_large_comp----------------------------------------
new_call <- ddfit$call
new_call$Q_0 <- diag(10000, 23)
ddfit_large_Q_0 <- eval(new_call)

## ----show_EKF_can_have_large_change_fac_names, echo = FALSE--------------
ddfit_large_Q_0 <- get_pretty_model_factors(ddfit_large_Q_0)

## ----show_EKF_can_have_large, par_3x3=TRUE, fig.cap = paste0("Predicted parameters with the EKF with one iteration in the correction used on the hard disk failure dataset with $\\mat{Q}_0$ has low entries (", ddfit$Q_0[1,1], ") and large entries (", ddfit_large_Q_0$Q_0[1,1], ") in the diagonal. The red curves are the estimates with large entries in the diagonal.", fig.cap_bar_expla),cache=FALSE----
for(i in 1:9){
  # Shown in figure \ref&fig:show_EKF_can_have_large;
  plot(ddfit_large_Q_0, cov_index = i, col = "red")
  plot(ddfit, cov_index = i, add = TRUE)
  add_hist(i)
}

## ----ukf_cleanup, echo = FALSE-------------------------------------------
rm(ddfit_large_Q_0, ddfit_ukf)

## ----echo=FALSE, cache = FALSE-------------------------------------------
set.seed(914587)

## ----fit_SMA_hd_fail-----------------------------------------------------
ddfit_SMA <- ddhazard(
  formula = frm,
  data = hd_dat_sub,
  by = 1,
  max_T = 60,
  id = hd_dat_sub$serial_number,
  Q_0 = diag(1, 23),
  Q = diag(0.01, 23),
  control = list(
    eps = 0.001,
    method = "SMA",                  # Use SMA
    posterior_version = "cholesky")) # The Cholesky method in algorithm \ref&alg:approxModeChol;

## ----fit_GMA_hd_fail-----------------------------------------------------
ddfit_GMA <- ddhazard(
  formula = frm,
  data = hd_dat_sub,
  by = 1,
  max_T = 60,
  id = hd_dat_sub$serial_number,
  Q_0 = diag(1, 23),
  Q = diag(0.01, 23),
  control = list(
    eps = 0.001,
    method = "GMA"))                 # Use GMA instead

## ----EKF_vs_SMA_n_GMA, par_3x3 = TRUE, fig.cap= paste0("Predicted parameters using the EKF, the GMA and the SMA for the hard disk failure data set. The gray lines are the parameters from the SMA, red lines are parameters from the GMA and the black lines are the parameters form the EKF.", fig.cap_bar_expla), echo = FALSE, cache = FALSE----
for(i in 1:9){
  plot(ddfit, cov_index = i)
  plot(ddfit_SMA, cov_index = i, col = "gray40", add = T)
  plot(ddfit_GMA, cov_index = i, col = "red", add = T)
  add_hist(i)
}

## ----SMA_seed_large_Q0, echo=FALSE, cache = FALSE------------------------
set.seed(914587)

## ----SMA_w_large_Q_0-----------------------------------------------------
new_call <- ddfit_SMA$call
new_call$Q_0 <- diag(100000000, 23)
ddfit_SMA_large_Q_0 <- eval(new_call)

## ----SMA_w_large_Q_0_plot, par_3x3=TRUE, fig.cap = paste0("Similar plot to figure~\\ref{fig:show_EKF_can_have_large} but with the SMA instead. The red curves are the estimates with the SMA with large entries in the diagonal of $\\mat{Q}_0$. The black curves are the same fit using the EKF as in figure~\\ref{fig:show_EKF_can_have_large}.", fig.cap_bar_expla), echo = FALSE----
ddfit_SMA_large_Q_0 <- get_pretty_model_factors(ddfit_SMA_large_Q_0)
for(i in 1:9){
  # Shown in figure \ref&fig:SMA_w_large_Q_0_plot;
  plot(ddfit_SMA_large_Q_0, cov_index = i, col = "red")
  plot(ddfit, cov_index = i, add = TRUE)
  add_hist(i)
}

## ----SMA_n_GMA_cleanup, echo = FALSE-------------------------------------
rm(ddfit_SMA_large_Q_0, ddfit_SMA, ddfit_GMA)

## ----order_2_est---------------------------------------------------------
# Define new formula
frm_fixed <-
  Surv(tstart, tstop, fails) ~ -1 + model +
  ddFixed(ns(smart_12, knots = c(3, (1:2) * 10), Boundary.knots = c(0, 30)))

ddfit_fixed_E <- ddhazard(
  formula = frm_fixed,
  data = hd_dat_sub,
  by = 1,
  order = 2,
  max_T = 60,
  id = hd_dat_sub$serial_number,
  Q_0 = diag(10, 32),
  Q = diag(0.01, 16),
  control = list(
    method = "GMA",
    NR_eps = .001,
    eps = 0.001,
    LR = .2,
    fixed_terms_method = "E_step")) # Use E-step method

## ----order_2_plot, echo = FALSE, par_3x3 = TRUE, fig.cap=paste0("Predicted parameters for some of the factor levels with the second order random walk. The red lines are parameters with fixed effects estimated in the E-step and the black lines is the first order random walk where all parameters are time-varying.", fig.cap_bar_expla),cache=FALSE----
ddfit_fixed_E <- get_pretty_model_factors(ddfit_fixed_E)

for(i in 1:9){
  # Shown in figure \ref&fig:order_2_plot;
  plot(ddfit_fixed_E, cov_index = i, col = "red")
  plot(ddfit, cov_index = i, add = TRUE)
  add_hist(i)
}

## ----order_2_plot_fixed, echo = FALSE, par_1x1 = TRUE, fig.cap="Fixed effects estimates for the SMART 12 attribute on the linear predictor scale using the E-step method.",cache=FALSE----
expr <- tail(colnames(attr(terms(ddfit_fixed_E$formula), "factors")), 1)
x <- with(
  list(smart_12 = (x_org <- seq(0, 70, length.out = 1000))),
  eval(parse(text = expr)))

spline_E <- drop(x %*% ddfit_fixed_E$fixed_effects)

plot(x_org, spline_E,
     type = "l", xlab = "SMART 12 attribute", ylab = "Linear predictor term")

## ----smart_12_through_time, echo = FALSE, par_1x1 = TRUE, fig.cap="Quantiles of the SMART 12 attribute through time. The quantiles are at 5\\%, 10\\%, \\ldots, 95\\%. This is the data set without the ST3000DM001 version.",cache=FALSE----
#####
# Find quantiles through time
tmp_dat <- get_survival_case_weights_and_data(
  Surv(tstart, tstop, fails) ~ smart_12,
  data = hd_dat_sub, by = 1, max_T = 60, use_weights = F,
  id = hd_dat_sub$serial_number)
quants <- tapply(
  tmp_dat$X$smart_12, tmp_dat$X$t, quantile, probs = (1:19)/20)
quants <- do.call(cbind, quants)

#####
# Plot
x <- as.numeric(colnames(quants))
plot(range(x), range(quants), type = "n", xlab = "Month", ylab = "SMART 12 attribute")
col <- rgb(0, 0, 0, alpha = .1)
for(i in 1:floor(nrow(quants) / 2)){
  lb <- quants[i, ]
  ub <- quants[nrow(quants) - (i - 1), ]
  polygon(c(x, rev(x)), c(ub, rev(lb)), col = col, border = NA)
}

#####
# Cleanup
rm(tmp_dat)

## ----hd_dat_ex_setup, echo = FALSE---------------------------------------
#####
# Alter digits for the next print
old_digits <- getOption("digits")
options(digits = 5)

## ----hd_dat_ex-----------------------------------------------------------
hd_dat_sub[
  1:10, c("serial_number", "model", "tstart", "tstop", "smart_12")]

## ----hd_dat_ex_cleanup, echo = FALSE-------------------------------------
options(digits = old_digits)

## ----binning_fig, echo=FALSE, results="hide", fig.cap = "Illustration of a data set with 7 individuals with time-varying covariates. Each horizontal line represents an individual. Each number indicates a start time and/or stop time in the initial data. A cross indicates that new covariates are observed while a filled circle indicates that the individual has an event. An open circle indicates that the individual is right censored. Vertical dashed lines are time interval borders. The symbols for the covariate vectors and stop times are shown for observation a.", fig.height=3.5, fig.width=6, par_1x1 = TRUE,cache=FALSE----
# Alter par values
par(mar = c(1, 5, 1, 2), cex = par()$cex * 1.33, xpd=TRUE)

# Make dummy plot
plot(c(0, 4), c(0, 1), type="n", xlab="", ylab="", axes = F)

# Add interval lines and interval text
abline(v = c(0.5, 1.5, 2.5), lty = 2)

text(1, 0.01, expression(paste("Interval ", d - 1)), adj = .5)
text(2, 0.01, expression(paste("Interval ", d)), adj = .5)

#####
# Setup for the observations
n_series = 7
y_pos = seq(0, 1, length.out = n_series + 2)[-c(1, n_series +2)]

# Each element is an observation
# First column is the x coordinates
# Second column is the pch value
x_vals_and_point_codes <-
  list(
    structure(c(0.017182620242238, 1.3,
                4, 1),
              .Dim = c(2L, 2L)),
    structure(c(0.0665940539911389, 0.8, 1.9, 2.20536063500913,
                4, 4, 4, 16),
              .Dim = c(4L, 2L)),
    structure(c(0.406624712632038, 1.13612816265784,
                4, 16),
              .Dim = c(2L, 2L)),
    structure(c(0.092221238068305, 0.8, 2.2, 3, 3.7,
                4, 4, 4, 4, 1),
              .Dim = c(5L, 2L)),
    structure(c(0.285507099120878, 1, 2.12480141618289,
                4, 4, 16),
              .Dim = c(3L, 2L)),
    structure(c(2.15463116460014,2.33,
                4, 16), .Dim = c(2L, 2L)),
    structure(c(0.254905348294415,1.4, 2.1, 3.68174365991727,
                4, 4, 4, 16),
              .Dim = c(4L, 2L)))

# Add the observations to the plot
for(i in seq_along(x_vals_and_point_codes)){
  vals <- x_vals_and_point_codes[[i]]
  y = y_pos[i]
  xs = vals[, 1]
  n_xs = length(xs)

  # add lines
  segments(xs[-n_xs], rep(y, n_xs - 1),
           xs[-1], rep(y, n_xs - 1))

  # Add point
  points(xs, rep(y, n_xs), pch = vals[, 2],
         cex = ifelse(vals[, 2] == 1, par()$cex * 2.5, par()$cex))

  # Add numbers
  text(xs, rep(y + .05, n_xs), as.character(1:n_xs))

  # Add extra text to first observation
  if(i == length(x_vals_and_point_codes)){
    for(j in 1:(n_xs - 1)){
      text(xs[j], y + +.12,
           substitute(paste(bold(x)[ajp], ", (", t[aj], ",",  t[ajp], "]"),
                      list(aj = paste0("a", j - 1),
                           ajp = paste0("a", j))),
           cex = par()$cex * 1.2)
    }
  }
}

# Add letters
x <- sapply(x_vals_and_point_codes, "[", 1, 1)
x <- pmin(x - .1, .4)
text(rev(x), rev(y_pos), letters[1:n_series], cex = par()$cex * 1.5)

## ----cont_fit------------------------------------------------------------
ddfit_cont <- ddhazard(
  formula = frm,
  data = hd_dat_sub,
  model = "exp_clip_time_w_jump", # Change model from default
  by = 1,
  max_T = 60,
  id = hd_dat_sub$serial_number,
  Q_0 = diag(1, 23),
  Q = diag(0.01, 23),
  control = list(
    NR_eps = 0.001,  # EKF with extra iterations
    eps = 0.001,
    LR = .8,
    method = "EKF")) # Use EKF

## ----cont_plot, echo=FALSE, par_3x3 = TRUE, fig.cap="Predicted factor levels parameters with the model using the right clipped time variable with a jump term.",cache=FALSE----
ddfit_cont <- get_pretty_model_factors(ddfit_cont)
plot(ddfit_cont, cov_index = 1:9)

## ----cont_cleanup, echo = FALSE------------------------------------------
rm(ddfit_cont)

## ----run_sim_exp, echo = FALSE,cache = FALSE-----------------------------
# Although caching is used, I also make an additional copy of the results
# as this part takes a while and to avoid re-computations in case of some
# minor change in the code here or similar
 library(tcltk)

# Check if results are already computed
result_file <- "results_from_simulation.Rds"
sim_env_file <- "sim_env.Rds"
do_not_compute <-
  file.exists(result_file) &&
  file.exists(sim_env_file) &&
  (!tclvalue(tkmessageBox(
    title = "Rerun", message = "Want to rerun the simulation?",
    type = "yesno")) == "yes")

if(!do_not_compute){
  with(sim_env <- new.env(), {
    #####
    # Function to make sampling go quicker
    get_exp_draw <- with(environment(ddhazard), get_exp_draw())
    get_unif_draw <- with(environment(ddhazard), get_unif_draw())
    get_norm_draw <- with(environment(ddhazard), get_norm_draw())

    #####
    # Define simulation function
    sim_func <- function(
      n_series, n_vars = 10L, t_0 = 0L, t_max = 30L, cov_params = 1,
      re_draw = T, beta_start = rnorm(n_vars), intercept_start,
      sds = rep(1, n_vars + 1), run_n = 1){
      # Make output matrix
      n_row_max <- n_row_inc <- 10^5
      res <- matrix(
        NA_real_, nrow = n_row_inc, ncol = 4 + n_vars,
        dimnames = list(NULL, c("id", "tstart", "tstop", "event",
                                paste0("x", 1:n_vars))))
      cur_row <- 1

      if(re_draw){
        get_unif_draw(re_draw = T)
        get_exp_draw(re_draw = T)
        get_norm_draw(re_draw = T)
      }

      if(length(beta_start) == 1)
        beta_start <- rep(beta_start, n_vars)

      # draw betas
      betas <- matrix(get_norm_draw((t_max - t_0 + 1) * (n_vars + 1)),
                      ncol = n_vars + 1, nrow = t_max - t_0 + 1)
      betas <- t(t(betas) * sds)
      betas[1, ] <- c(intercept_start, beta_start)
      betas <- apply(betas, 2, cumsum)

      # covariate sim expression
      cov_exp <- expression(cov_params * get_norm_draw(n_vars))

      #####
      # Simulate
      for(id in 1:n_series){
        interval_start <- tstart <- tstop <-
          max(floor(get_unif_draw(1) *  2 * t_max) - t_max, 0L)
        repeat{
          tstop <- tstop + 5L
          if(tstop >= t_max)
            tstop <- t_max

          x_vars <- eval(cov_exp)
          l_x_vars <- c(1, x_vars) # add intercept

          tmp_t <- tstart
          while(tmp_t <= interval_start &&  interval_start < tstop){
            # Plus 2 for the coefficient vector at time zero and interval_start
            # is t - 1
            exp_eta <- exp(.Internal(drop(
              betas[interval_start + 2, ] %*% l_x_vars)))
            event <- exp_eta / (1 + exp_eta) > get_unif_draw(1)

            interval_start <- interval_start + 1L
            if(event){
              tstop <- interval_start
              break
            }

            tmp_t <- tmp_t + 1L
          }

          res[cur_row, ] <- c(id, tstart, tstop, event, x_vars)

          if(cur_row == n_row_max){
            # We need to add more rows to the ouput matrix
            n_row_max <- n_row_max + n_row_inc
            res = rbind(res, matrix(NA_real_, nrow = n_row_inc, ncol = 4 + n_vars))
          }
          cur_row <- cur_row + 1

          if(event || tstop >= t_max)
            break

          tstart <- tstop
        }
      }

      list(res = as.data.frame(res[1:(cur_row - 1), ]), betas = betas)
    }

    sim_func <- compiler::cmpfun(sim_func)

    #####
    # Define parameters
    # Some have to elements; one for each simulation experiment

    # General parameters
    n_series <- 2^(9:20)
    n_vars <- c(20, 20)
    t_max <- 30L
    intercept_start <- -3.5
    beta_sd <- c(.33, .33)
    intercept_sd <- 0.1
    Q_0_arg <- 1e5
    Q_arg <- 0.01
    denom_term <- 0.00001
    LR <- 1
    n_max <- 9
    eps <- 0.01
    cov_params <- list(c("sigma" = 1), c("sigma" = .33))

    # SMA and EKF parameters
    NR_eps = 0.1
    SMA_meth = "woodbury"
    Q_0_small <- 1

    # Number of trials and seeds
    n_sims <- 11
    set.seed(4368560)
    seeds <- sample.int(n_sims)

    # UKF parameters
    ukf_alpha <- 1
    ukf_w0 <- 0.0001
    ukf_beta <- 0
    ukf_max <- 2^18
    ukf_kappa <- (2 * n_vars + 1) * (1 + ukf_alpha^2 * (ukf_w0 - 1)) / (ukf_alpha^2 * (1 - ukf_w0))
    ukf_Q_0 <- .01

    # Sanity check
    m <- 2 * n_vars + 1
    lambda <- ukf_alpha^2 * (m + ukf_kappa) - m
    stopifnot(all.equal(lambda / (m + lambda), rep(ukf_w0, 2)))
    rm(m, lambda)

    # Number of threads to use
    n_threads <- max(parallel::detectCores() - 1, 2)
    options(ddhazard_max_threads = n_threads)

    # Define result array
    results <- array(
      NA_real_, dim = c(2, length(seeds), length(n_series), 5, 3),
      dimnames = list(
        NULL, NULL, NULL,
        c("EKF", "EKFx", "UKF", "SMA", "GMA"),
        c("elapsed", "MSE", "niter")))

    #####
    # Function to get estimates
    get_fit <- eval(bquote(
      function(data, method, run_n = 1){
        gc() # Make sure garbage collection is run before

        try({
          time <- system.time(fit <- ddhazard(
            Surv(tstart, tstop, event) ~ . - tstart - tstop - event - id,
            data = data, by = 1L, max_T = .(t_max), id = data$id,
            Q_0 = diag(
              ifelse(method %in% c("GMA", "EKFx"),
                     .(Q_0_small),
                     ifelse(method == "UKF" ,
                            .(ukf_Q_0), .(Q_0_arg))),
              .(n_vars)[run_n] + 1),
            Q = diag(.(Q_arg), .(n_vars)[run_n] + 1),
            control = list(
              LR = if(method == "EKF") .(LR) else 1,
              method = stringr::str_replace(method, "x", ""),
              alpha = .(ukf_alpha), beta = .(ukf_beta),
              kappa = .(ukf_kappa)[run_n], denom_term = .(denom_term),
              save_risk_set = F, save_data = F,
              LR_max_try = 1, # we only try one learning rate
              n_max = .(n_max),
              eps = .(eps),
              posterior_version = .(SMA_meth),
              NR_eps =  if(method == "EKFx") NR_eps else NULL
            )))["elapsed"]

          return(list(fit = fit, time = time))
        })

        return(NULL)
      }))

    #####
    # Function to get MSE
    mse_func <- function(betas, fit)
      mean((betas[-1, ] - fit$state_vecs[-1, ])^2) # remove the first entry which
                                                   # is just a same as period
                                                   # one estimate

    ######
    # Run experiment
    for(run_n in 1:2){
      for(i in seq_along(seeds)){
        s <- seeds[i]
        for(j in seq_along(n_series)){
          print(paste0("Using seed ", i, " with number of series index ", j,
                       " in run ", run_n))
          n <- n_series[j]
          set.seed(s)

          # Simulate
          sims <- sim_func(
            n_series = n, n_vars = n_vars[run_n], t_max = t_max,
            intercept_start = intercept_start, cov_params = cov_params[[run_n]],
            sds = c(intercept_sd, rep(beta_sd[run_n], n_vars[run_n])),
            run_n = run_n)

          # EKF
          out <- get_fit(data = sims$res, "EKF", run_n)
          if(!is.null(out)){
            results[run_n, i, j, "EKF", "elapsed"] <- out$time
            results[run_n, i, j, "EKF", "MSE"] <- mse_func(sims$betas, out$fit)
            results[run_n, i, j, "EKF", "niter"] <- out$fit$n_iter
          }

          # EKFx
          out <- get_fit(data = sims$res, "EKFx", run_n)
          if(!is.null(out)){
            results[run_n, i, j, "EKFx", "elapsed"] <- out$time
            results[run_n, i, j, "EKFx", "MSE"] <- mse_func(sims$betas, out$fit)
            results[run_n, i, j, "EKFx", "niter"] <- out$fit$n_iter
          }

          # UKF
          if(n <= ukf_max){
            out <- get_fit(data = sims$res, "UKF", run_n)
            if(!is.null(out)){
              results[run_n, i, j, "UKF", "elapsed"] <- out$time
              results[run_n, i, j, "UKF", "MSE"] <- mse_func(sims$betas, out$fit)
              results[run_n, i, j, "UKF", "niter"] <- out$fit$n_iter
            }
          }

          # SMA
          out <- get_fit(data = sims$res, "SMA", run_n)
          if(!is.null(out)){
            results[run_n, i, j, "SMA", "elapsed"] <- out$time
            results[run_n, i, j, "SMA", "MSE"] <- mse_func(sims$betas, out$fit)
            results[run_n, i, j, "SMA", "niter"] <- out$fit$n_iter
          }

          # GMA
          out <- get_fit(data = sims$res, "GMA", run_n)
          if(!is.null(out)){
            results[run_n, i, j, "GMA", "elapsed"] <- out$time
            results[run_n, i, j, "GMA", "MSE"] <- mse_func(sims$betas, out$fit)
            results[run_n, i, j, "GMA", "niter"] <- out$fit$n_iter
          }

          print(results[run_n, i, j,,])
        }

        print(results[run_n, i,,,])
        rm(out, sims)
      }
    }
  })

  # Take copy of results and cleanup
  results <- sim_env$results
  rm(i, j, s, get_fit, mse_func, envir = sim_env)

  # Save for later
  saveRDS(sim_env, file = sim_env_file)
  saveRDS(results, file = result_file)
} else {
  sim_env <- readRDS(sim_env_file)
  results <- readRDS(result_file)
}

## ----sim_coefficients_ex, echo=FALSE, results="hide", fig.cap = "Example of parameters in the simulation experiment. The black curve is the intercept and the gray curves are the parameters for the covariates.", par_1x1 = TRUE, cache = FALSE----
# Plot one example
with(sim_env, {
  sims <- sim_func(
          n_series = 100, n_vars = n_vars[1], t_max = t_max,
          intercept_start = intercept_start,
          sds = c(intercept_sd, rep(beta_sd[1], n_vars[1])))

  matplot(sims$betas, type = "l", lty = 1,
          col = c("black", rep("gray40", n_vars[1])),
          xlab = "Time", ylab = "Parameter")
})

## ----tbl_sim_stats, results='asis',echo=FALSE,cache=FALSE----------------
# Compute means and medians
medians <- apply(results, c(1, 3:5), median, na.rm = T)
means <- apply(results, c(1, 3:5), mean, na.rm = T)

# Define output table
tbl_sum <- matrix(
  NA_real_, nrow = 2, ncol = dim(results)[4],
  dimnames = list(
    c("Run time", "Log-log slope"),
    dimnames(results)[[4]]))

# Enter run times
tbl_sum["Run time", ] <-
  apply(medians[1, , , "elapsed"], 2, max, na.rm = T)

# Compute log-log regression slope
log_reg_cut_off <- 2^14
for(n in dimnames(tbl_sum)[[2]]){
  tmp <- results[1, , , n, "elapsed"]
  tbl_sum["Log-log slope", n] <-
    lm(log(c(t(tmp[, log_reg_cut_off<= sim_env$n_series]))) ~
       log(rep(sim_env$n_series[sim_env$n_series >= log_reg_cut_off],
               sim_env$n_sims)))$coefficient[2]
}

# Change column name for EKFx
colnames(tbl_sum)[colnames(tbl_sum) == "EKFx"] <- "EKF with extra iterations"

# Make table and print
xtable(tbl_sum, digits = getOption("digits"),
       caption = paste0("Summary information of the computation time in the simulation study. The first row prints the median run time for largest number of individuals. The UKF is only up to $n=",
                        sim_env$ukf_max, "$. The second row shows the slope of the log computation time regressed on the log number of individuals for $n\\geq ",
                        log_reg_cut_off, "$."),
       label = "tab:runSummaryStats")

## ----sim_labels_setup_n_more, echo = FALSE-------------------------------
#####
# Plot settings
pchs <- structure(
  c(15, 4, 16, 17, 5),
  names = dimnames(medians)[[3]])
col_medians <- "black"
col_means <- rgb(0, 0, 0, alpha = .5)

#####
# Legend function
.legend <- names(pchs)
.legend[.legend == "EKFx"] <- "EKF w/ extra"

draw_legend <- function(x, y = NULL)
  legend(
    x = x, y = y,
    bty = "n",
    pch = pchs,
    legend = .legend)

#####
# fig.cap
sim_fig_cap <-
  "The EKF is the filled squares, the EKF with extra iteration is the crosses, the UKF is the circles, the SMA is the triangles and the GMA is the open square."

## ----sim_comp_time, echo=FALSE, results="hide", fig.cap = paste("Median computation times of the simulations for each method for different values of $n$. The gray symbols to the right are the means.", sim_fig_cap, "The scales are logarithmic so a linear trend indicates that computation time is a power of $n$."), par_1x1 = TRUE, cache = FALSE----
# Plot for computation time
# I save it as an epxression to re-use it
marg <- 1
time_plot_exp <- expression({
  par(xpd=TRUE)
  with(sim_env, {
    # Medians points
    matplot(
      n_series, medians[marg,,, "elapsed"], log  = "xy", col = col_medians,
      pch = pchs, type = "p", xaxt='n',
      xlab = "Number of individuals", ylab = "Computation time (seconds)")
    axis(1, at = sim_env$n_series)
    draw_legend("topleft")

    # Lines between points
    matlines(n_series, medians[marg,,, "elapsed"], lty = 2, col = col_medians)

    # Mean points
    matplot(
      n_series * 1.15, means[marg,,, "elapsed"], col = col_means,
      pch = pchs, type = "p", add = T)
  })})

eval(time_plot_exp)

## ----sim_MSE, echo=FALSE, results="hide", fig.cap = paste("Median mean square error of predicted parameters of the simulations for each method for different values of $n$. The gray symbols to the right are the means.", sim_fig_cap, "The axis are on the logarithmic scale."), par_1x1 = TRUE, cache = FALSE----
# Plot for MSE
marg <- 1

plot_exp <- expression({
  par(xpd=TRUE, yaxs = "i")
  with(sim_env, {
    # Median points
    matplot(
      n_series, medians[marg,,, "MSE"], log  = "xy", col = col_medians,
      pch = pchs, xaxt='n',
      xlab = "Number of individuals", ylab = "MSE of predicted parameters",
      ylim = c(min(medians[marg,,, "MSE"], na.rm = T),
               max(medians[marg,,, "MSE"], na.rm = T)))
    axis(1, at = sim_env$n_series)

    # Lines between points
    matlines(n_series,medians[marg,,, "MSE"], lty = 2, col = col_medians)

    # Mean points
    matplot(
      n_series * 1.15, means[marg,,, "MSE"], col = col_means,
      pch = pchs, type = "p", add = T)

    draw_legend("bottomleft")
  })
})

eval(plot_exp)

## ----n_iter_plot, echo = FALSE, fig.cap = paste0("Median number of iterations of the EM-algorithm. ", sim_fig_cap), par_1x1= TRUE, cache = FALSE----
par(xpd=TRUE)

with(sim_env, {
  # Median points
  matplot(n_series, medians[1,,, "niter"],
          pch = pchs, xaxt='n',
          ylim = c(0, 10.5), log  = "x",
          col = col_medians,
          xlab = "Number of individuals", ylab = "Iterations of the EM")
  axis(1, at = sim_env$n_series)

  # Lines between points
  matplot(n_series, medians[1,,, "niter"],
          type = "l", lty = 2,
          col = col_medians, add = T)

  draw_legend("bottomleft")
})

## ----hist_of_lp_dens, echo = FALSE, par_1x1=TRUE, fig.cap="Estimated density for the linear predictor in the last interval in the first simulation experiment."----
# The density is not Gaussian. See https://math.stackexchange.com/a/42761/253239
tmp_env <- new.env(parent = sim_env)
with(tmp_env,{
  set.seed(6339855)

  lps <- replicate(1e3, {
    # Draw state at time 30
    n_vars <- n_vars[1]
    a0 <- c(intercept_start, eval(formals(sim_func)$beta_start))
    a30 <- mapply(
      rnorm, n = 1, mean = a0,
      sd = sqrt(t_max) * c(intercept_sd, rep(beta_sd[1], n_vars[1])))

    # Draw linear predictors
    n <- 1e3
    q <- n_vars
    X <- rnorm(n * (q + 1), sd = formals(sim_func)$cov_params)
    dim(X) <- c(q + 1, n)

    drop(a30 %*% X)
    })

  # Compute density estimate
  # hist(lp, breaks = 50)
  dens <- density(c(lps))
  plot(dens$x, dens$y, type = "l", ylab = "Density",
       main = "", xlab = "Linear predictor", xlim = range(dens$x), yaxs="i")
})

# Cleanup
rm(tmp_env)

## ----define_altered_plot_text, echo = FALSE------------------------------
altered_plot_text <- paste0(
  "but where each element of the covariate vectors are drawn from $N\\Lparen{0, ", sim_env$cov_params[[2]]["sigma"], "^2}$.")

## ----sim_comp_time_altered, echo=FALSE, results="hide", fig.cap = paste("Similar plot to figure~\\ref{fig:sim_comp_time}", altered_plot_text), par_1x1 = TRUE, cache = FALSE----
marg <- 2
eval(time_plot_exp)

## ----sim_MSE_altered, echo=FALSE, results="hide", fig.cap = paste("Similar plot to figure~\\ref{fig:sim_MSE}", altered_plot_text), par_1x1 = TRUE, cache = FALSE----
marg <- 2
eval(plot_exp)
