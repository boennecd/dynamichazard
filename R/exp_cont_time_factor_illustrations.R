options(digits = 4)

##############
# Function 1: original formulation

g_org <- function(e, a_t){
  nominator <- exp(a_t * exp(e) + e) * (1 - exp(a_t * exp(e)) + a_t * exp(e))
  denominator <- -1 * (1 - exp(2 * a_t * exp(e)) + 2 * a_t * exp(a * exp(e) +e))

  c("nominator" = nominator, "denominator" = denominator,
    "result" = nominator / denominator)
}


##############
# Function 2: using same factors and inverse of factors
g_prev_1 <- function(e, a_t){
  exp_eta <- exp(e)
  inv_exp_eta <- exp_eta^-1

  exp_term <- exp(a_t * exp(e))
  inv_exp_term <- exp_term^-1

  nominator <- inv_exp_term * (1.0 + a_t * exp_eta) - 1
  denominator <- inv_exp_eta * (1.0 -  inv_exp_term * inv_exp_term)  - 2 * a_t * inv_exp_term

  c("nominator" = nominator, "denominator" = denominator,
    "result" = nominator / denominator)
}


##############
# Function 3: computing all exponential separately
g1 <- function(e, a_t){
  nominator <- exp(-a_t * exp(e)) + a_t * exp(e - a_t * exp(e)) - 1
  denominator <- exp(-e) - exp(-2 * a_t * exp(e) - e) - 2 * a_t * exp(-a_t * exp(e))

  c("nominator" = nominator, "denominator" = denominator,
    "result" = nominator / denominator)
}


##############
# Function 4: using Laurent series
# See: https://www.wolframalpha.com/input/?i=(1+%2B+v+-+exp(v))%2F(exp(v)+-++exp(-v)+-+2v)
g2 <- function(e, a_t){
  v = a_t * exp(e)
  exp(e) * (- 3 /(2 * v) - 1 / 2 - v / 20 + v^3 / 8400)
}




########
# Example 1: small eta

# Good
a <- 1
eta <- -3

g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)

# Less well
a <- 1e-4
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)

# very bad
a <- 1e-8
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)




########
# Example 2: abs small eta

# Goes well
a <- 1
eta <- 1e-8

g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)

# Ok
a <- 1e-4
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)

# Bad?
a <- 1e-8
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)



########
# Example 3: large eta

# Good
a <- 1
eta <- 3

g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a) # <-- bad though!

# Good
a <- 1e-4
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)

# Bad
a <- 1e-8
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)


# Example from debugging
a <- 1
eta <- 3.18634
g_org(eta, a)
g_prev_1(eta, a)
g1(eta, a)
g2(eta, a)

#
eps <- 1
for(eta in -3:3 + 1e-8){
  a <- eps / exp(eta)
  print(c("eta" = eta, "alternative" = g1(eta, a)[3], "Lauren" = g2(eta, a)))
}

eps <- 1e-3
for(eta in -3:3 + 1e-8){
  a <- eps / exp(eta)
  print(c("eta" = eta, "alternative" = g1(eta, a)[3], "Lauren" = g2(eta, a)))
}

options("scipen"=100, "digits"=4)

grid_vals <- expand.grid(
  eps = 10^(0:(-6)),
  eta = -5:5
)

tmp <- t(mapply(function(eps, eta){
  a <- eps / exp(eta)
  c("eps" = eps, "eta" = eta, "alternative" = g1(eta, a)[3], "Lauren" = g2(eta, a))
}, eta = grid_vals$eta, eps = grid_vals$eps))

tmp[order(tmp[, "eps"]), ]

#############
# Score term from binary

f1 <- function(e, a)
  a * exp(e) / (1 - exp(- a * exp(e)))

f2 <- function(e, a){
  v <- a * exp(e)
  1 + v / 2 + v^2 / 12 - v^4/720
}

# This is not an issue
eta <- 1
a <- 1
f1(eta, a) # <-- Should be right
f2(eta, a)

# This is a small issue
eta <- 1
a <- 1e-14
f1(eta, a) # <-- Should be right
f2(eta, a)

# This is a large issue
eta <- -40
a <- 1
f1(eta, a) # <-- Should be right
f2(eta, a)


# Taylor series from https://www.wolframalpha.com/input/?i=v%2F(1-exp(-v))
v <- 10^(0:-14)
cbind("v" = v,
      "v/(1 - exp(-v))" = v / (1 - exp(-v)),
      "Taylor" = 1 + v * (1 / 2 + v * (1/12 - v * v /720)))


#############
# Variance factor from time variable

f1 <- function(e, a){
  (1 + a * exp(e) - exp(a * exp(e)))^2 /
    (exp(2 * a * exp(e)) - 1 - 2 * a * exp(e + a * exp(e)))
}

# Issue
f1(10, 10)
f1(10, 1)
f1(-10.8169, 0.291538)

# Not issue
f1(1, 10)
f1(1, 1)
f1(1, 1e-10)
f1(1e-10, 1e-10)
f1(1e-14, 1e-10)
f1(1e-14, 1e-14)


# More numerical examples
# from https://www.wolframalpha.com/input/?i=(1+%2B+v+-+exp(v))%5E2+%2F+(exp(2*v)+-+1-+2+*+v*+exp(v))
v <- 10^(14:-14)
cbind("v" = v,
      "(1 + v - exp(v))^2/(exp(2*v) - 1 - 2 * v* exp(v))" =
        (1 + v - exp(v))^2/(exp(2*v) - 1 - 2 * v* exp(v)),
      "Taylor" =
        v * (3 / 4 - v * (1 / 4 - v * (11 / 240 - v * (1 / 240 - v / 16800)))))

# Have to keep this in mind
(v <- 0.291538 * exp(-10.8169))
(1 + v - exp(v))^2/(exp(2*v) - 1 - 2 * v* exp(v))

#########
# Variance factor for binary
v <- 10^(-2:-20)
cbind("v" = v,
      "v^2 * exp(-v)/(1 - exp(-v))" =
        v^2 * exp(-v)/(1 - exp(-v)),
      "Taylor" =
        v - v^2 / 2 + v^3 / 12 - v^5 / 720)

#########
# Survival prob
# Source: https://www.wolframalpha.com/input/?i=1-exp(-v)
v <- 10^(-2:-16)
cbind("v" = v,
      "1 - exp(-v)" =
        1 - exp(-v),
      "Taylor" =
        v - v^2 / 2 + v^3 / 6 - v^4 / 24 + v^5 / 120)


#########
# Expected truncated survival time
# Source: https://www.wolframalpha.com/input/?i=(1-exp(-v))v%5E-1
f1 <- function(e, a)
  (1 - exp(-exp(e) * a)) / exp(e)

f_taylor <- function(e, a){
  v <- a * exp(e)

  a * ((1 - v / 2 * (1 - v / 6 * (1 - v / 24 * (1 - v / 120)))))
}

grid_vals <- expand.grid(e = -1:-6, a = 10^(-2:-8))

tmp <- t(mapply(function(e, a)
  c("v" = a * exp(e), "a" = a, "eta" = e,
    "Org formula" = f1(e, a), "Taylor" = f_taylor(e, a)),
  e = grid_vals$e, a = grid_vals$a))

tmp[sort(tmp[, "eta"]), ]





