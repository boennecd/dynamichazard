#ifndef EXP_MODEL_FUNCS
#define EXP_MODEL_FUNCS

#include <cmath>

namespace exp_model_funcs {
// Namespace to avoid namespace pollution
// By definition:
//  eps                 epsilon in ridge regression like
//  a                   at risk length
//  eta                 linear predictor x^T * beta
//  v                   a * exp(eta)
//  exp_x               exp(x)
//  inv_x               inv(x)
//  expect_chance_die   1 - exp(- v)

inline double expect_time(const double v, const double a,
                          const double inv_exp_v, const double exp_eta){
  return((v >= 1e-8) ? (1.0 - inv_exp_v) / exp_eta :
           a + exp_eta*(-pow(a,2)/2 + exp_eta*(pow(a,3)/6 - (pow(a,4)*exp_eta)/24))
  );
}

inline double expect_time_w_jump(const double exp_eta, const double inv_exp_eta,
                                 const double inv_exp_v, const double a){
  return((a * exp_eta >= 1e-4) ?
           (inv_exp_v - 1)*(-1 + a*exp_eta)*inv_exp_eta :
           a - (3*pow(a,2)*exp_eta)/2);
}

inline double expect_chance_die(const double v, const double inv_exp_v){
  return((v >= 1e-5) ? 1.0 - inv_exp_v :
           v * (1.0 - v / 2.0 * (1.0 + v / 6.0 * (1.0 - v / 24 * (1.0 - v /120.0)))));
}

inline double inv_var_wait_time(const double v, const double exp_eta, const double inv_exp_v){
  return((v >= 1e-4) ?
           // Set v = delta * exp(eta)
           // Then: exp(2eta) / (1 - exp(-2 delta * exp(eta)) - 2 exp(-delta * exp(eta)) delta exp(eta)) =
           //               exp(2eta) / (1 - exp(-2v) - 2 * v * exp(-v))
           exp_eta * exp_eta / (1.0 - inv_exp_v * inv_exp_v - 2.0 * v * inv_exp_v) :
           // Laruent series from https://www.wolframalpha.com/input/?i=1%2F(1-exp(2v)-2v*exp(v))
           exp_eta * exp_eta *
             (-1 / v * (1 / 4 - v * (1 / 4 - v * (5 / 48 - v * (1/48 - v /1440)))))
  );
}

inline double var_wait_time(const double v, const double a,  const double exp_eta, const double inv_exp_v){
  // exp(eta)^(-2) * (1 - exp(-2v) - 2 * v * exp(-v))
  // v = a * exp(eta) => ... = a^2 * v^(-2) * (1 - exp(-2v) - 2 * v * exp(-v))
  //
  // Use Taylor series for the latter when v is small: https://www.wolframalpha.com/input/?i=(1-exp(-2v)-2v*exp(-v))%2F(v%5E2)

  return((v >= 1e-4) ?
           (1.0 - inv_exp_v * inv_exp_v - 2.0 * v * inv_exp_v) / (exp_eta * exp_eta) :
           a * a * v * (1/3 - v * (1/3 - v * (11/60 - v * (13/180 - v * 19/840)))));
}

inline double var_wait_time_w_jump(const double exp_eta,
                                   const double inv_exp_v, const double a){
  return((a * exp_eta >= 1e-4) ?
           (1 - pow(inv_exp_v,2) + pow(a*exp_eta,2)*(3*inv_exp_v - pow(inv_exp_v,2)) +
           a*exp_eta*(-4*inv_exp_v + 2*pow(inv_exp_v,2)))/ pow(exp_eta,2) :
           exp_eta*((7*pow(a,3))/3 + exp_eta*((-19*pow(a,4))/6 + exp_eta*(
               (34*pow(a,5))/15 - (407*pow(a,6)*exp_eta)/360)))
  );
}

inline double var_chance_die(const double v, const double inv_exp_v){
  // Taylor series from https://www.wolframalpha.com/input/?i=exp(-v)+*+(1+-+exp(-v))
  return((v >= 1e-4) ?
           inv_exp_v * (1 - inv_exp_v) :
           v * (1 - v * (3/2 - v * (7/6 - v * (5/8 - v * 31 /120)))));
}
}

#endif
