#include <cmath>

#ifndef EXP_MODEL_FUNCS
#define EXP_MODEL_FUNCS

namespace exp_model_funcs {
// Namespace to avoid namespace pollution and avoid error-40 for Taylor/Laruens series
// By definition:
//  eps                 epsilon in ridge regression like solution
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

inline double inv_var_chance_die(const double v, const double inv_exp_v){
  return((v >= 1e-4) ?
           // Set v = a exp(eta)
           // Then: 1 / exp(- a exp(eta))(1 - exp(-a exp(eta))) = 1 / exp(-v) (1 - exp(-v))
           1 / (inv_exp_v * (1 - inv_exp_v)) :
           //Lauren series from https://www.wolframalpha.com/input/?i=1%2F((1-exp(-v))exp(-v))
           1 / v * (1 + v * (3 / 2 + v * (13 / 12 + v * (1 / 2 + v * 119 / 720))))
  );

}

inline double EKF_fac_score_die(const double exp_eta, const double v,
                                const double exp_v, const double a,
                                const double eps){
  if(a * exp_eta >= 1e-4){
    return(-(((-1 + exp_eta + exp_v*(-1 + a*(-2 + exp_eta))*exp_eta +
           pow(exp_v,2)*(1 + eps*pow(exp_eta,2)))*v)/
             (-1 + pow(exp_v,2)*(-1 + 2*eps*v - eps*pow(exp_eta,2)) -
               pow(exp_v,3)*eps*(1 + eps*pow(exp_eta,2)) + exp_v*(2 + eps + pow(v,2) +
               eps*pow(exp_eta,2)))));
  } else {
    // Taylor expansion around exp_eta = 0
    const double v_eps_ratio = v / eps;

    return(v_eps_ratio * (
        1 + v_eps_ratio * (
            (-2 + a - 2*eps) / 2 - v_eps_ratio *
              (-12 + 6*a - 3*pow(a,2) + 2*pow(a,3) - 30*eps + 14*a*eps - 6*pow(eps,2)) / 12)));
  }
}

inline double EKF_fac_score_time(const double exp_eta, const double v,
                                 const double exp_v, const double a,
                                 const double eps){
  if(a * exp_eta >= 1e-3){
    return(-(((-1 + exp_v - v)*(1 - exp_eta + eps*exp_eta*pow(exp_v,2) + exp_v*(-1 + exp_eta + v)))/
           (1 - (2 + eps + pow(v,2) + eps*pow(exp_eta,2))*exp_v + eps*(1 + eps*pow(exp_eta,2))*pow(exp_v,3) +
             pow(exp_v,2)*(1 + eps*pow(exp_eta,2) - 2*eps*v))));
  } else{
    const double exp_eta_eps_ratio = exp_eta / eps;

    return(exp_eta_eps_ratio * (-pow(a,2)/2 + exp_eta_eps_ratio *
           (-3*pow(a,4) + 2*pow(a,5) + 4*pow(a,3)*eps)/12));
  }
}

inline double EKF_fac_var(const double exp_eta, const double v,
                          const double exp_v, const double a,
                          const double eps){
  if(a * exp_eta >= 1e-3){
    return((-1 + exp_v*(3 + pow(v,2) + eps*pow(exp_v,3) + exp_v*(-3 + eps + pow(v,2)*(-1 + eps) +
           pow(v,2)*eps*pow(exp_eta,2)+ 2*eps*v) + pow(exp_v,2)*(1 - 2*eps*(1 + v))))/
             (exp_v - (2 + eps + (pow(a,2) + eps)*pow(exp_eta,2))*pow(exp_v,2) + (1 + eps*exp_eta*(-2*a + exp_eta))*
               pow(exp_v,3) + eps*(1 + eps*pow(exp_eta,2))*pow(exp_v,4)));
  } else{
    return((4*pow(a,2) + pow(a,4))*exp_eta * (exp_eta / (4 * eps)));
  }
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

inline double covar(const double v,
                    const double inv_exp_v, const double exp_eta){
  // Taylor series from http://www.wolframalpha.com/input/?i=-1+*+exp(-2v)+*+(1+%2B+v+*+exp(v)+-+exp(v))
  return(
    (v >= 1e-4) ?
  -1 * inv_exp_v * (inv_exp_v + v - 1.) / exp_eta :
  -v * v * (1/2 - v * (2/3 - v * (11/24 - v * (13/60 - 19 * v / 240)))) / exp_eta);
}

inline double binary_score_fac(const double v, const double  inv_exp_v, double eps){
  return((inv_exp_v*v)/(eps + (1 - inv_exp_v)*inv_exp_v));
}

inline double binary_var_fac(const double v, const double  inv_exp_v, double eps){
  return(pow(inv_exp_v*v, 2)/(eps + (1 - inv_exp_v)*inv_exp_v));
}


inline double clip_time_score_fac(const double exp_eta, const double inv_exp_eta,
                                  const double inv_exp_v, const double a,
                                  const double eps){

  return((a * exp_eta >= 1e-4) ?
           (-1 + inv_exp_v + a*exp_eta*inv_exp_v)/(eps*exp_eta + inv_exp_eta - 2*a*inv_exp_v -
           inv_exp_eta*pow(inv_exp_v,2)) :
           ((exp_eta / eps) *(-3*pow(a,2)*eps + (pow(a,5) + 2*pow(a,3)*eps)*exp_eta))/(6*eps));
}

inline double clip_time_var_fac(const double exp_eta, const double inv_exp_eta,
                                const double inv_exp_v, const double a,
                                const double eps){
  return((a * exp_eta >= 1e-4) ?
           (1 - 2*inv_exp_v - 2*a*exp_eta*inv_exp_v + pow(inv_exp_v,2) + 2*a*exp_eta*pow(inv_exp_v,2) +
           pow(a*exp_eta*inv_exp_v, 2))/
             (1 + eps*pow(exp_eta,2) - 2*a*exp_eta*inv_exp_v - pow(inv_exp_v,2)) :
           (pow(exp_eta/eps,2)*(3*pow(a,4)*eps + (-pow(a,7) - 4*pow(a,5)*eps)*exp_eta))/12);
}

inline double clip_time_w_jump_score_fac(const double exp_eta, const double v,
                                         const double inv_exp_v, const double a,
                                         const double eps){
  if(v >= 1e-4){
    return((exp_eta*(-1 + inv_exp_v) + a*pow(exp_eta,2)*inv_exp_v - pow(a*exp_eta,2) *exp_eta*inv_exp_v)/
           (1 - pow(inv_exp_v,2) + a*exp_eta*(-4*inv_exp_v + 2*pow(inv_exp_v,2)) +
             pow(exp_eta,2)*(eps + pow(a,2)*(3*inv_exp_v - pow(inv_exp_v,2)))));
  } else{
    return((pow(v,2)*(-18*pow(a,2) - 9*eps + (7*pow(a,2) + 8*eps)*v))/(24*pow(a,4)*exp_eta + 24*pow(a,2)*eps*exp_eta +
           6*pow(eps,2)*exp_eta));
  }
}

inline double clip_time_w_jump_var_fac(const double exp_eta, const double v,
                                       const double inv_exp_v, const double a,
                                       const double eps){
  if(v >= 1e-4){
    return((1 - 2*inv_exp_v + pow(inv_exp_v,2) + (- 2*pow(a*exp_eta,3) + pow(a*exp_eta,4))*pow(inv_exp_v,2) +
           pow(a*exp_eta,2)*(2*inv_exp_v - pow(inv_exp_v,2)) + a*exp_eta*(-2*inv_exp_v + 2*pow(inv_exp_v,2)))/
             (1 - pow(inv_exp_v,2) + a*exp_eta*(-4*inv_exp_v + 2*pow(inv_exp_v,2)) +
               pow(exp_eta,2)*(eps + pow(a,2)*(3*inv_exp_v - pow(inv_exp_v,2)))));
  } else{
    return(9*(v/eps)*pow(v/exp_eta, 2) * (v/4));
  }
}}

#endif
