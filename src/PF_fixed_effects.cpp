#include "arma_n_rcpp.h"
#include "thread_pool.h"
#include "family.h"

static constexpr double zero_eps = 1e-100;

inline arma::vec set_offsets(const arma::vec &dts, const glm_base *family){
  arma::vec offsets;

  if(family->name() == POISSON){
    offsets = arma::log(dts);

  } else if(family->name() == BINOMIAL or family->name() == CLOGLOG){
    offsets = arma::vec(dts.n_elem, arma::fill::zeros);

  } else
    Rcpp::stop("family not implemented");

  return offsets;
}

class pf_fixed_it_worker {
  const arma::mat X;
  const arma::vec Y;
  const arma::vec dts;
  const arma::mat cloud;
  const arma::vec cl_weights;
  const arma::mat ran_vars;
  const arma::vec &fixed_params;
  std::unique_ptr<glm_base> family;

public:
  pf_fixed_it_worker(arma::mat &&X, arma::vec &&Y, arma::vec &&dts,
                     arma::mat &&cloud, arma::vec &&cl_weights,
                     arma::mat &&ran_vars, const arma::vec &fixed_params,
                     const std::string family, const bool debug):
  X(X), Y(Y), dts(dts), cloud(cloud), cl_weights(cl_weights),
  ran_vars(ran_vars), fixed_params(fixed_params),
  family(get_fam<glm_base>(family)) {
#ifdef DDHAZ_DEBUG
    auto throw_invalid_arg = [](const std::string &what){
      throw std::invalid_argument("pf_fixed_it_worker: invalid '" +
                                  what + "'");
    };

    const unsigned long n_obs = X.n_cols;
    if(n_obs < 1L)
      throw_invalid_arg("n_obs");
    if(X.n_rows != fixed_params.n_elem)
      throw_invalid_arg("X");
    if((cloud.n_cols != 0L and ran_vars.n_rows != cloud.col(0).n_elem) or
         ran_vars.n_cols != n_obs)
      throw_invalid_arg("ran_vars");
    if(Y.n_elem != n_obs)
      throw_invalid_arg("Y");
    if(dts.n_elem != n_obs)
      throw_invalid_arg("dts");
#endif
  }

  struct result {
    arma::mat Rmat;
    arma::vec XtWY;
  };

  result operator()(){
    arma::uword n = X.n_cols, p = X.n_rows;
    result res { arma::mat(), arma::vec(p, arma::fill::zeros) };
    arma::mat &Rmat = res.Rmat;
    arma::vec &XtWY = res.XtWY, ws_sum(n, arma::fill::zeros),
      res_sum(n, arma::fill::zeros);

    const arma::vec offset = set_offsets(dts, family.get());
    const arma::vec eta_wo_rng = X.t() * fixed_params;

    for(unsigned int i = 0; i < cloud.n_cols; ++i){
      double cl_w = cl_weights[i];
      if(cl_w <= 0)
        continue;

      arma::vec offset_it = ran_vars.t() * cloud.col(i) + offset;
      arma::vec eta = eta_wo_rng + offset_it;

      arma::vec mu(n), mu_eta_val(n);
      double *m = mu.begin(), *mev = mu_eta_val.begin(), *e = eta.begin();
      for(arma::uword j = 0; j < mu.n_elem; ++j, ++m, ++mev, ++e){
        *m = family->glm_linkinv(*e);
        *mev = family->glm_mu_eta(*e);
      }

      arma::uvec good = arma::find(
          ((-zero_eps < mu_eta_val) + (mu_eta_val < zero_eps) != 2));

      mu = mu(good);
      eta = eta(good);
      mu_eta_val = mu_eta_val(good);
      arma::vec var(mu.n_elem);
      m = mu.begin();
      for(auto v = var.begin(); v != var.end(); ++v, ++m)
        *v = family->glm_variance(*m);

      /* compute working responses and add to result */
      arma::vec z = (eta - offset_it(good)) + (Y(good) - mu) / mu_eta_val;
      arma::vec w = cl_w * arma::square(mu_eta_val) / var;

      double *pz = z.begin(), *pw = w.begin();
      for(auto j : good){
        ws_sum(j) += *pw;
        res_sum(j) += *(pz++) * *(pw++);
      }
    }

    arma::mat XtWsqrt = X;
    ws_sum.for_each( [](arma::vec::elem_type& val) { val = sqrt(val); } );
    arma::inplace_trans(XtWsqrt);
    XtWsqrt.each_col() %= ws_sum;

    /* QR factroize, get R, and pivot */
    QR_factorization qr(XtWsqrt);
    arma::uvec piv = qr.pivot();
    piv(piv) = arma::regspace<arma::uvec>(0, 1, piv.n_elem - 1);
    Rmat =  qr.R().cols(piv);

    for(arma::uword i = 0; i < n; ++i)
      XtWY += X.col(i) * res_sum(i);

    return res;
  }
};

// [[Rcpp::export]]
Rcpp::List pf_fixed_effect_get_QR(
    Rcpp::List clouds, Rcpp::List risk_obj, const arma::mat &ran_vars,
    const arma::mat &fixed_terms, const arma::mat &R_top,
    const arma::vec &tstart, const arma::vec &tstop,
    const arma::vec &fixed_params, const std::string family,
    const int max_threads, const bool debug){
  const unsigned int n_clouds = clouds.size();
  Rcpp::List risk_sets = Rcpp::as<Rcpp::List>(risk_obj["risk_sets"]);
  arma::ivec is_event_in = Rcpp::as<arma::ivec>(risk_obj["is_event_in"]);
  arma::vec event_times = Rcpp::as<arma::vec>(risk_obj["event_times"]);

  thread_pool pool(std::max((int)1L, std::min((int)n_clouds, max_threads)));
  std::vector<std::future<pf_fixed_it_worker::result> > futures;
  futures.reserve(n_clouds);

  for(unsigned int i = 0; i < n_clouds; ++i){
    Rcpp::List cl = Rcpp::as<Rcpp::List>(clouds[i]);
    arma::uvec risk_set = Rcpp::as<arma::uvec>(risk_sets[i]);
    if(risk_set.n_elem == 0L)
      continue;
    risk_set -= 1L;

    arma::mat ran_vars_i = ran_vars   .cols(risk_set);
    arma::mat X_i        = fixed_terms.cols(risk_set);
    arma::vec y_i;
    {
      arma::ivec tmp = is_event_in.elem(risk_set);
      tmp.for_each([i](arma::ivec::elem_type &val) { val = val == (int)i; } );
      y_i = arma::conv_to<arma::vec>::from(tmp);
    }

    arma::vec ws = Rcpp::as<arma::vec>(cl["weights"]);
    arma::uvec good = arma::find(ws >= 1e-7);

    double int_stop  = event_times[i + 1L],
           int_start = event_times[i     ];
    arma::vec sta = tstart.elem(risk_set), sto = tstop.elem(risk_set);
    sta.for_each([int_start](arma::vec::elem_type &val) {
      val = std::max(val, int_start); } );
    sto.for_each([int_stop ](arma::vec::elem_type &val) {
      val = std::min(val, int_stop ); } );
    arma::vec dts = sto - sta;

    ws = ws.elem(good);
    arma::mat particle_coefs = Rcpp::as<arma::mat>(cl["states"]);
    particle_coefs = R_top * particle_coefs.cols(good);

    pf_fixed_it_worker wo(
        std::move(X_i), std::move(y_i), std::move(dts),
        std::move(particle_coefs), std::move(ws), std::move(ran_vars_i),
        fixed_params, family, debug);
    futures.push_back(pool.submit(std::move(wo)));
  }

  const unsigned n_futures = futures.size();
  Rcpp::List out(n_futures);
  for(unsigned long i = 0; i < n_futures; ++i)
  {
    pf_fixed_it_worker::result res = futures[i].get();
    out[i] = Rcpp::List::create(
        Rcpp::Named("Rmat") = std::move(res.Rmat),
        Rcpp::Named("XtWY") = std::move(res.XtWY));
  }

  return out;
}
