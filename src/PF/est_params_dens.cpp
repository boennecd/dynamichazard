#include "est_params.h"
#include "../parallel_qr.h"

using particle_pairs = smoother_output::particle_pairs;
using pair_iterator = std::vector<const particle_pairs*>::const_iterator;

class update_parameters_data {
  using trans_like = std::vector<std::vector<particle_pairs>>;

  /* function to create a list of pointers to be used in loops later */
  std::vector<const smoother_output::particle_pairs*> set_pairs(){
    unsigned long int total_n_pairs = 0;
    for(auto i = tr.begin(); i != tr.end(); ++i)
      total_n_pairs += i->size();

    std::vector<const particle_pairs*> out(total_n_pairs);
    auto ptr = out.begin();
    for(auto i = tr.begin(); i != tr.end(); ++i)
      for(auto j = i->begin(); j != i->end(); ++j, ++ptr)
        *ptr = &(*j);

    return out;
  }

public:
  const trans_like &tr;
  select_mapper R;
  const unsigned int n_periods;
  const unsigned int n_elem_X;
  const unsigned int n_elem_Y;

  const std::vector<const particle_pairs*> pairs;

  update_parameters_data
    (const smoother_output &sm_output, const arma::mat &R):
    tr(*sm_output.get_transition_likelihoods().get()),
    R(R), n_periods(tr.size()),
    n_elem_X(tr[0][0].p->get_state().n_elem),
    n_elem_Y(this->R.map(tr[0][0].p->get_state(), trans).sv.n_elem),
    pairs(set_pairs())
    {}
};

/* generator for chunks for QR decomposition */
class generator : public qr_data_generator {
  const update_parameters_data &dat;
  pair_iterator i_start;
  pair_iterator i_end;
  const unsigned long int n_pairs;

public:

  generator
    (const update_parameters_data &dat, pair_iterator i_start,
     pair_iterator i_end, const unsigned long int n_pairs):
    dat(dat), i_start(i_start), i_end(i_end), n_pairs(n_pairs) { }

  qr_work_chunk get_chunk() const override {
    // setup objects
    arma::mat  Y   (dat.n_elem_Y, n_pairs);
    arma::mat  X   (dat.n_elem_X, n_pairs);
    /* keep defaults to not include the pair */
    arma::uvec keep(              n_pairs, arma::fill::zeros);

    unsigned int i = 0;
    for(auto pairs = i_start; pairs != i_end; ++pairs){
      const particle_pairs &pas = **pairs; /* pairs iterator for a pointer */

      if(exp(pas.log_weight) < 1e-8){ /* don't keep */
        i += pas.transition_pairs.size();
        continue;
      }

      arma::vec Y_val = dat.R.map(pas.p->get_state(), trans).sv;
      for(auto pair = pas.transition_pairs.begin();
          pair != pas.transition_pairs.end();
          ++pair, ++i){
        double w_sqrt = exp((pas.log_weight +  pair->log_weight) / 2);
        if(w_sqrt < 1e-8) /* don't keep */
          continue;

        keep[i]  = 1L;
        Y.col(i) = w_sqrt * Y_val;
        X.col(i) = w_sqrt * pair->p->get_state();
      }
    }

    keep = arma::find(keep);
    Y = Y.cols(keep);
    X = X.cols(keep);

    return { std::move(X), std::move(Y), Y.t() * Y};
  }
};

PF_parameters
  est_params_dens
  (const smoother_output &sm_output, const arma::vec &a_0, const arma::mat &Q,
   const arma::mat &Q_0, const arma::mat &R, int max_threads,
   const unsigned long int max_bytes){
    update_parameters_data dat(sm_output, R);

    /* setup generators */
    const unsigned long int max_per_block =
        max_bytes / ((dat.n_elem_X + dat.n_elem_Y) * 8L);

    qr_parallel qr_calc(
        std::vector<std::unique_ptr<qr_data_generator>>(), max_threads);
    for(auto i_start = dat.pairs.begin(); i_start != dat.pairs.end(); ){
      pair_iterator i_end = i_start;
      unsigned long int n_pairs = 0;
      do {
        n_pairs += (*i_end)->transition_pairs.size();
        ++i_end;

      } while (n_pairs < max_per_block && i_end != dat.pairs.end());

      qr_calc.submit(std::unique_ptr<qr_data_generator>(
        new generator(dat, i_start, i_end, n_pairs)));
      i_start = i_end;
    }

    /* compute new estimates. The notation is a bit confusing. The output from
     * qr_parallel.compute() is the R and F in
     *   (R^\top R)^{-1} X = R^\top F
     *
     * where X is the new coef matrix which is PF_parameters.F. */
    R_F res = qr_calc.compute();
    PF_parameters out;

    /* TODO: needs to be updated for higher order models since we need the
     *       full F matrix in a bit */
    arma::mat R_from_QR = res.R_rev_piv();
    out.R_top_F = arma::solve(R_from_QR.t(), R_from_QR.t() * res.F,
                              arma::solve_opts::no_approx);
    out.R_top_F = arma::solve(R_from_QR    , out.R_top_F,
                              arma::solve_opts::no_approx);

    /* use that
     *    Q = Y^\top Y - F^\top F
     *
     * where F is the output from qr_parallel.compute().
     */
    out.Q = res.dev - res.F.t() * res.F;

    // update a_0
    /* TODO: needs to be updated for higher order models. Need Full F matrix */
    inv_mapper inv_map(out.R_top_F);
    const std::vector<particle_pairs> &first_particles = dat.tr[0];
    out.a_0 = arma::zeros<arma::vec>(a_0.n_elem);

    for(auto i = first_particles.begin(); i != first_particles.end(); ++i)
        out.a_0 += exp(i->log_weight) * i->p->get_state();

    return out;
  }
