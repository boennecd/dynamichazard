#ifndef QR_PARALLEL
#define QR_PARALLEL

#include "arma_n_rcpp.h"
#include <memory>

struct qr_work_chunk {
  arma::mat X;
  arma::mat Y;
  arma::mat dev; /* used different. Either Y^\top Y in the multivaraite case
                  * and the deviance (a scalar) in the GLM case */
};

class qr_data_generator {
public:
  virtual qr_work_chunk get_chunk() const = 0;
};

/* Let x = QR. Then this is a data holder for the R matrix and f = Q^\top y.
 * dev is deviance computed for the chunk computed with the current
 * coefficient vector                                                        */
struct R_F {
  const arma::mat R;
  const arma::uvec pivot;
  const arma::mat F;
  const arma::mat dev;

  arma::mat R_rev_piv() const;
};

class qr_parallel {
  using ptr_vec = std::vector<std::unique_ptr<qr_data_generator> >;

  class worker {
    std::unique_ptr<qr_data_generator> my_generator;

  public:
    worker(std::unique_ptr<qr_data_generator>);

    R_F operator()();
  };

  ptr_vec generators;
  const unsigned int max_threads;

public:
  qr_parallel(ptr_vec, const unsigned int);

  R_F compute();
};

#endif // QR_PARALLEL

