#include "parallel_qr.h"
#include "arma_BLAS_LAPACK.h"

arma::mat R_F::R_rev_piv() const {
  arma::uvec piv = pivot;
  piv(piv) = arma::regspace<arma::uvec>(0, 1, piv.n_elem - 1);
  return R.cols(piv);
}

qr_parallel::worker::worker
  (std::unique_ptr<qr_data_generator> generator):
  my_generator(std::move(generator)) {}

R_F qr_parallel::worker::operator()(){
  qr_work_chunk my_chunk = my_generator->get_chunk();
  QR_factorization qr(my_chunk.X);
  arma::mat F = qr.qy(my_chunk.Y, true).rows(0, my_chunk.X.n_cols - 1);

  return R_F { qr.R(), qr.pivot(), std::move(F), my_chunk.dev };
}

qr_parallel::qr_parallel(
  ptr_vec generators, const unsigned int max_threads):
  pool(std::max((unsigned int)1L, max_threads))
  {
    futures.reserve(generators.size());
    while(!generators.empty()){
      submit(std::move(generators.back()));
      generators.pop_back();
    }
  }

void qr_parallel::submit(std::unique_ptr<qr_data_generator> generator){
  futures.push_back(pool.submit(worker(std::move(generator))));
}

R_F qr_parallel::compute(){
  /* gather results */
  bool is_first = true;
  arma::mat R_stack;
  arma::mat F_stack;
  arma::mat dev;

  unsigned int num_blocks = futures.size();
  for(unsigned int i = 0; i < num_blocks; ++i){
    R_F R_Fs_i = futures[i].get();
    if(is_first){
      R_stack = R_Fs_i.R_rev_piv();
      F_stack = std::move(R_Fs_i.F);
      dev = R_Fs_i.dev;
      is_first = false;
      continue;
    }

    R_stack = arma::join_cols(R_stack, R_Fs_i.R_rev_piv());
    F_stack = arma::join_cols(F_stack, std::move(R_Fs_i.F));
    dev += R_Fs_i.dev;
  }

  /* make new QR decomp and compute new F */
  QR_factorization qr(R_stack);
  arma::mat F = qr.qy(F_stack, true).rows(0, R_stack.n_cols - 1);

  return { qr.R(), qr.pivot(), std::move(F), dev };
}
