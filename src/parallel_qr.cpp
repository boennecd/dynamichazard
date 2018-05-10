#include "parallel_qr.h"
#include "thread_pool.h"
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
  ptr_vec generators, const signed int max_threads):
  generators(std::move(generators)), max_threads(max_threads) {}

R_F qr_parallel::compute(){
  std::vector<R_F> R_Fs;
  {
    /* compute each QR */
    signed int num_blocks = generators.size();
    std::vector<std::future<R_F> > futures(num_blocks);
    thread_pool pool(std::min(num_blocks, max_threads));

    signed int i = 0;
    while(!generators.empty()){
      futures[i] = pool.submit(worker(std::move(generators.back())));
      generators.pop_back();
      ++i;

    }

    /* gather results */
    for(i = 0; i < num_blocks; ++i)
      R_Fs.push_back(futures[i].get());
  }

  /* combine results */
  R_F *R_Fs_i = &R_Fs.back();
  arma::mat R_stack = R_Fs_i->R_rev_piv();
  arma::vec F_stack = std::move(R_Fs_i->F);
  double  dev = R_Fs_i->dev;
  R_Fs.pop_back();

  while(!R_Fs.empty()){
    R_Fs_i = &R_Fs.back();
    R_stack = arma::join_cols(R_stack, R_Fs_i->R_rev_piv());
    F_stack = arma::join_cols(F_stack, std::move(R_Fs_i->F));
    dev += R_Fs_i->dev;

    R_Fs.pop_back();
  }

  /* make new QR decomp and compute new F */
  QR_factorization qr(R_stack);
  arma::mat F = qr.qy(F_stack, true).rows(0, R_stack.n_cols - 1);

  return { qr.R(), qr.pivot(), std::move(F), dev };
}
