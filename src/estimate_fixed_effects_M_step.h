#ifndef DDHAZARD_EST_M_STEP_H
#define DDHAZARD_EST_M_STEP_H
#include "family.h"
#include "problem_data.h"
#include "arma_n_rcpp.h"

// Method to estimate fixed effects like in biglm::bigglm
void estimate_fixed_effects_M_step(
    ddhazard_data * const, arma::uword, family_base&);
#endif
