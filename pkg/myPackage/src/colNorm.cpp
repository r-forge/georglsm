#include <RcppGSL.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>

extern "C" SEXP colNorm(SEXP sM) {
try {
RcppGSL::matrix<double> M = sM;
int k = M.ncol();
Rcpp::NumericVector n(k);
// create gsl data structures from SEXP
// to store results
for (int j = 0; j < k; j++) {
RcppGSL::vector_view<double> colview = gsl_matrix_column (M, j);
n[j] = gsl_blas_dnrm2(colview);
}
M.free() ;
return n;
// return vector
} catch( std::exception &ex ) {
forward_exception_to_r( ex );
} catch(...) {
::Rf_error( "c++ exception (unknown reason)" );
}
return R_NilValue; // -Wall
}

