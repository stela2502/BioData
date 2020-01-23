// [[Rcpp::depends(RcppEigen)]]
#include <RcppEigen.h>
using namespace Rcpp;
typedef Eigen::MappedSparseMatrix<double> MSpMat;

// [[Rcpp::export]]
void mat_loop1 (MSpMat X, int j, ){
    Rcout << "Standard looping over a sparse matrix" << std::endl;
    for (int i=0; i<X.rows(); ++i){
        Rcout << " i,j=" << i << "," << j << " value=" << X.coeff(i,j) << std::endl;
    }
}
