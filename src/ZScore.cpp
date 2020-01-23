// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
/*#include <progress.hpp> */
#include <math.h>
using namespace Rcpp;
typedef Eigen::SparseVector<double> SpVec;
typedef SpVec::InnerIterator InIterVec;

//[[Rcpp::export]]
Eigen::SparseMatrix<double> ZScore (Eigen::SparseMatrix<double> data, bool display_progress=true){
	/* Progress p(data.outerSize(), display_progress); */
	data = data.transpose();
	for (int k=0; k < data.outerSize(); ++k){
		/*p.increment();*/
		double sum = 0.0;
		int c = 0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			c++;
			sum += it.value();
		}
		double mean = sum / c;
		sum = 0.0;
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			double entry = (it.value() - mean);
			sum += entry * entry;
			it.valueRef() = entry;
		}
		double sd = sqrt(sum/c);
		/*Rcout << k << " mean " << mean << " and sd " << sd << "with count "<< c<< std::endl;*/
		for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
			double entry = (it.value() / sd) +13.0 ;
			it.valueRef() = entry;
		}
	}
	data = data.transpose();
    return (data);
}
