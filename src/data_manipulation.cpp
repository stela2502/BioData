#include <RcppEigen.h>
#include <progress.hpp>
#include <cmath>
#include <unordered_map>
#include <fstream>
#include <string>

using namespace Rcpp;
// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]



// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISampling(Eigen::SparseMatrix<double> data, int sample_val, bool upsample = false, bool display_progress=true){
    Progress p(data.outerSize(), display_progress);
    Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.cols());
    for (int k=0; k < data.outerSize(); ++k){
      p.increment();
      for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
        double entry = it.value();
        if( (upsample) || (colSums[k] > sample_val)){
          entry = entry * double(sample_val) / colSums[k];
          if (fmod(entry, 1) != 0){
            double rn = R::runif(0,1);
            if(fmod(entry, 1) <= rn){
              it.valueRef() = floor(entry);
            }
            else{
              it.valueRef() = ceil(entry);
            }
          }
          else{
            it.valueRef() = entry;
          }
        }
      }
    }
  return(data);
}

// [[Rcpp::export]]
Eigen::SparseMatrix<double> RunUMISamplingPerCell(Eigen::SparseMatrix<double> data, NumericVector sample_val, bool upsample = false, bool display_progress=true){
  Progress p(data.outerSize(), display_progress);
  Eigen::VectorXd colSums = data.transpose() * Eigen::VectorXd::Ones(data.cols());
  return(data);
  for (int k=0; k < data.outerSize(); ++k){
    p.increment();
    for (Eigen::SparseMatrix<double>::InnerIterator it(data, k); it; ++it){
      double entry = it.value();
      if( (upsample) || (colSums[k] > sample_val[k])){
        entry = entry * double(sample_val[k]) / colSums[k];
        if (fmod(entry, 1) != 0){
          double rn = R::runif(0,1);
          if(fmod(entry, 1) <= rn){
            it.valueRef() = floor(entry);
          }
          else{
            it.valueRef() = ceil(entry);
          }
        }
        else{
          it.valueRef() = entry;
        }
      }
    }
  }
  return(data);
}

