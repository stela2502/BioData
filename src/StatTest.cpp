// [[Rcpp::depends(RcppEigen)]]
// [[Rcpp::depends(RcppProgress)]]
#include <RcppEigen.h>
/*#include <progress.hpp> */
#include <math.h>
#include <stat_rank.h>
using namespace Rcpp;
#include <vector>
#include <stdexcept>
typedef Eigen::MappedSparseMatrix<double> MSpMat;
#include <Rcpp/Dimension.h>

#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>

#define MIN(x,y) ((x) > (y) ? (y) : (x))
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#define ABSLOG(x) fabs(log10( (x) ))




// [[Rcpp::export]]
double logFC ( std::vector<double> A, std::vector<double> B ) {
	double res = 0.0;
	double Asum = A[0] * 1.0;
	double Bsum = B[0] * 1.0;
	for ( int i=1; i<A.size(); i++ ){
		Asum = log( exp( Asum - A[i]) + 1.0)  + A[i];
	}
	for ( int i=1; i<B.size(); i++ ){
		Bsum = log( exp( Bsum - B[i]) + 1.0 ) + B[i];
	}
	/*Rcout << "Int values A; a size; B; b size:" << Asum <<";"<< A.size()<<";"<< Bsum <<";"<< B.size() << std::endl;*/
	return (Asum - log(A.size()))-(Bsum - log(B.size())) ;
}

std::vector<int> minusOne ( std::vector<int>  X ){
	for ( int i = 0; i < X.size(); i ++) {
		X[i] --;
	}
	return X;
}

std::vector<int> plusOne ( std::vector<int>  X ){
	for ( int i = 0; i < X.size(); i ++) {
		X[i] ++;
	}
	return X;
}

/* direct copy from BioQC/src/wmw_test.c */
/*	greater=0,
	less=1,
	twoSided=2,
	U=3,
	abslog10greater=4,
	log10less=5,
	abslog10twoSided=6,
	Q=7
 */
double wmw_test_stat(double rankSum, int nInds, int nTotal, double tieCoef, int type) {

	double uStat, mu, sigma2, zval, pgt, plt;
	double res;
	int nBg = nTotal-nInds;

	uStat = nInds*nBg+nInds*(nInds+1.0)*0.5-rankSum;

	if(type == 3) {
		res = uStat;
	} else {
		mu = (double)nInds*nBg*0.5; // NOT mu=n1*n2*0.5
		sigma2 = nInds*nBg*(nTotal+1.0)/12.0*tieCoef; //NOT sigma2 = n1*n2*(n+1)/12*tieCoef

		if(type == 0 || type == 4) { /* greater */
			zval = (uStat+0.5-mu)/sqrt(sigma2); // z lower tail
			R::pnorm_both(zval, &pgt, &plt, 0, 0);
			res = type==0 ? pgt : ABSLOG(pgt);
		} else if (type == 1 || type == 5) { /* less */
			zval = (uStat-0.5-mu)/sqrt(sigma2); // z higher tail
			R::pnorm_both(zval, &pgt, &plt, 1, 0);
			res = type==1 ? plt : log10(plt);
		} else if (type == 2 || type == 6 || type == 7) { /* two sided*/
			zval = (uStat-mu-(uStat>mu ? 0.5 : -0.5))/sqrt(sigma2);
			R::pnorm_both(zval, &pgt, &plt, 2, 0);
			res = mu==0.0 ? 1.0 : 2.0*MIN(pgt, plt);
			if(type == 4) {
				res = ABSLOG(res);
			} else if (type == 7) {
				res = pgt<=plt ? ABSLOG(res) : -ABSLOG(res);
			}
		} else {
			/* error("Unrecognized type %d. Should not happen\nPossible only int values  0=greater, 1=less, 2=twoSided, 3=U, 4=abslog10greater, 5=log10less, 6=abslog10twoSided, 7=Q",
            type); */
		}
	}
	return(res);
}

//’ StatTest runs wilcox test on the columns of the sparse matrix
//’
//’ @param X the sparse matrix
//’ @param interest row IDs for the group of interest
//’ @param backgound row IDS for the background
//’ @param logFCcut data is meant to be log() transformed and only columns passing a logFCcut of (default 1) are tested
//’ @param display_progress unused
//’ @return a matrix with tested column ids, logFC and p.value
// [[Rcpp::export]]
extern SEXP StatTest (Eigen::MappedSparseMatrix<double> X, std::vector<int> interest,
		std::vector<int> backgound, double logFCcut = 1.0, bool display_progress=true ){

	Rcout << "Standard looping over a sparse matrix initializing" << std::endl;

	std::vector<double> logFCpass(X.cols(), 0.0);
	std::vector<double> A(interest.size(), 0.0);
	std::vector<double> B(backgound.size(), 0.0 );
	int pass = 0;
	std::vector<int> itA = minusOne( interest );
	std::vector<int> itB = minusOne( backgound );

	Rcout << "Standard looping over a sparse matrix calculating logFC" << std::endl;
	for ( int c_=0; c_ < X.cols(); ++c_ ){
		for (int i = 0; i< itA.size(); i++ ) {
			if ( itA.at(i) < 0 || itA.at(i) >= X.rows() ) {
				throw std::invalid_argument( "test out of bounds" );
			}
			A[i] = X.coeff(itA.at(i),c_);
		}
		for ( int i = 0; i< itB.size(); i++ ) {
			if (itB.at(i) < 0 || itB.at(i) >= X.rows() ) {
				throw std::invalid_argument( "itB out of bounds" );
			}
			B.at(i) = X.coeff(itB[i],c_);
		}
		logFCpass[c_] = logFC( A, B );
		if ( logFCpass.at(c_) > logFCcut ) {
			pass++;
		}
	}

	Rcout << "Standard looping over a sparse matrix calculating wilcox test" << std::endl;
	/* allocate a result 'matrix' */
	NumericMatrix res(pass, 4);
	int n = X.rows();
	double *total = new double[ itA.size() + itB.size() ];
	int id = 0;
	DRankList list;
	int j;
	int nInd;
	double tie;
	double indRankSum;

	for ( int c_=0; c_ < X.cols(); c_++ ){
		if ( logFCpass[c_] > logFCcut ) {

			/*Test stats copied from the BioOC package */
			j = 0;
			for (unsigned int i = 0; i< itA.size(); i++ ) {
				total[j++] = X.coeff(itA.at(i) ,c_);
			}
			for (unsigned int i = 0; i< itB.size(); i++ ) {
				total[j++] = X.coeff(itB.at(i) ,c_);
			}
			n = j;
			list = createDRankList(total, n);
			prepareDRankList(list);

			tie = tieCoef(list);

			nInd=itA.size();

			indRankSum = 0.0;
			for(j=0; j<nInd; ++j) {
				if(!(itA.at(j)>=0 && itA.at(j)<=n-1))
					::Rf_error("Index out of range: gene set %d, gene %d\n", c_+1, j+1);
				//if ( total[itA.at(j)] > 0 ){
					indRankSum += list->list[itA.at(j)]->rank;
				//}
			}
			/* store the results */
			res(id,0) = c_ + 1;
			res(id,1) = logFCpass.at(c_);
			res(id,2) = indRankSum;
			/* store the higher p value as we do drop all lower anyhow. */
			res(id,3) = wmw_test_stat(indRankSum, nInd, n, tie, 0);
			id ++;
		}
	}
	colnames(res) = CharacterVector::create("colID", "logFC", "rank.sum", "p.value");
	Rcout << "n return values: " << pass <<std::endl;
	return res;
}


