#include <R.h>
#include <Rdefines.h>
#include <Rmath.h>
#include <math.h>

#ifndef __APPLE__
#include "omp.h"
#endif

#include "stat_rank.h"
/*#include "wmw_test.h"*/

#define MIN(x,y) ((x) > (y) ? (y) : (x))
#define NROW(x) INTEGER(GET_DIM((x)))[0]
#define NCOL(x) INTEGER(GET_DIM((x)))[1]
#define ABSLOG(x) fabs(log10( (x) ))

// [[Rcpp::interfaces(r, cpp)]]

/*#include <Rcpp.h>

using namespace Rcpp;*/

/*	greater=0,
	less=1,
	twoSided=2,
	U=3,
	abslog10greater=4,
	log10less=5,
	abslog10twoSided=6,
	Q=7
*/

/*
 * void pnorm_both(double x, double *cum, double *ccum, int i_tail, int log_p)

 * i_tail in {0,1,2} means: "lower", "upper", or "both" :
 * if(lower) return  *cum := P[X <= x]
 * if(upper) return *ccum := P[X >  x] = 1 - P[X <= x]
 */
