#include <Rconfig.h>
#include <R_ext/BLAS.h>
#include <R_ext/Lapack.h>
#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <numeric>
#include "mvphi.h"

#ifndef  USE_FC_LEN_T
# define USE_FC_LEN_T
#endif
#ifndef FCONE
# define FCONE
#endif


using namespace Rcpp;
using namespace std;

//' Univariate variable reordering, described in Genz and Bretz (2009)
//' If failed due to PD singularity, the unfinished order will be returned 
//' and a warning will be issued
//'
//' @param a lower integration limits
//' @param b upper integration limits
//' @param sigma covariance matrix
//' @return the new order
//' @export
// [[Rcpp::export]]
IntegerVector univar_order(NumericVector a, NumericVector b, 
    NumericMatrix sigma)
{
    int n = a.size();
	NumericVector a_cp = clone(a);
	NumericVector b_cp = clone(b);
	NumericMatrix sigma_cp = clone(sigma);
    double *a_ptr = a_cp.begin();
    double *b_ptr = b_cp.begin();
    double *sigma_ptr = sigma_cp.begin();
    double *sd = new double[n];
    double *aScaled = new double[n];  // work+n;
	double *bScaled = new double[n]; 
	double *aScaledCDF = new double[n];
	double *bScaledCDF = new double[n];
	double *pAll = new double[n];
    double *y = new double[n];
    IntegerVector order(n); 
    for(int i = 0; i < n; i++) order[i] = i;
    int *oldIdx = order.begin();

    for(int i = 0; i < n; i++)
	{
		for(int j = i; j < n; j++) 
            sd[j-i] = sqrt(sigma_ptr[n*j+j]);
		copy_n(a_ptr+i, n-i, aScaled);
		copy_n(b_ptr+i, n-i, bScaled);
		transform(aScaled, aScaled+n-i, sd, aScaled, 
            [](double aCoef, double sdCoef){return aCoef/sdCoef;});
		transform(bScaled, bScaled+n-i, sd, bScaled, 
            [](double bCoef, double sdCoef){return bCoef/sdCoef;});
		lc_vdCdfNorm(n-i, aScaled, aScaledCDF);
		lc_vdCdfNorm(n-i, bScaled, bScaledCDF);
		transform(bScaledCDF, bScaledCDF+n-i, aScaledCDF, pAll, [](double p1,
			double p2){return p1 - p2;});
		double *minCoef = min_element(pAll, pAll+n-i);
		int minCoefIdx = minCoef - pAll + i;
		// swap sigma_ptr, a_ptr, b_ptr
		iter_swap(oldIdx+i, oldIdx+minCoefIdx);
		iter_swap(sigma_ptr+i+n*i, sigma_ptr+minCoefIdx+n*minCoefIdx);
		for(int j = 0; j < i; j++)
			iter_swap(sigma_ptr+i+n*j, sigma_ptr+minCoefIdx+n*j);
		for(int j = i+1; j < minCoefIdx; j++) 
			iter_swap(sigma_ptr+j+n*i, sigma_ptr+minCoefIdx+n*j);
		for(int j = minCoefIdx+1; j < n; j++)
			iter_swap(sigma_ptr+j+n*i, sigma_ptr+j+n*minCoefIdx);
		iter_swap(a_ptr+i, a_ptr+minCoefIdx);
		iter_swap(b_ptr+i, b_ptr+minCoefIdx);
		// Cholesky
		if(*(sigma_ptr+i+i*n) <= 0) {
            warning("Singularity occurred at univariate reordering\n");
            break;
        }	
		*(sigma_ptr+i+i*n) = sqrt(*(sigma_ptr+i+i*n));
		if(i < n-1)
		{
			for_each(sigma_ptr+i+1+i*n, sigma_ptr+n+i*n, [&quo = *(sigma_ptr+i+i*n)](
				double &Bcoef){Bcoef /= quo;});
			int nrowLeft = n-i-1;
			double alpha = -1.0;
			int step = 1;
			F77_CALL(dsyr)("L", &nrowLeft, &alpha, sigma_ptr+i+1+i*n, &step, 
				sigma_ptr+i+1+(i+1)*n, &n FCONE);
		} // i < n-1
		// compute y
		double aScaledCoef = a_ptr[i] / *(sigma_ptr+i+i*n);
		double bScaledCoef = b_ptr[i] / *(sigma_ptr+i+i*n);
		double aScaledCoefCDF;
		double bScaledCoefCDF;
		lc_vdCdfNorm(1, &aScaledCoef, &aScaledCoefCDF);
		lc_vdCdfNorm(1, &bScaledCoef, &bScaledCoefCDF);
		double pCoef = bScaledCoefCDF - aScaledCoefCDF;
		y[i] = ((exp(-aScaledCoef*aScaledCoef/2.0) - exp(-bScaledCoef*
			bScaledCoef/2.0)) / sqrt(2.0*M_PI)) / pCoef;
		// update a_ptr, b_ptr
		if(i < n-1)
		{
			transform(a_ptr+i+1, a_ptr+n, sigma_ptr+i+1+n*i, a_ptr+i+1, [&yCoef = y[i]](
				double aCoef, double BCoef){return aCoef - BCoef*
				yCoef;});
			transform(b_ptr+i+1, b_ptr+n, sigma_ptr+i+1+n*i, b_ptr+i+1, [&yCoef = y[i]](
				double bCoef, double BCoef){return bCoef - BCoef*
				yCoef;});
		}
	}

    delete[] sd;
    delete[] aScaled;
    delete[] bScaled;
    delete[] aScaledCDF;
    delete[] bScaledCDF;
    delete[] pAll;
    delete[] y;

	for(int i = 0; i < n; i++) order[i]++;
    return order;
}