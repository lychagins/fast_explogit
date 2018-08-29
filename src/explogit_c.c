#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include "explogit.h"

void dgemv(double *X, size_t nX, double *b, size_t nb, double *Xb)
{
	size_t i, j;
	double u;
	
	for (i=0; i<nX; i++) {
		
		u = 0.0;
		for (j=0; j<nb; j++){
			u += b[j]*(*X++);
		}
		
		Xb[i] = u;
	}
}

void dgemv_t(double *X, size_t nX, double *b, size_t nb, double *Xb)
{
	size_t i, j;
	double u;
	
	for (j=0; j<nb; j++) {
		
		u = b[j];
		for (i=0; i<nX; i++){
			Xb[i] += u*(*X++);
		}
		
	}
}

double explogit(double *beta, size_t num_covariates, size_t num_students,\
	double *X, uint16_t *nskipped, uint16_t *nlisted, double *grad)
{
	/* Pointer to the head of X */
	double *X_first;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;

	/* Various dimensions of data */
	size_t num_choices, cssize, csmax;
	
	/* Pointers for the vector of mean values, their exponentials */
	double *u, *u_first, *v;

	/* Workspace for computing components of the gradient and the logit shares */
	double *numer, denom, expu, xb, x;
	
	/* Generic indices */
	size_t i, j, l, l_last;
	
	/*-------------------------------------------------------------------------
	 * Read and transform parameters
	 *-----------------------------------------------------------------------*/
	
	num_choices = 0;
	csmax = 0;
	for (i=0; i<num_students; i++) {
		cssize = nlisted[i] + nskipped[i];
		num_choices += cssize;
		if(csmax < cssize) {
			csmax = cssize;
		}
	}
	
	/*-------------------------------------------------------------------------
	 * Allocate workspace here
	 *-----------------------------------------------------------------------*/	
	/* To compute the gradient of loglikelihood, we need to know unconditional 
	 * probabilities. Therefore, we have to store some of the gradient's components */
	u_first = (double *)malloc(num_choices*sizeof(double));
	v = (double *)malloc(csmax*sizeof(double));
	numer = (double *)malloc(num_covariates*sizeof(double));
	
	u = u_first;
	dgemv(X, num_choices, beta, num_covariates, u);
	
	loglik = 0.0;
	memset(grad, 0, num_covariates*sizeof(double));
	
	/* Index for the "short stack" array */
	l = 0;

	/* Loop over students */
	for (i=0; i<num_students; i++){
		
		/* Initialize accumulators */
		denom = 0.0;
		memset(numer, 0, num_covariates*sizeof(double));

		/* Accumulate logit denominator and gradient's numerator over skipped choices */
		l_last = nskipped[i];
		for (l=0; l<l_last; l++){
			v[l] = exp(u[l]);
			denom += v[l];
		}
		u += l_last;
		
		dgemv_t(X, num_covariates, v, l_last, numer);
		X += num_covariates*l_last;
		
		/* Go over the preference list in reverse order */
		l_last = nlisted[i];
		for (l=0; l<l_last; l++) {
			xb = *u++;
			expu = exp(xb);
			denom += expu;
			for (j=0; j<num_covariates; j++) {
				x = *X++;
				numer[j] += expu*x;
				grad[j] += x - numer[j]/denom;
			}
			loglik += xb - log(denom);
		}
		
	}
	free(numer);
	free(v);
	free(u_first);

	return loglik;

}