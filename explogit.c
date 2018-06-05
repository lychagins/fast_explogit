#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "explogit.h"

double explogit(double *raw_param, int num_types, int num_covariates, int num_students,\
	double *X, int *nskipped, int *nlisted, double *grad)
{
	/* raw_param -- raw parameter vector as passed by fminunc */
	

	/* Pointers for the matrix of covariates: first/last/current elements */
	double *X_first, *X_last;
	
	/* Same parameters in a more convenient shape: shares of types */
	double *w_t;
	/* Ditto: pointers to the current/first of the type utility parameters */
	double *beta, *beta_first, *beta_cur, *beta_last;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;

	/* Pointers for the vector of mean values: first/last/current */
	double *u, *u_first, *u_last;

	/* Probability of observing the preference list; unconditional */
	double pr;
	
	/* Workspace for computing components of the gradient and the logit shares */
	double *numer, *numer_first, *numer_last, *dpr_mult, *dpr_mult_first, *dpr_mult_cur,\
		denom, expu, x, xb, p_ratio;
		
	/* Generic indices */
	size_t i, j, t, cmax;
	
	double *dldw, *dldb, *dldb_first;

	/* In case we want to accumulate individual likelihood in the log space */
	double *logpr_type, *logpr_type_first;
	
	
	/*-------------------------------------------------------------------------
	 * Read and transform parameters
	 *-----------------------------------------------------------------------*/	

	 w_t = (double *)malloc(num_types*sizeof(double));
	w_t[0] = 1;
	denom = 1;
	for (t=1; t<num_types; t++) {
		expu = exp(*raw_param++);
		w_t[t] = expu;
		denom += expu;
	}

	for (t=0; t<num_types; t++) {
		w_t[t] /= denom;
	}

	/* Type-specific utility parameters */
	beta_first = raw_param;

	loglik = 0;

	/*-------------------------------------------------------------------------
	 * Allocate workspace here
	 *-----------------------------------------------------------------------*/
	
	cmax = nskipped[0] + nlisted[0];
	for (i=1; i<num_students; i++) {
		t = nskipped[i] + nlisted[i];
		if (cmax < t) {
			cmax = t;
		}
	}
	
	u_first = (double *)malloc(cmax*sizeof(double));
	logpr_type_first = (double *)malloc(num_types*sizeof(double));
	/* To compute the gradient of loglikelihood, we need to know unconditional 
	 * probabilities. Therefore, we have to store some of the gradient's components */
	dpr_mult_first = (double *)malloc(num_covariates*num_types*sizeof(double));
	/* One component of the gradient is a sum of fractions. This is a workspace for 
	 * computing the numerator for these fractions. */
	numer_first = (double *)malloc(num_covariates*sizeof(double));
		
	/* This is where we accumulate the gradient w.r.t. the type weight parameters */
	dldw = (double *)calloc(num_types, sizeof(double));

	/* This is where we accumulate the gradient w.r.t. to the type-specific utility parameters */
	dldb_first = (double *)calloc(num_types*num_covariates, sizeof(double));
	
	/* Initialize pointers */
	numer_last = numer_first + num_covariates;	
	
	for (i=0; i<num_students; i++) {
		
		/* Student i's frame in the matrix of covariates */
		X_first = X;
		X_last = X_first + (*nskipped + *nlisted)*num_covariates;
		
		pr = 0.0;
		
		/* Reset pointers to the first type's frame */
		beta_cur = beta_first;
		dpr_mult_cur = dpr_mult_first;
		logpr_type = logpr_type_first;

		for (t=0; t<num_types; t++) {

			/* Inner product for every feasible choice, store in u[] */
			X = X_first;
			u = u_first;
			beta_last = beta_cur + num_covariates;
			while (X<X_last) {
				xb = 0.0;
				beta = beta_cur;
				while (beta < beta_last) {
					xb += (*beta++)*(*X++);
				}
				*u++ = xb;
			}

			/* Initialize accumulators */
			denom = 0.0;
			numer = numer_first;
			dpr_mult = dpr_mult_cur;
			while (numer<numer_last) {
				*numer++ = 0.0;
				*dpr_mult++ = 0.0;
			}
			
			/* Skipped choices */
			u = u_first;
			u_last = u_first + *nskipped;
			X = X_first;
			while (u < u_last) {
				expu = exp(*u++);
				denom += expu;
				numer = numer_first;
				while (numer<numer_last) {
					*numer += (*X++)*expu;
					numer++;
				}
			}

			/* Listed choices */
			u_last += *nlisted;
			*logpr_type = 0.0;
			while (u < u_last) {
				xb = *u++;
				expu = exp(xb);
				denom += expu;
				numer = numer_first;
				dpr_mult = dpr_mult_cur;
				while (numer<numer_last) {
					x = *X++;
					*numer += x*expu;
					*dpr_mult += x - (*numer++)/denom;
					dpr_mult++;
				}
				*logpr_type += xb - log(denom);
			}
			
			pr += w_t[t]*exp(*logpr_type++);
			
			/* Move pointers to the next type's frame */
			beta_cur += num_covariates;
			dpr_mult_cur += num_covariates;
			
		}
		
		nskipped++;
		nlisted++;
		
		/* We've got the likelihood function */
		loglik -= log(pr);

		dpr_mult = dpr_mult_first;
		dldb = dldb_first;
		logpr_type = logpr_type_first;
		for (t=0; t<num_types; t++) {
			p_ratio = w_t[t]*exp(*logpr_type++)/pr;
			dldw[t] += p_ratio - w_t[t];
			for (j=0; j<num_covariates; j++) {
				*dldb += (*dpr_mult++)*p_ratio;
				dldb++;
			}
		}
				
	}
	
    /*-------------------------------------------------------------------------
	 * Prepare the results for output
	 *-----------------------------------------------------------------------*/
	for (t=1; t<num_types; t++) {
		*grad++ = -dldw[t];
	}
	
	dldb = dldb_first;
	for (j=0; j<num_covariates*num_types; j++) {
		*grad++ = -(*dldb++);
	}

	/* Clean up */
	free(w_t);	
	free(u_first);
	free(logpr_type_first);
	free(dpr_mult_first);
	free(numer_first);
	free(dldw);
	free(dldb_first);

	return loglik;

}
