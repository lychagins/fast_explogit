#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "explogit.h"


double explogit(double *raw_param, int num_types, int num_covariates, int num_students,\
	double *X, int *nskipped, int *nlisted, double *grad)
{

	/* Pointers for the matrix of covariates: first/last/current elements */
	double *X_first;

	/* Same parameters in a more convenient shape: shares of types */
	double *w_t;
	/* Ditto: pointers to the current/first of the type utility parameters */
	double *beta;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;

	/* Various dimensions of data */
	size_t num_choices;
	
	/* Pointers for the vector of mean values: first/last/current/pref. list boundary */
	double *u;

	/* Probability of observing the preference list; unconditional. First element and
	 *  a movable pointer. */
	double *pr, *pr_first;
	double *pr_type, *pr_type_first; /* Same probability, but conditional on type */
	
	/* Workspace for computing components of the gradient and the logit shares */
	double *dpr_type_db, *dpr_type_db_first, *numer, *dpr_mult, \
		denom, expu, xb, x, p_ratio;
	
	/* Generic indices */
	size_t i, j, k, l, l_last, l_first;

	double *dldw, *dldb, *dldb_first;

	/* Workspace for fast access to type weights, mean utility (perhaps, compiler 
	 * optimizations make using these variables) */
	double w_cur;
	/* In case we want to accumulate individual likelihood in the log space */
	double logpr_type;
	
	
	/*-------------------------------------------------------------------------
	 * Read and transform parameters
	 *-----------------------------------------------------------------------*/
	
	num_choices = 0;
	for (i=0; i<num_students; i++) {
		num_choices += nlisted[i] + nskipped[i];
	}
	
	w_t = (double *)malloc(num_types*sizeof(double));
	w_t[0] = 1;
	denom = 1;
	for (i=1; i<num_types; i++) {
		expu = exp(*raw_param++);
		w_t[i] = expu;
		denom += expu;
	}

	for (i=0; i<num_types; i++) {
		w_t[i] /= denom;
	}

	/*-------------------------------------------------------------------------
	 * Allocate workspace here
	 *-----------------------------------------------------------------------*/	
	/* To compute the gradient of loglikelihood, we need to know unconditional 
	 * probabilities. Therefore, we have to store some of the gradient's components */
	pr_first = (double *)calloc(num_students, sizeof(double));
	pr_type_first = (double *)malloc(num_students*num_types*sizeof(double));
	
	/* Derivatives of the conditional choice probability, by student and type */
	dpr_type_db_first = (double *)calloc(num_students*num_types*num_covariates,\
		sizeof(double));
	
	/* This is where we accumulate the gradient w.r.t. the type weight parameters */
	dldw = (double *)calloc(num_types, sizeof(double));

	/* This is where we accumulate the gradient w.r.t. to the type-specific utility parameters */
	dldb_first = (double *)calloc(num_types*num_covariates, sizeof(double));
	
	X_first = X;
	
	#pragma omp parallel for \
		private(u, xb, beta, denom, logpr_type, l, l_first, l_last, \
			numer, dpr_mult, expu, i, j, k, pr_type, X, x, dpr_type_db) \
		shared(X_first, pr_type_first, num_choices, num_covariates,\
			num_types, num_students, nlisted, nskipped, dpr_type_db_first, raw_param) \
		default(none)
	for (i=0; i<num_types; i++) {
		
		/* Initialize pointers for type i */
		pr_type = pr_type_first + i*num_students;
		beta = raw_param + i*num_covariates;
		dpr_type_db = dpr_type_db_first + i*num_students*num_covariates;
		
		u = (double *)malloc(num_choices*sizeof(double));
		numer = (double *)malloc(num_covariates*sizeof(double));
		dpr_mult = (double *)malloc(num_covariates*sizeof(double));
				
		/* Dot-product; naive algorithm */
		X = X_first;
		for (l=0; l<num_choices; l++) {
			
			xb = 0;
			for (j=0; j<num_covariates; j++) {
				xb += (*X++)*beta[j];
			}			
			u[l] = xb;
		}	

		/* Reset X pointer */
		X = X_first;
		
		/* Index for the "short stack" array */
		l = 0;

		/* Loop over students */
		for (k=0; k<num_students; k++){
			
			/* Initialize accumulators */
			denom = 0.0;
			logpr_type = 0.0;
			for (j=0; j<num_covariates; j++) {
				numer[j] = 0.0;
				dpr_mult[j] = 0.0;
			}


			/* Accumulate logit denominator and gradient's numerator over skipped choices */

			l_last = l + nskipped[k];
			l_first = l;
			for (l=l_first; l<l_last; l++){
				expu = exp(u[l]);
				for (j=0; j<num_covariates; j++) {
					numer[j] += expu*(*X++);
				}
				denom += expu;
			}

			/* Go over the preference list in reverse order */			
			
			l_last = l + nlisted[k];
			l_first = l;
			for (l=l_first; l<l_last; l++) {
				xb = u[l];
				expu = exp(xb);
				denom += expu;
				for (j=0; j<num_covariates; j++) {
					x = *X++;
					numer[j] += expu*x;
					dpr_mult[j] += x - numer[j]/denom;
				}
				logpr_type += xb - log(denom);
			}
			
			
			for (j=0; j<num_covariates; j++) {
				*dpr_type_db++ = dpr_mult[j];
			}
			
			pr_type[k] = exp(logpr_type);
		}
				
		free(dpr_mult);
		free(numer);
		free(u);
	}
	
	pr = pr_first;
	pr_type = pr_type_first;
	for (i=0; i<num_types; i++) {
		w_cur = w_t[i];
		for (k=0; k<num_students; k++) {
			pr[k] += w_cur*(*pr_type++);
		}
	}

	/*-------------------------------------------------------------------------
	 * Accumulate the gradients
	 *-----------------------------------------------------------------------*/
	/* Reset the pointers */
	pr_type = pr_type_first;
	dpr_type_db = dpr_type_db_first;
	dldb = dldb_first;
	
	for (i=0; i<num_types; i++) {

		pr = pr_first;
		w_cur = w_t[i];
		
		for (k=0; k<num_students; k++) {
			
			/* Compute conditional-to-marginal probability ratio. This is a
			 * common multiplier in all gradient terms */
			p_ratio = w_cur*(*pr_type++)/(*pr++);
			
			/* Accumulate parts of the gradient responsible for type mix */
			dldw[i] += p_ratio - w_cur;
			
			for (j=0; j<num_covariates; j++) {
				
				dldb[j] += p_ratio*(*dpr_type_db++);
				
			}
			
		}
		
		dldb += num_covariates;
		
	}
	
	
    /*-------------------------------------------------------------------------
	 * Prepare the results for output
	 *-----------------------------------------------------------------------*/
	/* Save the negative loglikelihood */
	pr = pr_first;
	loglik = 0;
	for (i=0; i<num_students; i++) {
		loglik -= log(pr[i]);
	}
	
	/* Save the gradient */
	
	for (i=1; i<num_types; i++) {
		*grad++ = -dldw[i];
	}
	dldb = dldb_first;
	for (j=0; j<num_types*num_covariates; j++) {
		*grad++ = -dldb[j];
	}
	
	free(w_t);
	free(pr_first);
	free(pr_type_first);
	free(dpr_type_db_first);
	free(dldw);
	free(dldb);

	return loglik;

}