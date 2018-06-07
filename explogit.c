#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <omp.h>
#include "explogit.h"


double explogit(double *raw_param, int num_types, int num_covariates, int num_students,\
	double *X, int *nskipped, int *nlisted, double *u_first, double *grad)
{

	/* Pointers for the matrix of covariates: first/last/current elements */
	double *X_first, *X_last;

	/* Same parameters in a more convenient shape: shares of types */
	double *w_t;
	/* Ditto: pointers to the current/first of the type utility parameters */
	double *beta, *beta_type;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;

	/* Various dimensions of data */
	size_t num_choices;

	/* Pointers for the vector of mean values: first/last/current/pref. list boundary */
	double *u_last, *u, *u_bound;

	/* Probability of observing the preference list; unconditional. First element and
	 *  a movable pointer. */
	double *pr, *pr_first;
	double *pr_type, *pr_type_first; /* Same probability, but conditional on type */
	
	/* Workspace for computing components of the gradient and the logit shares */
	double *dpr_type_db, *dpr_type_db_first, *numer, *dpr_mult, \
		denom, expu, xb, pr_cur, x, p_ratio;

	/* Generic indices */
	size_t i, j, k;

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

	/* Type-specific utility parameters */
	beta_type = raw_param;

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
	
	/* Initialize pointers */
	
	pr = pr_first;
	X_first = X;
	X_last = X_first + num_covariates*num_choices;
	pr_type = pr_type_first;
	dpr_type_db = dpr_type_db_first;
	
	#pragma omp parallel for \
		private(u, u_bound, u_last, w_cur, xb, beta, denom, logpr_type, numer, dpr_mult, expu, i, j, pr_cur, k) \
		shared(X_last, X_first, num_choices, num_covariates, u_first) \
		firstprivate(nskipped, nlisted, X, pr, beta_type, dpr_type_db)
	for (i=0; i<num_types; i++) {
		
		k = 0;
		beta_type += i*num_covariates;
		dpr_type_db += i*num_students*num_covariates;
		
		numer = (double *)malloc(num_covariates*sizeof(double));
		dpr_mult = (double *)malloc(num_covariates*sizeof(double));

		w_cur = w_t[i];

		/* Reset pointers */
		u = u_first + i*num_choices;
		u_last = u + num_choices;
		X = X_first;
		u_bound = u;

		/* Loop over students */
		while (u<u_last) {
			
			u_bound += *nskipped++;
			
			/* Initialize accumulators */
			denom = 0;
			logpr_type = 1;
			for (j=0; j<num_covariates; j++) {
				numer[j] = 0;
				dpr_mult[j] = 0;
			}
			
			
			/* Accumulate logit denominator and gradient's numerator over skipped choices */
			while (u < u_bound) {
				expu = *u++;
				for (j=0; j<num_covariates; j++) {
					numer[j] += expu*(*X++);
				}
				denom += expu;
			}
			
			/* Go over the preference list in reverse order */
			u_bound += *nlisted++;
			while (u < u_bound) {
				expu = *u++;
				denom += expu;
				for (j=0; j<num_covariates; j++) {
					x = *X++;
					numer[j] += expu*x;
					dpr_mult[j] += x - numer[j]/denom;
				}
				logpr_type *= expu/denom;
			}
			
			
			for (j=0; j<num_covariates; j++) {
				*dpr_type_db++ = dpr_mult[j];
			}
			
			pr_cur = logpr_type;
			*pr_type++ = pr_cur;
			#pragma omp atomic
			pr[k] += pr_cur*w_cur;
			k++;
		}
		
		free(dpr_mult);
		free(numer);
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