#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "explogit.h"


double explogit(double *raw_param, int num_types, int num_covariates, int num_students,\
	double *X, int *nskipped, int *nlisted, double *grad)
{

	/* Pointers for the matrix of covariates: first/last/current elements */
	double *X_first, *X_last;

	/* Same parameters in a more convenient shape: shares of types */
	double *w_t;
	/* Ditto: pointers to the current/first of the type utility parameters */
	double *beta, *beta_type;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;

	/* Two arrays mapping the structure of the short stack and related pointers */
	int *nskipped_first, *nlisted_first;

	/* Various dimensions of data */
	size_t num_choices;

	/* Pointers for the vector of mean values: first/last/current/pref. list boundary */
	double *u_first, *u_last, *u, *u_bound;

	/* Probability of observing the preference list; unconditional. First element and
	 *  a movable pointer. */
	double *pr, *pr_first;
	double *pr_type, *pr_type_first; /* Same probability, but conditional on type */
	
	/* Workspace for computing components of the gradient and the logit shares */
	double *numer, *dpr_type_db, *dpr_type_db_first, *dpr_mult,\
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
	u_first = (double *)malloc(num_choices*sizeof(double));
	
	dpr_mult = (double *)malloc(num_covariates*sizeof(double));
	
	/* Derivatives of the conditional choice probability, by student and type */
	dpr_type_db_first = (double *)calloc(num_students*num_types*num_covariates,\
		sizeof(double));
	
	/* One component of the gradient is a sum of fractions. This is a workspace for 
	 * computing the numerator for these fractions. */
	numer = (double *)malloc(num_covariates*sizeof(double));
	
	/* This is where we accumulate the gradient w.r.t. the type weight parameters */
	dldw = (double *)calloc(num_types, sizeof(double));

	/* This is where we accumulate the gradient w.r.t. to the type-specific utility parameters */
	dldb_first = (double *)calloc(num_types*num_covariates, sizeof(double));
	
	/* Initialize pointers */
	u_last = u_first + num_choices;
	nskipped_first = nskipped;
	nlisted_first = nlisted;
	X_first = X;
	X_last = X_first + num_covariates*num_choices;
	pr_type = pr_type_first;
	dpr_type_db = dpr_type_db_first;

	for (i=0; i<num_types; i++) {
		
		/* Reset pointers */
		nskipped = nskipped_first;
		nlisted = nlisted_first;
		X = X_first;
		pr = pr_first;
		u = u_first;
		
		w_cur = w_t[i];
		
		/* Dot-product; naive algorithm */
		while (X<X_last) {
			
			xb = 0;
			beta = beta_type;
			for (j=0; j<num_covariates; j++) {
				xb += (*X++)*(*beta++);
			}			
			*u++ = xb;
		}
		/* Move beta vector pointer to the next type's position */
		beta_type += num_covariates;

		/* Reset pointers */
		u = u_first;
		X = X_first;
		u_bound = u;

		/* Loop over students */
		while (u<u_last) {
			
			u_bound += *nskipped++;
			
			/* Initialize accumulators */
			denom = 0;
			logpr_type = 0;
			for (j=0; j<num_covariates; j++) {
				numer[j] = 0;
				dpr_mult[j] = 0;
			}
			
			
			/* Accumulate logit denominator and gradient's numerator over skipped choices */
			while (u < u_bound) {
				expu = exp(*u++);
				for (j=0; j<num_covariates; j++) {
					numer[j] += expu*(*X++);
				}
				denom += expu;
			}
			
			
			
			/* Go over the preference list in reverse order */
			u_bound += *nlisted++;
			while (u < u_bound) {
				xb = *u++;
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
			
			pr_cur = exp(logpr_type);
			*pr_type++ = pr_cur;
			*pr += pr_cur*w_cur;
			pr++;
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
	free(u_first);
	free(dpr_type_db_first);
	free(dpr_mult);
	free(numer);
	free(dldw);
	free(dldb);

	return loglik;

}