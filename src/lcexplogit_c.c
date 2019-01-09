#include <stdlib.h>
#include <stdio.h>
#include <math.h>
//#include <omp.h>
#include <cblas.h>
#include <string.h>
#include <stdint.h>
#include "lcexplogit.h"

double lcexplogit(double *raw_param, int num_types, int num_covariates, int num_agents,\
	double *X, uint16_t *nskipped, uint16_t *nlisted, double *grad)
{

	/* Pointers for the matrix of covariates: first/last/current elements */
	double *X_first;

	/* Same parameters in a more convenient shape: shares of types */
	double *w_t;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;

	/* Various dimensions of data */
	size_t num_choices, cssize, csmax;
	
	/* Pointers for the vector of mean values: first/last/current/pref. list boundary */
	double *u, *u_first;

	/* Probability of observing the preference list; unconditional. First element and
	 *  a movable pointer. */
	double *pr, *pr_first;
	double *pr_type, *pr_type_first; /* Same probability, but conditional on type */
	
	/* Workspace for computing components of the gradient and the logit shares */
	double *dpr_type_db, *dpr_type_db_first, *numer, *dpr_mult, \
		denom, expu, xb, x, p_ratio, max_u;
	
	/* Generic indices */
	size_t i, j, k, l, l_last;

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
	csmax = 0;
	for (i=0; i<num_agents; i++) {
		cssize = nlisted[i] + nskipped[i];
		num_choices += cssize;
		if(csmax < cssize) {
			csmax = cssize;
		}
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
	pr_first = (double *)calloc(num_agents, sizeof(double));
	pr_type_first = (double *)malloc(num_agents*num_types*sizeof(double));
	
	/* Derivatives of the conditional choice probability, by student and type */
	dpr_type_db_first = (double *)calloc(num_agents*num_types*num_covariates,\
		sizeof(double));
	
	/* This is where we accumulate the gradient w.r.t. the type weight parameters */
	dldw = (double *)calloc(num_types, sizeof(double));

	/* This is where we accumulate the gradient w.r.t. to the type-specific utility parameters */
	dldb_first = (double *)calloc(num_types*num_covariates, sizeof(double));
	
	u_first = (double *)malloc(num_choices*num_types*sizeof(double));
	
	X_first = X;
 	/* openblas_set_num_threads(omp_get_num_procs()); */
	cblas_dgemm(CblasColMajor, CblasTrans, CblasNoTrans, \
		num_choices, num_types, num_covariates, 1, X, num_covariates, \
		raw_param, num_covariates, 0, u_first, num_choices);
	/* omp_set_num_threads(num_types); */
	#pragma omp parallel for \
		private(u, xb, denom, logpr_type, l, l_last, \
			numer, dpr_mult, expu, i, j, k, pr_type, X, x, dpr_type_db, max_u) \
		shared(X_first, pr_type_first, num_choices, num_covariates, u_first, csmax, \
			num_types, num_agents, nlisted, nskipped, dpr_type_db_first) \
		default(none)
	for (i=0; i<num_types; i++) {
		/* Initialize pointers for type i */
		pr_type = pr_type_first + i*num_agents;
		dpr_type_db = dpr_type_db_first + i*num_agents*num_covariates;
		u = u_first + i*num_choices;
		double *v = (double *)malloc(csmax*sizeof(double));
		
		numer = (double *)malloc(num_covariates*sizeof(double));
		dpr_mult = (double *)malloc(num_covariates*sizeof(double));

		/* Reset X pointer */
		X = X_first;
		
		/* Index for the "short stack" array */
		l = 0;

		/* Loop over students */
		for (k=0; k<num_agents; k++){
			
			/* Initialize accumulators */
			denom = 0.0;
			logpr_type = 0.0;
			memset(dpr_mult, 0, num_covariates*sizeof(double));
			memset(numer, 0, num_covariates*sizeof(double));
			
			/* Normalize max to 1 */
			l_last = nskipped[k] + nlisted[k];
			max_u = u[0];
			for (l=0; l<l_last; l++) {
				if (max_u < u[l]){
					max_u = u[l];
				};
			}
			for (l=0; l<l_last; l++) {
				u[l] -= max_u;
			}

			/* Accumulate logit denominator and gradient's numerator over skipped choices */
			l_last = nskipped[k];
			for (l=0; l<l_last; l++){
				v[l] = exp(u[l]);
				denom += v[l];
			}
			u += l_last;

			cblas_dgemv(CblasColMajor, CblasNoTrans, num_covariates, l_last, \
				1.0, X, num_covariates, v, 1, 0.0, numer, 1);
			X += num_covariates*l_last;
			
			/* Go over the preference list in reverse order */
			l_last = nlisted[k];
			for (l=0; l<l_last; l++) {
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
			
			memcpy(dpr_type_db, dpr_mult, num_covariates*sizeof(double));
			dpr_type_db += num_covariates;
			
			pr_type[k] = exp(logpr_type);
		}
				
		free(dpr_mult);
		free(numer);
		free(v);
	}
	
	cblas_dgemv(CblasColMajor, CblasNoTrans, num_agents, num_types, \
				1.0, pr_type_first, num_agents, w_t, 1, 0.0, pr_first, 1);

	/*-------------------------------------------------------------------------
	 * Accumulate gradients
	 *-----------------------------------------------------------------------*/
	/* Reset the pointers */
	pr_type = pr_type_first;
	dpr_type_db = dpr_type_db_first;
	dldb = dldb_first;
	
	for (i=0; i<num_types; i++) {

		pr = pr_first;
		w_cur = w_t[i];
		
		for (k=0; k<num_agents; k++) {
			
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
	/* Save the loglikelihood */
	pr = pr_first;
	loglik = 0;
	for (i=0; i<num_agents; i++) {
		loglik += log(pr[i]);
	}
	
	/* Save the gradient */
	
	for (i=1; i<num_types; i++) {
		*grad++ = dldw[i];
	}
	dldb = dldb_first;
	for (j=0; j<num_types*num_covariates; j++) {
		*grad++ = dldb[j];
	}
	
	free(w_t);
	free(pr_first);
	free(pr_type_first);
	free(dpr_type_db_first);
	free(dldw);
	free(dldb);
	free(u_first);

	return loglik;

}