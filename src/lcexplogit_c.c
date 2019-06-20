#include <math.h>
#include <string.h>
#include <stdint.h>

#include "mex.h"
#include "lcexplogit.h"

double lcexplogit(double *raw_param,\
        size_t num_types,\
        size_t num_covar,\
        size_t num_agents,\
        double *X,\
        uint16_t *nskipped,\
        uint16_t *nlisted,\
        double *weight,\
        double *grad)
{

	/* First element in the matrix of covariates */
	double *X_first;

	/* Same parameters in a more convenient shape: shares of types */
	double *w_t;
	
	/* Return value -- loglikelihood (matrix of one element) */
	double loglik;
	
	/* Declare constants for BLAS calls. Fortran BLAS only accepts 
     * constants by reference */
	double one = 1.0, zero = 0.0;
	char *chn = "N";
	char *cht = "T";
	size_t onei = 1;

	/* Various dimensions of data */
	size_t num_choices, cssize, csmax;
	
	/* Mean values: current, first element, exponential */
	double *xbeta, *xbeta_first, *expxbeta;

	/* Probability of observing the preference list; unconditional. First 
     * element and a movable pointer. */
	double *pr, *pr_first;
	double *prt, *prt_first; /* Same probability, but conditional on type */
	
	/* Workspace for computing components of the gradient and the logit shares */
	double *dlogprt_db, *dlogprt_db_first, *numer, \
		denom, expxb, xb, x, p_ratio, max_xbeta;
	
	/* Generic indices */
	size_t i, j, k, l, l_last;

	double *dl_dalpha, *dl_dbeta, *dl_dbeta_first;

	/* Workspace for fast access to type weights, mean utility (perhaps, compiler 
	 * optimizations make using these variables) */
	double w_cur;
	/* In case we want to accumulate individual likelihood in the log space */
	double logprt;
	
	/* Number of choices. Should match one dimension of X: total number of 
     * elements in X = num_covar*num_choices. */
	num_choices = 0;
    /* Maximum number of choices faced by one agent */
	csmax = 0;
	for (i=0; i<num_agents; i++) {
		cssize = nlisted[i] + nskipped[i];
		num_choices += cssize;
		if(csmax < cssize) {
			csmax = cssize;
		}
	}

    /*-------------------------------------------------------------------------
	 * Read and transform the parameters
	 *-----------------------------------------------------------------------*/
	w_t = (double *)malloc(num_types*sizeof(double));
	w_t[0] = 1;
	denom = 1;
	for (i=1; i<num_types; i++) {
		expxb = exp(*raw_param++);
		w_t[i] = expxb;
		denom += expxb;
	}

	for (i=0; i<num_types; i++) {
		w_t[i] /= denom;
	}
	
	/*-------------------------------------------------------------------------
	 * Allocate workspace here
	 *-----------------------------------------------------------------------*/	
	/* Array for the individual likelihood function,
     * Pr{L_i|X_i} */
	pr_first = (double *)calloc(num_agents, sizeof(double));
    
    /* Array for the individual likelihood conditional on type, 
     * Pr{L_i|X_i, t_i} */
	prt_first = (double *)malloc(num_agents*num_types*sizeof(double));
	
	/* Derivative of ln(Pr{L_i|X_i, t_i}) with respect to beta_t_i,
     * the coefficient on X in type t_i's utility function */
	dlogprt_db_first = (double *)calloc(num_agents*num_types*num_covar,\
		sizeof(double));
	
	/* Derivatives of L = sum_i ln(Pr{L_i|X_i}) with respect to the type 
     * weight parameters, \alpha_i */
	dl_dalpha = (double *)calloc(num_types, sizeof(double));

	/* Derivatives of L with respect to beta_t, type-specific coefficients 
     * on X from the utility function */
	dl_dbeta_first = (double *)calloc(num_types*num_covar, sizeof(double));
	
    /* Type-specific mean utility, beta_t*X_i, for each agent, choice and 
     * unobservable agent type */
	xbeta_first = (double *)malloc(num_choices*num_types*sizeof(double));
	
	/* First element of X */
    X_first = X;
    
    /* Variables ending with "_first" are pointers to the first element 
     * in the respective array. These pointers are not supposed to move. */
 	
    #ifdef EXTERNAL_BLAS
            openblas_set_num_threads(omp_get_num_procs());
    #endif
	
    dgemm(cht, chn, \
		&num_choices, &num_types, &num_covar, \
        &one, X, &num_covar, raw_param, &num_covar, &zero, \
        xbeta_first, &num_choices);
    
	#ifdef _OPENMP
            omp_set_num_threads(num_types);
    #endif
    
	#pragma omp parallel for \
		private(xbeta, expxbeta, max_xbeta, xb, denom, logprt, l, l_last, \
			numer, expxb, i, j, k, prt, X, x, dlogprt_db) \
		shared(X_first, prt_first, num_choices, num_covar,\
            xbeta_first, csmax, num_types, num_agents, nlisted, nskipped,\
            dlogprt_db_first, onei, one, zero, chn) \
		default(none)
	for (i=0; i<num_types; i++) {
		
        /* Initialize pointers for type i */
		prt = prt_first + i*num_agents;
		dlogprt_db = dlogprt_db_first + i*num_agents*num_covar;
		xbeta = xbeta_first + i*num_choices;
		/* Reset the X pointer */
		X = X_first;
        
        /* Allocate workspace */
        expxbeta = (double *)malloc(csmax*sizeof(double));
		numer = (double *)malloc(num_covar*sizeof(double));
		
		/* Index for the array of agent-specific choices */
		l = 0;

		/* Loop over agents */
		for (k=0; k<num_agents; k++){
			
			/* Initialize accumulators */
			denom = 0.0;
			logprt = 0.0;
			memset(numer, 0, num_covar*sizeof(double));
			
			/* Normalize max to 1 */
			l_last = nskipped[k] + nlisted[k];
			max_xbeta = xbeta[0];
			for (l=0; l<l_last; l++) {
				if (max_xbeta < xbeta[l]){
					max_xbeta = xbeta[l];
				};
			}
			for (l=0; l<l_last; l++) {
				xbeta[l] -= max_xbeta;
			}

			/* Accumulate logit denominator and gradient's numerator over 
             * skipped choices */
			l_last = nskipped[k];
			for (l=0; l<l_last; l++){
				expxbeta[l] = exp(xbeta[l]);
				denom += expxbeta[l];
			}
			dgemv(chn, &num_covar, &l_last, &one, \
				X, &num_covar, expxbeta, &onei, &zero, numer, &onei);
            
            /* Advance the pointers */
			xbeta += l_last;
            X += num_covar*l_last;
			
			/* Go over the preference list in reverse */
			l_last = nlisted[k];
			for (l=0; l<l_last; l++) {
				xb = *xbeta++;
				expxb = exp(xb);
				denom += expxb;
				for (j=0; j<num_covar; j++) {
					x = *X++;
					numer[j] += expxb*x;
					dlogprt_db[j] += x - numer[j]/denom;
				}
				logprt += xb - log(denom);
			}
			
			dlogprt_db += num_covar;
			prt[k] = exp(logprt);
		}
        
		free(numer);
		free(expxbeta);
	}
	
	dgemv(chn, &num_agents, &num_types, &one,\
            prt_first, &num_agents, w_t, &onei, &zero, pr_first, &onei);

	/*-------------------------------------------------------------------------
	 * Accumulate gradients
	 *-----------------------------------------------------------------------*/
	/* Reset the pointers */
	prt = prt_first;
	dlogprt_db = dlogprt_db_first;
	dl_dbeta = dl_dbeta_first;
	
	for (i=0; i<num_types; i++) {

		pr = pr_first;
		w_cur = w_t[i];
		
		for (k=0; k<num_agents; k++) {
			
			/* Compute conditional-to-marginal probability ratio. This is a
			 * common multiplier in all gradient terms */
			p_ratio = weight[k]*w_cur*(*prt++)/(*pr++);
			
			/* Accumulate parts of the gradient responsible for type mix */
			dl_dalpha[i] += p_ratio - weight[k]*w_cur;
			
			for (j=0; j<num_covar; j++) {
				
				dl_dbeta[j] += p_ratio*(*dlogprt_db++);
				
			}
			
		}
		
		dl_dbeta += num_covar;
		
	}
	
	
    /*-------------------------------------------------------------------------
	 * Prepare the results for output
	 *-----------------------------------------------------------------------*/
	/* Save the loglikelihood */
	pr = pr_first;
	loglik = 0;
	for (i=0; i<num_agents; i++) {
		loglik += weight[i]*log(pr[i]);
	}
	
	/* Save the gradient */
	for (i=1; i<num_types; i++) {
		*grad++ = dl_dalpha[i];
	}
	dl_dbeta = dl_dbeta_first;
	for (j=0; j<num_types*num_covar; j++) {
		*grad++ = dl_dbeta[j];
	}
	
	free(w_t);
	free(pr_first);
	free(prt_first);
	free(dlogprt_db_first);
	free(dl_dalpha);
	free(dl_dbeta);
	free(xbeta_first);

	return loglik;

}