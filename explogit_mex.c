/* This is a MEX interface for the exploded logit likelihood function */

#include <stdint.h>
typedef uint16_t char16_t;

#include "mex.h"
#include "explogit.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	/*
	* RHS[0] = parameters, compatible with fminunc
	* RHS[1] = number of unobservable types.
	* RHS[2] = demand shifters in the short stack form
	* RHS[3] = NSKIPPED; number of unlisted items from the choice set, by student
	* RHS[4] = NLISTED; length of the preference list, by student
	*/

	double *beta, *x, *loglik, *grad;
	int num_types, num_covariates, num_students;
	unsigned int *nskipped, *nlisted;
	
	beta = mxGetPr(prhs[0]);
	num_types = mxGetScalar(prhs[1]);
	x = mxGetPr(prhs[2]);
	num_covariates = mxGetM(prhs[2]);
	nskipped = (unsigned int *)mxGetData(prhs[3]);
	nlisted = (unsigned int *)mxGetData(prhs[4]);
	num_students = mxGetM(prhs[3]);

	
	/* Output */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	loglik = mxGetPr(plhs[0]);

	/* Gradient */
	plhs[1] = mxCreateDoubleMatrix(num_covariates*num_types \
		+ num_types - 1, 1, mxREAL);
	grad = mxGetPr(plhs[1]);
	
	*loglik = explogit(beta, num_types, num_covariates, num_students, x, \
		nskipped, nlisted, grad);
	return;
	
}
