/* MEX interface for the exploded logit likelihood function */
#include <stdint.h>
#include "mex.h"
#include "explogit.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	/*
	* RHS[0] = parameters, compatible with fminunc
	* RHS[1] = covariates in the short stack form
	* RHS[2] = NSKIPPED; number of unlisted items from the choice set, by agent
	* RHS[3] = NLISTED; length of the preference list, by agent
	* RHS[4] = weight>=0; by student
	*/

	double *beta, *x, *loglik, *grad, *weight;
	size_t num_covariates, num_students, i;
	uint16_t *nskipped, *nlisted;
	
	beta = mxGetPr(prhs[0]);
	
	x = mxGetPr(prhs[1]);
	num_covariates = mxGetM(prhs[1]);
	nskipped = (uint16_t *)mxGetData(prhs[2]);
	nlisted = (uint16_t *)mxGetData(prhs[3]);
	num_students = mxGetM(prhs[2]);
	
	if(nrhs == 5) {
		weight = mxGetPr(prhs[4]);
	} else {
		weight = (double *)malloc(num_students*sizeof(double));
		for(i=0; i<num_students; i++){
			weight[i] = 1.0;
		}
	}
	
	/* Output */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	loglik = mxGetPr(plhs[0]);

	/* Gradient */
	plhs[1] = mxCreateDoubleMatrix(num_covariates, 1, mxREAL);
	grad = mxGetPr(plhs[1]);
	
	*loglik = explogit(beta, num_covariates, num_students, x, nskipped, nlisted, weight, grad);
	
	/* Clean up */
	if(nrhs < 5 && weight!=NULL){
		free(weight);
	}
	
	return;
}
