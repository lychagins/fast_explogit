/* MEX interface for the latent class exploded logit likelihood function */

#include <stdint.h>
#include "mex.h"
#include "lcexplogit.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
	
	/*
	* RHS[0] = parameters, compatible with fminunc
	* RHS[1] = number of latent classes.
	* RHS[2] = demand shifters in the short stack form
	* RHS[3] = NSKIPPED; number of unlisted items from the choice set, by agent
	* RHS[4] = NLISTED; length of the preference list, by agent
    * RHS[5] = weight>=0, by student
	*/

	double *beta, *x, *loglik, *grad, *weight;
	int num_types, num_covariates, num_agents, i;
	uint16_t *nskipped, *nlisted;
	
	beta = mxGetPr(prhs[0]);
	num_types = mxGetScalar(prhs[1]);
	
	x = mxGetPr(prhs[2]);
	num_covariates = mxGetM(prhs[2]);
	nskipped = (uint16_t *)mxGetData(prhs[3]);
	nlisted = (uint16_t *)mxGetData(prhs[4]);
	num_agents = mxGetM(prhs[3]);

    if(nrhs == 6) {
		weight = mxGetPr(prhs[5]);
	} else {
		weight = (double *)malloc(num_agents*sizeof(double));
		for(i=0; i<num_agents; i++){
			weight[i] = 1.0;
		}
	}
    
	/* Output */
	plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
	loglik = mxGetPr(plhs[0]);

	/* Gradient */
	plhs[1] = mxCreateDoubleMatrix(num_covariates*num_types \
		+ num_types - 1, 1, mxREAL);
	grad = mxGetPr(plhs[1]);
	
	*loglik = lcexplogit(beta, num_types, num_covariates, num_agents, x, \
		nskipped, nlisted, weight, grad);
	return;
	
}
