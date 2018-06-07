#ifndef EXPLOGIT_HEADER
#define EXPLOGIT_HEADER

double explogit(double *raw_param, int num_types, int num_covariates, int num_students, double *X,\
	int *nskipped, int *nlisted, double *u_first, double *grad);

#endif
