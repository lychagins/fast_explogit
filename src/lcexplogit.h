#ifndef LCEXPLOGIT_HEADER
#define LCEXPLOGIT_HEADER

double lcexplogit(double *raw_param, int num_types, int num_covariates, int num_students, double *X, uint16_t *nskipped, uint16_t *nlisted, double *grad);

#endif
