#ifndef LCEXPLOGIT_H
#define LCEXPLOGIT_H

double lcexplogit(double *raw_param, \
	size_t num_types, \
	size_t num_covariates, \
	size_t num_agents, \
	double *X, \
	uint16_t *nskipped, \
	uint16_t *nlisted, \
	double *grad);


#ifdef EXTERNAL_BLAS

extern void dgemm_(char*, char*, size_t*, size_t*, size_t*, double*, double*, size_t*, double*, size_t*, double*, double*, size_t*);
extern void dgemv_(char *, size_t *, size_t *, double *, double *, size_t *, double *, size_t *, double *, double *, size_t *);
#define dgemm dgemm_
#define dgemv dgemv_

#else
#include "blas.h"
#endif


#endif
