#ifndef EXPLODED_LOGIT
#define EXPLODED_LOGIT

void dgemv(double *X, size_t nX, double *b, size_t nb, double *Xb);
void dgemv_t(double *X, size_t nX, double *b, size_t nb, double *Xb);
double explogit(double *beta, size_t num_covariates, size_t num_students, double *X, uint16_t *nskipped, uint16_t *nlisted, double *weight, double *grad);

#endif
