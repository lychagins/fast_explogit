#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include "gperftools/profiler.h"
#include "explogit.h"

#define NUM_CLASSES 20
#define NUM_X 140526176
#define NUM_COVAR 14
#define NUM_STUDENTS 7875
#define NUM_BETA 299

int main(){
	
	int *nlisted, *nskipped;
	clock_t t;
	double *x, *beta, *grad, sec, logl;
	FILE *fp;

	ProfilerStart("explogit.profile");
	
	nlisted = (int *)malloc(NUM_STUDENTS*sizeof(int));
	nskipped = (int *)malloc(NUM_STUDENTS*sizeof(int));
	x = (double *)malloc(NUM_X*sizeof(double));
	beta = (double *)malloc(NUM_BETA*sizeof(double));
	// Space for storing the gradient
	grad = (double *)malloc((NUM_COVAR*NUM_CLASSES + NUM_CLASSES - 1)*sizeof(double));
	
	/* Load data, time */
	t = clock();

	fp = fopen("../Data/nlisted.bin", "rb");
	fread(nlisted, sizeof(int), NUM_STUDENTS, fp);
	fclose(fp);

	fp = fopen("../Data/nskipped.bin", "rb");
	fread(nskipped, sizeof(int), NUM_STUDENTS, fp);
	fclose(fp);

	fp = fopen("../Data/x.bin", "rb");
	fread(x, sizeof(double), NUM_X, fp);
	fclose(fp);

	fp = fopen("../Data/beta.bin", "rb");
	fread(beta, sizeof(double), NUM_BETA, fp);
	fclose(fp);
	
	t = clock() - t;
	sec = ((double)t)/CLOCKS_PER_SEC;

	printf("Finished reading the data. Time spent = %.2lf sec.\n", sec);

	t = clock();
	logl = explogit(beta, NUM_CLASSES, NUM_COVAR, NUM_STUDENTS, x, nskipped, nlisted, grad);
	t = clock() - t;
	sec = ((double)t)/CLOCKS_PER_SEC;

	printf("Log-likelihood = %lf. Time spent = %.2lf sec.\n", logl, sec);
	ProfilerStop();
	return 0;
}
