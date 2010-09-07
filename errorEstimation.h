

/*
 functions contained in errorEstimation.c
 */
int getCovarianceMatrix();
static int choldc (double **a, int N, double *p);
static void cholsl(double **a, int N, const double *p, double *b, double *x);

