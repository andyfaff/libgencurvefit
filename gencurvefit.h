/*
 *  genetic.h
 *  GeneticOptimisation
 *
 *  Created by andrew on 27/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#define NO_MEMORY -1
#define INCORRECT_LIMITS -2
#define HOLDVECTOR_COEFS_MISMATCH -3
#define NO_VARYING_PARAMS -4
#define WRONG_NUMBER_OF_PARAMS -5

#define PI 3.14159265358979323846

/*function type for a typical fitfunction
return 0 for no error.
return an integer != 0 for an error.
if you return an error the fit will finish immediately.
This error will be returned from genetic_optimisation
*/
typedef int (*fitfunction)(void *userdata, const double *coefs, double *model, const double **xdata, long numpnts, int numDataDims);

typedef double (*costfunction)(void *userdata, const double *params, int numparams, const double *data, const double *model, const double *errors, long numpnts);

/*
 an (optional) user defined hook function to keep themselves of the fit progress.  If the user wishes to halt the fit early, then they should return a non
 zero value.  To keep the fit going return 0.  This will be called after each lowering of the best chi2 value.
 */
typedef int (*updatefunction)(void *userdata, const double *coefs, unsigned int numcoefs, unsigned int iterations, double cost);

int genetic_optimisation(fitfunction fitfun,
						 costfunction costfun,
						 unsigned int numcoefs,
						 double* coefs,
						 long datapoints,
						 const int *holdvector,
						 const double** limits,
						 int numDataDims,
						 const double* ydata,
						 const double** xdata,
						 const double *edata,
						 double *chi2,
						 int iterations,
						 int popsizeMultiplier,
						 double k_m,
						 double recomb,
						 double tolerance,
						 unsigned int strategy,
						 updatefunction updatefun,
						 void* userdata
						 );

double gnoise(double sd);

double chisquared(void *userdata, const double *params, int numparams, const double *data, const double *model, const double *errors, long numpnts);
double robust(void *userdata, const double *params, int numparams, const double *data, const double *model, const double *errors, long numpnts);



