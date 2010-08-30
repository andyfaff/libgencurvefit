/*
 *  genetic.h
 *  GeneticOptimisation
 *
 *  Created by andrew on 27/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#ifdef __cplusplus
extern "C" {
#endif
	
#define NO_MEMORY -1
#define INCORRECT_LIMITS -2
#define HOLDVECTOR_COEFS_MISMATCH -3
#define NO_VARYING_PARAMS -4
#define WRONG_NUMBER_OF_PARAMS -5

#define PI 3.14159265358979323846

	/*
	 Create a two-dimensional array in a single allocation
	 *
	 * The effect is the same as an array of "element p[ii][jj];
	 * The equivalent declaration is "element** p;"
	 * The array is created as an array of pointer to element, followed by an array of arrays of elements.
	 * \param ii first array bound
	 * \param jj second array bound
	 * \param sz size in bytes of an element of the 2d array
	 * \return NULL on error or pointer to array
	 *
	 * assign return value to (element**)
	 */
	
	/* to use this in practice one would write 
	 
	 double **pp = NULL;
	 pp = (double**)malloc2d(5, 11, sizeof(double));
	 if(pp==NULL)
	 return NOMEM;
	 
	 <use pp as required>
	 free(pp);
	 
	 Note you can access elements by
	 *(*(p+i)+j) is equivalent to p[i][j]
	 In addition *(p+i) points to a whole row.
	 */
	
	void* malloc2d(int ii, int jj, int sz);
	
/*
 a function that calculates the dependent variable, given input parameters and independent variables. 
 If you return a non-zero value from this function the fit will stop, returning the same error code from genetic_optimisation.

	userdata				- an (optional) pointer that is passed to the fitfunction, costfunction and updatefunction.  Use this pointer to give extra
								information to your functions.
 
	coefs[numcoefs]			- an array containing all the parameters for calculating the model data.
 
	numcoefs				- total number of fit parameters.
 
	model[numpnts]			- the fitfunction should populate this array with the model data, calculated using the coefficients.
 
	xdata[numDataDims][numpnts] - a 2D array containing the independent variables that correspond to each of the datapoints.
									 One can fit multidimensional data, e.g. y = f(n, m).  In this case numDataDims = 2.
									 You can allocate a 2D dataset with m points using malloc2D(2, m, sizeof(double)) (2 rows, m columns)
									 If you want to pass in a 1D dataset simply pass a pointer to the array.
									 e.g. if your array is:
									 double *xP;
									 then pass in:
									 &xP
									BUT YOU HAVE TO REMEMBER TO DEREFERENCE THE POINTER IN THE FIT FUNCTION BEFORE YOU USE THE ARRAY.
									model[ii] = (*xP)[ii]
 
 numpnts					- the number of datapoints to be calculated.
	
	numDataDims				- the number of independent variables in the fit. For y = f(x) numDataDims = 1.  For y = f(n, m), numDataDims = 2, etc.
*/
typedef int (*fitfunction)(void *userdata, const double *coefs, int numcoefs, double *model, const double **xdata, long numpnts, int numDataDims);
	
	
typedef double (*costfunction)(void *userdata, const double *params, int numparams, const double *data, const double *model, const double *errors, long numpnts);

/*
 an (optional) user defined hook function to keep themselves of the fit progress.  If the user wishes to halt the fit early, then they should return a non
 zero value.  To keep the fit going return 0.  This will be called after each lowering of the best chi2 value.
 */
typedef int (*updatefunction)(void *userdata, const double *coefs, unsigned int numcoefs, unsigned int iterations, double cost);

/*
 genetic_optimisation - perform curvefitting with differential evolution.  Fitting is not limited to 1 independent variable,
  you can have as many as you like.  The function is threadsafe as long as you supply unique copies of the inputs to each instance.
 
	fitfun					- a function that calculates the dependent variable, given input parameters and independent variables. 
								If you return a non-zero value from this function the fit will stop. 
 
    costfun					- a function that calculates the costfunction to be minimised.  This is normally a chi2 type function.
								i.e. sum (((model[i] - data[i]) / dataerrors[i])^2 )
								If costfun == NULL then a default chi2 function is used.
 
	numcoefs				- total number of fit parameters.
 
	coefs[numcoefs]			- an array containing all the parameters for the fit.  After genetic_optimisation this is populated by the parameters
								that best fit the data.
 
	holdvector[numcoefs]	- an array (with numcoefs elements) that specifies which parameters are going to be held during the fit. 
								 0 = vary
								 1 = hold

	limits[numcoefs][2]		- a 2D array which contains the lower and upper limits for each parameter. The lower limit must be lower than the upper limit,
								but only for those parameters that are being varied.
 
	datapoints				- the total number of data points in the fit.
 
	ydata[datapoints]		- an array containing the dependent variable (i.e. the data one is trying to fit).

	xdata[numDataDims][datapoints]  - a 2D array containing the independent variables that correspond to each of the datapoints.
										One can fit multidimensional data, e.g. y = f(n, m).  In this case numDataDims = 2.
										You can allocate a 2D dataset with m points using malloc2D(2, m, sizeof(double)).
										If you want to pass in a 1D dataset simply pass a pointer to the array.
										e.g. if your array is:
										 double *xP;
										 then pass in:
										 &xP
										 BUT YOU HAVE TO REMEMBER TO DEREFERENCE THE POINTER IN THE FIT FUNCTION BEFORE YOU USE THE ARRAY.
										 model[ii] = (*xP)[ii]
 
	edata[datapoints]		- an array containing the experimental uncertainties for each of the datapoints.  If you use the default chi2 costfunction
								then it should contain standard deviations.  Set each element to 1 if you do not wish to weight the fit by the experimental
								uncertainties.  
 
	numDataDims				- the number of independent variables in the fit. For y = f(x) numDataDims = 1.  For y = f(n, m), numDataDims = 2, etc.
 
	chi2					- the final value of the cost function.
 
	iterations				- the maximum number of times the population is evolved during the fit (unless convergence is reached).
 
	popsizemultiplier		- the total size of the genetic population is popsizemultiplier multiplied by the number of varying parameters.
 
	k_m						- the mutation constant 0 < k_m < 2.  A typical value is 0.7.  Make larger to get more mutation.
 
	recomb					- the recombination constant, 0 < recomb < 1.  A typical value is 0.5.  Make smaller to get more exploration of parameter space.
 
	tolerance				- specifies the stopping tolerance for the fit, which is when the standard deviation of the chi2 values of the entire population
								divided by its mean is less than tolerance.
 
	strategy				- Choose the Differential Evolution strategy (see http://www.icsi.berkeley.edu/~storn/code.html#prac)
								 0 = Best1Bin;
								 1 = Best1Exp;
								 2 = Rand1Exp;
								 3 = RandToBest1Exp;
								 4 = Best2Exp;
								 5 = Rand2Exp;
								 6 = RandToBest1Bin;
								 7 = Best2Bin;
								 8 = Rand2Bin;
								 9 = Rand1Bin;
								 Try Best1Bin to start with.
 
	MCtemp					- Normally if the chi2 value of the trial vector is lower than vector i from the population then the trial vector replaces vector i. 
								However, if you specify temp  is specified then the probability of the trial vector being accepted is now done on a Monte Carlo basis. I.e.:
								accept if
									chi2(trial) < chi2(i)
								or accept if
								exp(-chi2(trial) / chi2(i) / MCtemp) < enoise(1) 
								This has the effect of exploring wider parameter space, and is more likely to find a global minimum, but may take longer to converge. 
								One should use more iterations with MCtemp. If one records the history of the fit using updatefun, then one can use the history for use in calculating
								a covariance matrix or use as the posterior probability distribution for Bayesian model selection.  
								IF YOU DON'T WANT THIS TEMPERING SET MCtemp TO NAN (e.g. sqrt(-1)) .
 
	updatefun				- an (optional) function that is called each time the costfunction improves.  Use this function to keep track of the fit.
								If you return a non-zero value from this function the fit will stop. This function will also be called if a move is accepted on a monte carlo basis (see MCtemp). 
 
	userdata				- an (optional) pointer that is passed to the fitfunction, costfunction and updatefunction.  Use this pointer to give extra
								information to your functions.
 */

//contains options for the genetic optimisation
struct gencurvefitOptions {
	int iterations;
	int popsizeMultiplier;
	double k_m;
	double recomb;
	double tolerance;
	unsigned int strategy;
	double temp;
	updatefunction updatefun;
};
typedef struct gencurvefitOptions gencurvefitOptions;
	
int genetic_optimisation(fitfunction fitfun,
						 costfunction costfun,
						 unsigned int numcoefs,
						 double* coefs,
						 const int *holdvector,
						 const double** limits,
						 long datapoints,
						 const double* ydata,
						 const double** xdata,
						 const double *edata,
						 int numDataDims,
						 double *chi2,
						 const gencurvefitOptions *gco,
						 void* userdata
						 );

double gnoise(double sd);

double chisquared(void *userdata, const double *params, int numparams, const double *data, const double *model, const double *errors, long numpnts);
double robust(void *userdata, const double *params, int numparams, const double *data, const double *model, const double *errors, long numpnts);


#ifdef __cplusplus
}
#endif
