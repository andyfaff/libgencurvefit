/*
 *  genetic.h
 *  GeneticOptimisation
 *
 *  Created by andrew on 27/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

/*
 we only ever want one version of the symbols in this file
 */

#ifndef GENCURVEFIT_H
#define GENCURVEFIT_H

#ifdef __cplusplus
extern "C" {
#endif


/**\mainpage
libgencurvefit is a lightweight library for data regression using differential
evolution.
To use you will have to:
 1) include gencurvefit.h.
 2) write your own fitfunction().
 3) call genetic_optimisation().

 Please note that you can create your own costfunction(), it doesn't have to be
 chisquared().

 In addition, you can specify many more options, see gencurvefitOptions.
 */

/**
 the error codes returned by this library.  They are all negative, allowing for
user errorcodes >0 to be returned from genetic_optimisation
 */


#define NO_MEMORY -1
#define INCORRECT_LIMITS -2
#define HOLDVECTOR_COEFS_MISMATCH -3
#define NO_VARYING_PARAMS -4
#define WRONG_NUMBER_OF_PARAMS -5
#define COEFS_MUST_BE_WITHIN_LIMITS -6
#define PROBLEM_CALCULATING_COVARIANCE -7
#define NO_FIT_FUNCTION_SPECIFIED -8
#define NO_Y_ARRAY -9
#define NO_X_ARRAY -10
#define NO_E_ARRAY -11
#define NO_COEFS_ARRAY -12
#define NO_LIMITS_ARRAY -13
#define	SINGULAR_MATRIX_ERROR -14

/**
 The mathematical constant Pi
 */
#define PI 3.14159265358979323846


	/**
	 Create a two-dimensional array in a single allocation
	 *
	 * The effect is the same as an array of "element p[ii][jj];
	 * The equivalent declaration is "element** p;"
	 * The array is created as an array of pointer to element, followed by an
	  array of arrays of elements.
	 * @param ii first array bound
	 * @param jj second array bound
	 * @param sz size in bytes of an element of the 2d array
	 * @return NULL on error or pointer to array
	 *
	 * assign return value to (element**)


	to use this in practice one would write

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

/**
 A function that calculates the dependent variable, given input parameters and
 independent variables.
 If you return a non-zero value from this function the fit will stop, returning
 the same error code from genetic_optimisation.

	@param userdata				- an (optional) pointer that is passed to the
	 fitfunction, costfunction and updatefunction.  Use this pointer to give
	 extra information to your functions.

	@param coefs[numcoefs]		- an array containing all the parameters for
	 calculating the model data.

	@param numcoefs				- total number of fit parameters.

	@param model[datapoints]	- the fitfunction should populate this array
	 with the model data, calculated using the coefficients.

	@param xdata[numDataDims][datapoints] - a 2D array containing the
	 independent variables that correspond to each of the datapoints. One can
	 fit multidimensional data, e.g. y = f(n, m).  In this case numDataDims = 2.
	 You can allocate a 2D dataset with m points using
	 malloc2D(2, m, sizeof(double)) (2 rows, m columns). If you want to pass in
	 a 1D dataset simply pass a pointer to the array.
         e.g. if your array is:
         double *xP;
         then pass in:
         &xP
        BUT YOU HAVE TO REMEMBER TO DEREFERENCE THE POINTER IN THE FIT FUNCTION
        BEFORE YOU USE THE ARRAY.
        model[ii] = (*xP)[ii]

	@param datapoints       - the number of datapoints to be calculated.

	@param numDataDims		- the number of independent variables in the fit.
	 For y = f(x) numDataDims = 1.  For y = f(n, m), numDataDims = 2, etc.
*/
typedef int (*fitfunction)(void *userdata,
						   const double *coefs,
						   unsigned int numcoefs,
						   double *model,
						   const double **xdata,
						   long datapoints,
						   unsigned int numDataDims);


/**
	 a function that calculates the cost function to be minimised (typically
	 chi2).
	 @param userdata				- an (optional) pointer that is passed to
	 the fitfunction, costfunction and updatefunction.  Use this pointer to give
	 extra information to your functions.

	 @param coefs[numcoefs]			- an array containing all the parameters for
	 calculating the model data.

	 @param numcoefs				- total number of fit parameters.

	 @param data[datapoints]		- the dependent variable you are trying to
	 fit

	 @param model[datapoints]		- the fitfunction will have populated this
	 array with the model data, calculated using the coefficients.

	 @param edata[datapoints]		- the error bars (assumed to be standard
	 deviation) on each of the datapoints you are trying to fit

	 @param datapoints				- the number of datapoints you are trying to
	 fit.
	 */

typedef double (*costfunction)(void *userdata,
							   const double *coefs,
							   const unsigned int numcoefs,
							   const double *data,
							   const double *model,
							   const double *errors,
							   long datapoints);

/**
 an (optional) user defined hook function to keep themselves of the fit
 progress.  If the user wishes to halt the fit early, then they should return a
 non zero value.  To keep the fit going return 0.  This will be called after
 each lowering of the best chi2 value.

 @param userdata			- an (optional) pointer that is passed to the
  fitfunction, costfunction and updatefunction.  Use this pointer to give extra
  information to your functions.

 @param coefs[numcoefs]		- an array containing all the parameters for
  calculating the model data.

 @param numcoefs			- total number of fit parameters.

 @param iterations			- how many iterations have passed.

 @param cost				- the value of the cost function (typically chi2)

 @param updatetime			- corresponds to the bitwise settings of
  gencurvefitOptions.updatefrequency

 @param convergenceNumber	- corresponds to how close the fit is to finishing
 (> 1 = finished)

 @param population	    	- the current population.  This is a 2D array having
 dimensions population[popsize][numvaryingparams].

 @param varyingparams		- those parameters that are varying. This 1D array
  is numvarying parameters long.

 @param numvaryingparams	- the number of varying parameters.

 @param totalpopsize		- how many members in the population.

 @param costmap		    	- 1D array containing all the values of the cost
  function for the population. popsize elements long.
 */

typedef int (*updatefunction)(void *userdata,
							  const double *coefs,
							  unsigned int numcoefs,
							  unsigned int iterations,
							  double cost,
							  unsigned int updatetime,
							  double convergenceNumber,
							  const double **population,
							  const unsigned int *varyingparams,
							  int numvaryingparams,
							  int totalpopsize,
							  const double *costmap);

/**
 gencurvefitOptions contains options for the genetic optimisation
 */
	struct gencurvefitOptions {
		/**
		 iterations		- the maximum number of times the population is evolved
		 during the fit (unless convergence is reached). Default = 200.
		*/
		unsigned int iterations;

		/**
		 popsizemultiplier	- the total size of the genetic population is
		 popsizemultiplier multiplied by the number of varying parameters.
		 */
		unsigned int popsizeMultiplier;

		/**
		 k_m			- the mutation constant 0 < k_m < 2.  A typical value is
		 0.7.  Make larger to get more mutation.
		 */
		double k_m;

		/**
		 dither         - If the dither array is set to two values between 0 and
		 2, then dithering is used to vary the mutation constant k_m at the
		 start of of each iteration. This can lead to improved convergence
		 properties.
		*/
		double dither[2];

		/**
		 recomb			- the recombination constant, 0 < recomb < 1.  A typical
		  value is 0.5.  Make smaller to get more exploration of parameter
		  space.
		 */
		double recomb;

		/**
		 tolerance		- specifies the stopping tolerance for the fit, which is
		  when the standard deviation of the chi2 values of the entire
		  population divided by its mean is less than tolerance.  The default
		  tolerance is 0.02
		 */
		double tolerance;

		/**
		 strategy		- Choose the Differential Evolution strategy (see
		  http://www.icsi.berkeley.edu/~storn/code.html#prac)
		 0 = Best1Bin (default);
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
		 */
		unsigned int strategy;

		/**
		 updatefun	- an (optional) function that is called each time the
		 costfunction improves.  Use this function to keep track of the fit.
		 If you return a non-zero value from this function the fit will stop.
		 */
		updatefunction updatefun;

		/**
		 updatefrequency		- Bitwise operator that specifies how often the
		 update function is called.
		 Bit No:
		 0 = everytime the fitfunction improves (default).
		 1 = everytime a monte carlo tempering move is accepted.
		 2 = after the initialisation, but before the optimisation loop starts.
		 3 = after each iteration finishes.
		 4 = after the fit has finished.
		*/
		unsigned int updatefrequency;

		/**
		seed	- seed the random number generator (must be an integer > 0).  If
		you don't want to seed put -1 in here
		 */
		int seed;

		/**
		useinitialguesses   - uses the initial guesses as a starting point
		for the fit.  If you specify useinitialguesses = 1 then the starting
		coefficients must lie in between the limits.
		 Default = useinitialguesses = 0.
		 */
		unsigned int useinitialguesses;
        
        /**
         initial_population - specifies the initial population for the fit.
         The array is assumed to be contiguous, and have total length
         `initial_population_rows * numcoefs`. The shape of the array is
         assumed to be (initial_population_rows, numcoefs), i.e. the coefficients
         for a given population member are sequential.
         This array is only used if initial_population is not NULL.
         */
        const double *initial_population;
        unsigned int initial_population_rows;

		/**
		monteCarlo          - when a fit is initialised a dataset is synthesised
		 by adding gaussian deviates to each of the y datapoints.  The gaussian
		 deviates for each point are taken from a distribution whose standard
		 deviation is equal to the error bar for that point.  This means you
		 aren't fitting the data you inputted, but a slightly altered version.
         Default = 0 (no adjustment)
		*/
		unsigned int monteCarlo;

		/**
		 polishUsingLM		- after a fit finishes the differential evolution
		  may not be at the lowest costfunction that is reachable, due to the
		  stochastic nature of the process.  If you set polishUsingLM = 1, then
		  a polishing fit using Levenberg Marquardt is carried out.
		 */
		unsigned int polishUsingLM;
	};
	typedef struct gencurvefitOptions gencurvefitOptions;


/**
 genetic_optimisation - perform curvefitting with differential evolution.
   Fitting is not limited to 1 independent variable, you can have as many as you
   like.  The function is threadsafe as long as you supply unique copies of the
   inputs to each instance. The function returns a non-zero error code (<0) if
   something goes wrong.  However, if you return a non-zero value from your fit
   function then the optimisation will stop and that value will be returned.

	@param fitfun		- a function that calculates the dependent variable,
	 given input parameters and independent variables. If you return a non-zero
	 value from this function the fit will stop.

    @param costfun		- a function that calculates the costfunction to be
     minimised.  This is normally a chi2 type function.
        i.e. sum (((model[i] - data[i]) / dataerrors[i])^2 )
        If costfun == NULL then a default chi2 function is used.

	@param numcoefs		- total number of fit parameters.

	@param coefs[numcoefs]		- an array containing all the parameters for the
	 fit. After genetic_optimisation this is populated by the parameters that best fit the data.

	@param holdvector[numcoefs]	- an array (with numcoefs elements) that
	 specifies which parameters are going to be held during the fit.
								 0 = vary
								 1 = hold

	@param limits[2][numcoefs]	- a 2D array which contains the lower and upper
	 limits for each parameter.

	@param datapoints	- the total number of data points in the fit.

	@param ydata[datapoints]	- an array containing the dependent variable
	 (i.e. the data one is trying to fit).  This array _needs_ to be datapoints
	 long.

	@param xdata[numDataDims][datapoints]  - a 2D array containing the
	 independent variables that correspond to each of the datapoints. This
	 pointer is passed directly to your fitfunction without alteration.
	 gencurvefit does not use _any_ of the values in it.  Thus, the xdata array
	 could be simply:  double **xdata = NULL;

	One can fit multidimensional data, e.g. y = f(n, m).  In this case
	numDataDims = 2. You can allocate a 2D dataset with m points using
	malloc2D(2, m, sizeof(double)). If you want to pass in a 1D dataset simply
	 pass a pointer to the array. e.g. if your array is:
                             double *xP;
		then pass in:
                             &xP
     BUT YOU HAVE TO REMEMBER TO DEREFERENCE THE POINTER IN THE FIT FUNCTION
     BEFORE YOU USE THE ARRAY.
                             model[ii] = (*xP)[ii]

	@param edata[datapoints]	- an array containing the experimental
	 uncertainties for each of the datapoints.  If you use the default chi2
	 costfunction then it should contain standard deviations.  Set each element
	 to 1 if you do not wish to weight the fit by the experimental
	 uncertainties.  This array _needs_ to be datapoints long.

	@param numDataDims		- the number of independent variables in the fit.
	 For y = f(x) numDataDims = 1.  For y = f(n, m), numDataDims = 2, etc.

	@param chi2				- the final value of the cost function.

	@param gco				- options for the genetic optimisation.  (see
	 above). If gco == NULL, then a default set of options are used.

	@param userdata			- an (optional) pointer that is passed to the
	 fitfunction, costfunction and updatefunction.  Use this pointer to give
	  extra information to your functions.
 */

int genetic_optimisation(fitfunction fitfun,
						 costfunction costfun,
						 unsigned int numcoefs,
						 double* coefs,
						 const unsigned int *holdvector,
						 const double** limits,
						 long datapoints,
						 const double* ydata,
						 const double** xdata,
						 const double *edata,
						 unsigned int numDataDims,
						 double *chi2,
						 const gencurvefitOptions *gco,
						 void* userdata);

	/**
	 does a levenberg marquardt fit to the data, instead of differential
	  evolution.  It returns a non-zero error code if something goes wrong.
	  However, it will also stop if your fitfunction returns a non-zero value.
	  As with genetic optimisation you can supply your own cost function.

	 @param fitfun		- a function that calculates the dependent variable,
	 given input parameters and independent variables. If you return a non-zero
	 value from this function the fit will stop.

	 @param costfun		- a function that calculates the costfunction to be
	 minimised.  This is normally a chi2 type function.
	 i.e. sum (((model[i] - data[i]) / dataerrors[i])^2 )
	 If costfun == NULL then a default chi2 function is used.

	 @param numcoefs	- total number of fit parameters.

	 @param coefs[numcoefs]		- an array containing all the parameters for the
	  fit. After genetic_optimisation this is populated by the parameters that
	  best fit the data.

	 @param holdvector[numcoefs]	- an array (with numcoefs elements) that
	  specifies which parameters are going to be held during the fit.
                         0 = vary
                         1 = hold

	 @param datapoints	- the total number of data points in the fit.

	 @param ydata[datapoints]	- an array containing the dependent variable
	  (i.e. the data one is trying to fit).

	 @param xdata[numDataDims][datapoints]  - a 2D array containing the
	  independent variables that correspond to each of the datapoints. One can
	  fit multidimensional data, e.g. y = f(n, m).  In this case numDataDims =
	   2. You can allocate a 2D dataset with m points using
	   malloc2D(2, m, sizeof(double)). If you want to pass in a 1D dataset
	   simply pass a pointer to the array.
       e.g. if your array is:
                    double *xP;
       then pass in:
                    &xP
       BUT YOU HAVE TO REMEMBER TO DEREFERENCE THE POINTER IN THE FIT FUNCTION
       BEFORE YOU USE THE ARRAY.
                	 model[ii] = (*xP)[ii]

	 @param edata[datapoints]	- an array containing the experimental
	  uncertainties for each of the datapoints.  If you use the default chi2
	 costfunction then it should contain standard deviations.  Set each element
	 to 1 if you do not wish to weight the fit by the experimental
	 uncertainties.

	 @param numDataDims		- the number of independent variables in the
	  fit. For y = f(x) numDataDims = 1.  For y = f(n, m), numDataDims = 2, etc.

	 @param chi2			- the final value of the cost function.

	 @param gco				- options for the genetic optimisation.  (see
	  above). If gco == NULL, then a default set of options are used.

	 @param userdata		- an (optional) pointer that is passed to the
	 fitfunction, costfunction and updatefunction.  Use this pointer to give
	 extra information to your functions.
	 */

	int levenberg_marquardt(fitfunction fitfun,
							costfunction costfun,
							unsigned int numcoefs,
							double* coefs,
							const unsigned int *holdvector,
							long datapoints,
							const double* ydata,
							const double** xdata,
							const double *edata,
							unsigned int numDataDims,
							double *chi2,
							const gencurvefitOptions* gco,
							void* userdata
							);

	/**
	 in levenbergMarquardt.c.  Calculates a hessian gradient matrix based
	 covariance matrix. The covariance matrix is returned via the
	 covarianceMatrix pointer and must be freed afterwards.

	 @param covarianceMatrix	-	a pointer to a 2D array is returned.  The
	 array (*covarianceMatrix) must be free'd afterwards.

	 @param hessianDeterminant  -	the determinant of the hessian matrix is
	 placed in this value.

	 @param userdata		-	pass in user specific information to the
	  fitfunction with this pointer.

	 @param fitfun			-	your fitfunction.
     
     @param costfun         -   your costfunction

	 @param coefs[numcoefs]	-	an array containing the coefficients.  The
	  covariance matrix is assessed for these values.

	 @param numcoefs		-	the number of coefficients.

	 @param holdvector[numcoefs]	-	an array specifying which parameters
	  were held (=1) or varied (=0) during the fit.

	 @param ydata[datapoints]	-	an array of the data being fitting.

	 @param edata[datapoints]	-	an array for the error bars for the data
	  being fitted.

	 @param datapoints			-	the number of datapoints being fitted.

	 @param xdata[numDataDims][datapoints]	-	an array containing the
	  independent variables for the fit.

	 @param numDataDims 		-	how many independent variables do you have?

	 @param unitSD				-	specify as 1 if the datapoints were unit
	  weighted.

	 */
	int getCovarianceMatrix(double ***covarianceMatrix,
							double *hessianDeterminant,
							void *userdata,
							fitfunction fitfun,
                            costfunction costfun,
							double *coefs,
							int numcoefs,
							unsigned int *holdvector,
							const double *ydata,
							const double *edata,
							const double **xdata,
							long datapoints,
							int numDataDims,
							int unitSD);

	/**
	 a default chi2 cost function

	 @param userdata	- an (optional) pointer that is passed to the
	  fitfunction, costfunction and updatefunction.  Use this pointer to give
	  extra information to your functions.

	 @param coefs[numcoefs]		- an array containing the coefficients for the
	  fit.

	 @param numcoefs			- the number of parameters being fitted.

	 @param data[datapoints]	- the data points being fitted.

	 @param model[datapoints]	- the model values calculated by the
	  fitfunction.

	 @param errors[datapoints]	- the error bars (standard deviation)
	  corresponding to each of the datapoints.

	 @param datapoints			- the number of datapoints being fitted.
	 */


double chisquared(void *userdata,
                  const double *coefs,
                  unsigned int numcoefs,
                  const double *data,
                  const double *model,
                  const double *errors,
                  long datapoints);

	/**
	 a default robust cost function

	 @param userdata	- an (optional) pointer that is passed to the
	  fitfunction, costfunction and updatefunction.  Use this pointer to give
      extra information to your functions.

	 @param coefs[numcoefs]		- an array containing the coefficients for the
	  fit.

	 @param numcoefs			- the number of parameters being fitted.

	 @param data[datapoints]	- the data points being fitted.

	 @param model[datapoints]	- the model values calculated by the
	  fitfunction.

	 @param errors[datapoints]	- the error bars (standard deviation)
	  corresponding to each of the datapoints.

	 @param datapoints			- the number of datapoints being fitted.

	 */
double robust(void *userdata,
              const double *coefs,
              unsigned int numcoefs,
              const double *data,
              const double *model,
              const double *errors,
              long datapoints);


#ifdef __cplusplus
}
#endif

#endif
