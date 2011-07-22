
/*
 *  errorEstimation.c
 * 
 *	This code works out the covariance matrix, and therefore the fit errors for the genetic curvefit.  It does so by a matrix method
 *	THis matrix method does not affect the fit parameters in any way, but it is a gradient technique.
 *
 *  Created by andrew on 24/09/07.
 *  Copyright 2007 __Andrew Nelson and The Australian Nuclear Science and Technology Organisation__. All rights reserved.
 *
 */
#include "gencurvefit.h"
#include "math.h"
#include "stdlib.h"
#include "string.h"

#define TINY 1.0e-20
#define EPSILON 1.0e-6

static double factorial(double num){
	int ii;
	double result = 0;
	
	if( num<70){
		result = 1;
		for(ii = 1 ; ii < num+1 ; ii+=1){
			result *= (double)ii;
		}
	} else {
		result = sqrt(2 * 3.14159 * num) * pow(num/2.71828,num);
	}
	return result;
}

/**
 perform LU decomposition
 */

static int ludcmp(double **a, int n, int *indx, double *d){
	int i, imax = 0, j, k, err = 0;
	double big, dum, sum, temp;
	double *vv = NULL;
	
	vv = (double*)malloc(sizeof(double)*n);
	if(vv == NULL){
		err = 0;
		goto done;
	}
	
	*d = 1.0;
	
	for(i=0 ; i<n ; i++){
		big = 0.0;
		for(j=0 ; j<n ; j++)
			if((temp = fabs(a[i][j])) > big) big = temp;
		if(big == 0.0){
			err = SINGULAR_MATRIX_ERROR;
			goto done;
		}
		vv[i] = 1.0/big;	
	} 
	
	for(j=0 ; j<n ; j++){
		for(i=0 ; i<j ; i++){
			sum = a[i][j];
			for (k=0 ; k<i ; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
		}
		big = 0.0;
		for(i=j ; i<n ; i++){
			sum = a[i][j];
			for(k=0; k<j ; k++) sum -= a[i][k]*a[k][j];
			a[i][j] = sum;
			if( (dum=vv[i]*fabs(sum)) >= big) {
				big = dum;
				imax = i;
			}
		}
		if(j != imax){
			for(k=0 ; k<n ; k++){
				dum = a[imax][k];
				a[imax][k] = a[j][k];
				a[j][k] = dum;
			}
			*d = -(*d);
			vv[imax] = vv[j];
		}
		indx[j] = imax;
		if(a[j][j] == 0.0) a[j][j] = TINY;
		
		if(j != n-1){
			dum = 1.0/(a[j][j]);
			for(i=j+1; i<n ; i++) a[i][j] *=dum;
		}
	}
	
	
done:
	if(vv != NULL)
		free(vv);
	
	return err;
	
}

/**
 perform the backsubstitution for LU decomposition
 */
static void lubksb(double **a, int n, int *indx, double b[]){
	int i, ii=0, ip, j;
	double sum;
	
	for(i=1 ; i<=n ; i++){
		ip = indx[i-1];
		sum = b[ip];
		b[ip] = b[i-1];
		if(ii)
			for(j=ii ; j<=i-1 ; j++) sum -= a[i-1][j-1]*b[j-1];
		else if (sum) ii=i;
		b[i-1] = sum;
	}
	for(i=n ; i>=1 ; i--){
		sum = b[i-1];
		for(j=i+1 ; j<=n ; j++) sum -= a[i-1][j-1]*b[j-1];
		b[i-1] = sum/a[i-1][i-1];
	}	
}

/**
 Cholesky Decomposition, only useful for Symmetric positive definite matrices.
 */
static int choldc (double **a, int N, double *p){
	int err = 0;	
	int ii, jj, kk;
	double sum = 0;
	
	
	for(ii = 0; ii < N ; ii++){
		for(jj = ii ; jj < N ; jj++){
			for(sum = a[ii][jj] , kk = ii - 1 ; kk >= 0; kk--) sum -= a[ii][kk] * a[jj][kk];
			if(ii == jj){
				if(sum <= 0.0)
					return PROBLEM_CALCULATING_COVARIANCE;
				p[ii] = sqrt(sum);
			} else a[jj][ii] = sum/p[ii];
		}
	}
	
	return err;
}


/**
 Cholesky back substitution
 */
static void cholsl(double **a, int N, const double *p, double *b, double *x){
	int ii, kk;
	double sum = 0;
	for(ii = 0 ; ii < N ; ii++){
		for(sum = b[ii], kk = ii-1 ; kk >=0 ; kk-- ) sum -=a[ii][kk] * x[kk];
		x[ii] = sum/p[ii];
	}
	for(ii = N-1 ; ii >= 0 ; ii--){
		for(sum = x[ii], kk = ii+1 ; kk < N ; kk++ ) sum -=a[kk][ii] * x[kk];
		x[ii] = sum/p[ii];
	}
}

/**
 LU matrix inversion
 */
int matrixInversion_lu(double **a, int N){
	int err=0;
	int i,j;
	int *indx = NULL;
	double *col = NULL;
	double **tempA = NULL;
	double d;
	
	
	indx = (int*)malloc(sizeof(int)*N);
	if(indx == NULL){
		err = 0;
		goto done;
	}
	
	tempA = (double**)malloc2d(N,N,sizeof(double));
	if(tempA == NULL){
		err = 0;
		goto done;
	}
	for(i=0; i<N; i++){
		for(j=0 ; j<N ; j++){
			tempA[i][j] = a[i][j];
		}	
	}
	
	col = (double*)malloc(sizeof(double)*N);
	if(col == NULL){
		err = 0;
		goto done;
	}
	
	if((err = ludcmp(tempA,N,indx,&d)))
		goto done;
	
	for(j=0 ; j<N ; j++){
		for(i=0 ; i<N ; i++) col[i] = 0.0;
		col[j] = 1.0;
		lubksb(tempA,N,indx,col);
		for(i=0 ; i<N ; i++){
			d = col[i];
			a[i][j] = col[i];
		};
	}
done:
	return err;
}

/**
 cholesky matrix inversion (only useful for symmetric positive definite matrices)
 */
static int matrixInversion_chol(double **a, int N, double *detA){
	int err=0;
	int i,j;
	double *x = NULL;
	double *b = NULL;
	double *p = NULL;
	double **tempA = NULL;
	if(detA)
		*detA = 1;
	
	x = (double*)malloc(sizeof(double) * N);
	if(x == NULL){
		err = NO_MEMORY;
		goto done;
	}
	
	tempA = (double**)malloc2d(N, N, sizeof(double));
	if(tempA == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memcpy(tempA, a, sizeof(double) * N * N);
	
	p = (double*)malloc(sizeof(double) * N);
	if(p == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memset(p, 0, sizeof(double) * N);
	
	b = (double*)malloc(sizeof(double) * N);
	if(b == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memset(b, 0, sizeof(double) * N);
	
	//perform the cholesky decomposition
	if((err = choldc(tempA, N, p))) goto done;
	
	//now do the back substitution
	for(j = 0 ; j < N ; j++){
		memset(b, 0, sizeof(double) * N);
		b[j] = 1.0;
		memset(x, 0, sizeof(double) * N);
		cholsl(tempA, N, p, b, x);
		
		for(i = 0 ; i < N ; i++){
			a[i][j] = x[i];
		}
	}
	
	//make the covariance matrix symmetric
	for(i=0 ; i < N ; i++)
		for(j = N-1 ; j > i ; j--)
			a[i][j] = a[j][i];
	
	
	//the determinant of the original matrix is the square of the products of the elements in the cholesky diagonal
	for(i = 0 ; i < N && detA ; i+=1)
		*detA *= tempA[i][i] * tempA[i][i];
	
done:
	if(p)
		free(p);
	if(b)
		free(b);
	if(x)
		free(x);
	if(tempA)
		free(tempA);
	return err;
}


static int partialDerivative(void *userdata, fitfunction fitfun, double** derivativeMatrix, int derivativeMatrixRow, int parameterIndex, const double* coefs, int numcoefs, const double **xdata, long datapoints, int numDataDims){
	int err = 0;
	double param, diff;
	int jj;
	double *dataTemp = NULL;
	double *coefs_temp = NULL;
	
	coefs_temp = (double*)malloc(sizeof(double) * numcoefs);
	 if(!coefs_temp){
		err = NO_MEMORY;
		goto done;
	 }
	memcpy(coefs_temp, coefs, numcoefs * sizeof(double));

	param = coefs[parameterIndex];	
	diff = EPSILON * param;
	coefs_temp[parameterIndex] = param + diff;
	
	if((err = fitfun(userdata, coefs_temp, numcoefs, *(derivativeMatrix + derivativeMatrixRow), (const double**)xdata, datapoints, numDataDims)))
		goto done;
	
	coefs_temp[parameterIndex] = param - diff;
	
	dataTemp = (double*)malloc(sizeof(double) * datapoints);
	if(!dataTemp){
		err = NO_MEMORY;
		goto done;
	}
	
	if((err = fitfun(userdata, coefs_temp, numcoefs, dataTemp, (const double**)xdata, datapoints, numDataDims)))
		goto done;
	
	for(jj = 0 ; jj < datapoints ; jj++)
		derivativeMatrix[derivativeMatrixRow][jj] = (derivativeMatrix[derivativeMatrixRow][jj] - dataTemp[jj]) / (2 * diff);
	
done:
	if(coefs_temp)
		free(coefs_temp);
	if(dataTemp)
		free(dataTemp);
	
	return err;
}


static int updatePartialDerivative(void *userdata, fitfunction fitfun, double **derivativeMatrix, double *coefs, int numcoefs, unsigned int *varparams, int numvarparams, const double **xdata, long datapoints, int numDataDims){
	int err = 0;
	int ii;
	for(ii = 0 ; ii < numvarparams ; ii++){
		if((err = partialDerivative(userdata, fitfun, derivativeMatrix, ii, varparams[ii], coefs, numcoefs, xdata, datapoints, numDataDims)))
			return err;
	} 
	return err;
}


void calculateAlphaElement(int row, int col, double **alpha, double **derivativeMatrix, const double *edata, long datapoints, double lambda) {
	int ii;
	double result = 0;
	double num = 0;
	
	for (ii = 0; ii < datapoints ; ii++) {
		num = derivativeMatrix[row][ii]	* derivativeMatrix[col][ii];
		num /= edata[ii] * edata[ii];
		
		result += num;
	}
	if(row == col)
		result *= 1 + lambda;
	
	alpha[row][col] = result;
}


/** packs the upper right elements of the alpha matrix, because the alpha matrix should be symmetrical*/
void packAlphaSymmetric(double** alpha, unsigned int numvarparams){  
	unsigned int ii,jj;
	
	for(ii = 0 ; ii < numvarparams ; ii++)
		for(jj = numvarparams - 1 ; jj > ii ; jj--)
			alpha[ii][jj] = alpha[jj][ii];	
}

/** Calculates the lower left elements for <code>alpha</code>. */
void updateAlpha(double **alpha, double **derivativeMatrix,  unsigned int numvarparams, const double *edata, long datapoints, double lambda) {
	unsigned int ii, jj;
	for (ii = 0; ii < numvarparams; ii++) {
		for (jj = 0; jj < ii+1 ; jj++)
			calculateAlphaElement(ii, jj, alpha, derivativeMatrix, edata, datapoints, lambda);
	}
	packAlphaSymmetric(alpha, numvarparams);
}

/** 
 * @return An calculated element for the beta-matrix.
 * NOTE: Does not change the value of beta-matrix.
 */
double calculateBetaElement(double **derivativeMatrix, int row, const double *ydata, const double *model, const double *edata, long datapoints) {
	int ii;
	double result = 0;
	for (ii = 0 ; ii < datapoints ; ii++){
		result +=  (ydata[ii] - model[ii]) * derivativeMatrix[row][ii];
		result /= edata[ii] * edata[ii];
	}

	return result;
}

/** Calculates all elements for <code>beta</code>. */
void updateBeta(double *b, double **derivativeMatrix, int numvarparams, const double *ydata, const double *model, const double *edata, long datapoints) {
	int ii;
	for (ii = 0; ii < numvarparams; ii++)
		b[ii] = calculateBetaElement(derivativeMatrix, ii, ydata, model, edata, datapoints);
}


int getCovarianceMatrix(double ***covarianceMatrix,
						double *hessianDeterminant,
						void *userdata,
						fitfunction fitfun,
						double cost,
						double *coefs,
						int numcoefs,
						unsigned int *holdvector,
						const double *ydata,
						const double *edata,
						const double **xdata,
						long datapoints,
						int numDataDims,
						int unitSD){
	int err;
	double **derivativeMatrix = NULL;
	double **reducedCovarianceMatrix = NULL;
	double hess = 0;
	unsigned int *varparams = NULL;
	int ii,jj, numvarparams = 0;
	err = 0;
	
	//fit function must exist
	if(!fitfun)
		return NO_FIT_FUNCTION_SPECIFIED;
	
	for(ii = 0 ; ii < numcoefs ; ii++)
		if(holdvector[ii] == 0)
			numvarparams++;

	varparams = (unsigned int*) malloc (sizeof(unsigned int) * numvarparams);
	if(!varparams){
		err = NO_MEMORY;
		goto done;
	}

	jj = 0;
	for(ii = 0 ; ii < numcoefs ; ii++)
		if(holdvector[ii] == 0){
			varparams[jj] = ii;
			jj++;
		}

	reducedCovarianceMatrix = (double**) malloc2d(numvarparams, numvarparams, sizeof(double));
	if(reducedCovarianceMatrix == NULL){
		err = NO_MEMORY;
		goto done;
	}

	derivativeMatrix = (double**) malloc2d(numvarparams, datapoints, sizeof(double));
	if(derivativeMatrix == NULL){
		err = NO_MEMORY;
		goto done;
	}
	
	if((err = updatePartialDerivative(userdata, fitfun, derivativeMatrix, coefs, numcoefs, varparams, numvarparams, xdata, datapoints, numDataDims)))
	   goto done;
	   	
	updateAlpha(reducedCovarianceMatrix, derivativeMatrix, numvarparams, edata, datapoints, 0);

	if((err = matrixInversion_chol(reducedCovarianceMatrix, numvarparams, &hess)))
		goto done;
		
	if(unitSD)
		for(ii = 0; ii < numvarparams ; ii++)
			for(jj = 0 ; jj < numvarparams ; jj += 1)
				reducedCovarianceMatrix[ii][jj] *= cost/(datapoints - numvarparams);
	
	*covarianceMatrix = (double**) malloc2d(numcoefs, numcoefs, sizeof(double));
	if(!covarianceMatrix){
		err = NO_MEMORY;
		goto done;
	}
	
	for (ii = 0; ii < numvarparams; ii++)
		for(jj = 0 ; jj < numvarparams ; jj++)
			(*covarianceMatrix)[varparams[ii]][varparams[jj]] =  reducedCovarianceMatrix[ii][jj];
	
	if(hessianDeterminant)
		*hessianDeterminant = hess;
			
done:
	if(varparams != NULL)
		free(varparams);
	if(reducedCovarianceMatrix != NULL)
		free(reducedCovarianceMatrix);
	if(derivativeMatrix != NULL)
		free(derivativeMatrix);
	
	return err;
}

/*
 insertVaryingParams inserts the current pvector into an array copy of the coefficients,
 then into a temporary wave
 returns 0 if no error
 returns errorcode otherwise
 */
static void
insertVaryingParams(double *coefs, const unsigned int* varparams, unsigned int numvarparams, double *vector){
	unsigned int ii;
	
	for(ii = 0 ; ii < numvarparams ; ii += 1)
		*(coefs + *(varparams + ii)) =  *(vector + ii);
		
}

/**
 does a levenberg marquardt fit to the data, instead of differential evolution.  It returns a
 non-zero error code if something goes wrong.  However, it will also stop if your fitfunction 
 returns a non-zero value.  As with genetic optimisation you can supply your own cost function.
 
 @param fitfun					- a function that calculates the dependent variable, given input parameters and independent variables. 
 If you return a non-zero value from this function the fit will stop. 
 
 @param costfun					- a function that calculates the costfunction to be minimised.  This is normally a chi2 type function.
 i.e. sum (((model[i] - data[i]) / dataerrors[i])^2 )
 If costfun == NULL then a default chi2 function is used.
 
 @param numcoefs				- total number of fit parameters.
 
 @param coefs[numcoefs]			- an array containing all the parameters for the fit.  After genetic_optimisation this is populated by the parameters
 that best fit the data.
 
 @param holdvector[numcoefs]	- an array (with numcoefs elements) that specifies which parameters are going to be held during the fit. 
 0 = vary
 1 = hold
 
 @param datapoints				- the total number of data points in the fit.
 
 @param ydata[datapoints]		- an array containing the dependent variable (i.e. the data one is trying to fit).
 
 @param xdata[numDataDims][datapoints]  - a 2D array containing the independent variables that correspond to each of the datapoints.
 One can fit multidimensional data, e.g. y = f(n, m).  In this case numDataDims = 2.
 You can allocate a 2D dataset with m points using malloc2D(2, m, sizeof(double)).
 If you want to pass in a 1D dataset simply pass a pointer to the array.
 e.g. if your array is:
 double *xP;
 then pass in:
 &xP
 BUT YOU HAVE TO REMEMBER TO DEREFERENCE THE POINTER IN THE FIT FUNCTION BEFORE YOU USE THE ARRAY.
 model[ii] = (*xP)[ii]
 
 @param edata[datapoints]		- an array containing the experimental uncertainties for each of the datapoints.  If you use the default chi2 costfunction
 then it should contain standard deviations.  Set each element to 1 if you do not wish to weight the fit by the experimental
 uncertainties.  
 
 @param numDataDims				- the number of independent variables in the fit. For y = f(x) numDataDims = 1.  For y = f(n, m), numDataDims = 2, etc.
 
 @param chi2					- the final value of the cost function.
 
 @param gco						- options for the genetic optimisation.  (see above).  If gco == NULL, then a default set of options are used.
 
 @param userdata				- an (optional) pointer that is passed to the fitfunction, costfunction and updatefunction.  Use this pointer to give extra
 information to your functions.
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
					 void* userdata){

	int err = 0;
	unsigned int ii, jj, numvarparams = 0, iterations;
	unsigned int *varparams = NULL;
	double cost = -1, incrementedCost = -1, lambda = 0.001;
	double **derivativeMatrix = NULL;
	double **alpha = NULL;
	double *reducedParameters = NULL;
	double *beta = NULL;
	double *model = NULL;
	double *temp_coefs = NULL;
	double *incrementedParameters = NULL;
	gencurvefitOptions lgco;
	costfunction mycostfun;
	
	//fit function must exist
	if(!fitfun)
		return NO_FIT_FUNCTION_SPECIFIED;
	
	//y, x, e arrays must exist
	if(!ydata)
		return NO_Y_ARRAY;
	
	if(!xdata)
		return NO_X_ARRAY;
	
	if(!edata)
		return NO_E_ARRAY;
	
	if(!coefs)
		return NO_COEFS_ARRAY;
	
	if(gco == NULL){
		lgco.iterations = 100;
		lgco.tolerance = 1.0e-8;
	} else {
		lgco.iterations = gco->iterations;
		lgco.tolerance = gco->tolerance;
	}
	
	if(costfun == NULL)
		mycostfun = &chisquared;
	else
		mycostfun = costfun;

		
	for(ii = 0; ii < numcoefs ; ii++)
		if(holdvector[ii] == 0) numvarparams ++;

	varparams = (unsigned int*) malloc(sizeof(unsigned int) * numvarparams);
	if(!varparams){
		err = NO_MEMORY;
		goto done;
	}
	
	derivativeMatrix = (double**) malloc2d(numvarparams, datapoints, sizeof(double));
	if(!derivativeMatrix){
		err = NO_MEMORY;
		goto done;
	}
	
	alpha = (double**) malloc2d(numvarparams, numvarparams, sizeof(double));
	if(!alpha){
		err = NO_MEMORY;
		goto done;
	}
	
	reducedParameters = (double*) malloc(sizeof(double) * numvarparams);
	if(!reducedParameters){
		err = NO_MEMORY;
		goto done;
	}
	
	beta = (double*) malloc(sizeof(double) * numvarparams);
	if(!beta){
		err = NO_MEMORY;
		goto done;
	}
	
	incrementedParameters = (double*) malloc(sizeof(double) * numvarparams);
	if(!incrementedParameters){
		err = NO_MEMORY;
		goto done;
	}
	
	temp_coefs = (double*) malloc(sizeof(double) * numcoefs);
	if(!temp_coefs){
		err = NO_MEMORY;
		goto done;
	}
	memcpy(temp_coefs, coefs, sizeof(double) * numcoefs);
	
	jj = 0;
	for(ii = 0 ; ii < numcoefs ; ii++){
		if(holdvector[ii] == 0){
			varparams[jj] = ii;
			reducedParameters[jj] = (double) coefs[ii];
			jj++;
		}
	}
	
	model = (double*) malloc(sizeof(double) * datapoints);
	if(!model){
		err = NO_MEMORY;
		goto done;
	}
	
	iterations = 0;
	do {
		insertVaryingParams(temp_coefs, varparams, numvarparams, reducedParameters);
		
		if((err = fitfun(userdata, temp_coefs, numcoefs, model, xdata, datapoints, numDataDims)))
			goto done;
		

		cost = mycostfun(userdata, temp_coefs, numcoefs, model, ydata, edata, datapoints);
		
		if((err = updatePartialDerivative(userdata, fitfun, derivativeMatrix, temp_coefs, numcoefs, varparams, numvarparams, xdata, datapoints, numDataDims)))
			goto done;
		
		updateAlpha(alpha, derivativeMatrix, numvarparams, edata, datapoints, lambda);
		   
		updateBeta(beta, derivativeMatrix, numvarparams, ydata, model, edata, datapoints);
		
		if((err = matrixInversion_lu(alpha, numvarparams)))
			goto done;
		
		for (ii = 0; ii < numvarparams ; ii++){
			double val = 0;
			for(jj = 0 ; jj < numvarparams ; jj++)
				val += alpha[ii][jj] * beta[jj];
			
			incrementedParameters[ii] = reducedParameters[ii] + val;
		}
		
		insertVaryingParams(temp_coefs, varparams, numvarparams, incrementedParameters);

		if((err = fitfun(userdata, temp_coefs, numcoefs, model, xdata, datapoints, numDataDims)))
			goto done;
		
		incrementedCost = mycostfun(userdata, temp_coefs, numcoefs, model, ydata, edata, datapoints);
		
		// The guess results to worse chi2 - make the step smaller
		if (incrementedCost >= cost) {
			lambda *= 10;
		}
		// The guess results to better chi2 - move and make the step larger
		else {
			lambda /= 10;
			memcpy(reducedParameters, incrementedParameters, sizeof(double) * numvarparams);
			insertVaryingParams(coefs, varparams, numvarparams, reducedParameters);
		}
		iterations++;
	} while (iterations < lgco.iterations && lgco.tolerance < fabs(incrementedCost - cost));
	
	if(chi2)
		*chi2 = cost;
	
done:
	if(temp_coefs)
		free(temp_coefs);
	if(alpha)
		free(alpha);
	if(model)
		free(model);
	if(derivativeMatrix)
		free(derivativeMatrix);
	if(beta)
		free(beta);
	if(reducedParameters)
		free(reducedParameters);
	if(incrementedParameters)
		free(incrementedParameters);
	if(varparams)
		free(varparams);
	
	return 0;
}

