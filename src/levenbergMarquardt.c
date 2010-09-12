
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

//Cholesky Decomposition
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
	
	
done:
	return err;
}


//Cholesky back substitution
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


static int matrixInversion(double **a, int N, double *detA){
	int err=0;
	int i,j;
	double *x = NULL;
	double *b = NULL;
	double *p = NULL;
	double **tempA = NULL;
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
	if(err = choldc(tempA, N, p)) goto done;
	
	//now do the back substitution
	for(j = 0 ; j < N ; j++){
		memset(b, 0, sizeof(double) * N);
		b[j] = 1.0;
		memset(x, 0, sizeof(double) * N);
		cholsl(tempA, N, p, b, x);
		
		for(i = 0 ; i < N ; i++)
			a[i][j] = x[i];
	}
	
	//make the covariance matrix symmetric
	for(i=0 ; i < N ; i++)
		for(j = N-1 ; j > i ; j--)
			a[i][j] = a[j][i];
	
	
	//the determinant of the original matrix is the square of the products of the elements in the cholesky diagonal
	for(i = 0 ; i < N ; i+=1)
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

static int partialDerivative(void *userdata, fitfunction fitfun, double** derivativeMatrix, int derivativeMatrixRow, int parameterIndex, double* coefs, int numcoefs, double **xdata, long datapoints, int numDataDims){
	int err = 0;
	double param, diff;
	int jj;
	double *dataTemp = NULL;
	
	param = coefs[parameterIndex];	
	diff = 1.e-6 * param;
	coefs[parameterIndex] = param + diff;
	
	if(err = fitfun(userdata, coefs, numcoefs, *(derivativeMatrix + derivativeMatrixRow), (const double**)xdata, datapoints, numDataDims))
		goto done;
	
	coefs[parameterIndex] = param - diff;
	
	dataTemp = (double*)malloc(sizeof(double) * datapoints);
	if(!dataTemp){
		err = NO_MEMORY;
		goto done;
	}

	if(err = fitfun(userdata, coefs, numcoefs, dataTemp, (const double**)xdata, datapoints, numDataDims))
		goto done;
	
	for(jj = 0 ; jj < datapoints ; jj++)
		derivativeMatrix[derivativeMatrixRow][jj] = (derivativeMatrix[derivativeMatrixRow][jj] - dataTemp[jj]) / (2 * diff);
	
	coefs[parameterIndex] = param;	

done:
	if(dataTemp)
		free(dataTemp);
	
	return err;
}


static int updatePartialDerivative(void *userdata, fitfunction fitfun, double **derivativeMatrix, double *coefs, int numcoefs, unsigned int *varparams, int numvarparams, double **xdata, long datapoints, int numDataDims){
	int err = 0;
	int ii;
	for(ii = 0 ; ii < numvarparams ; ii++){
		if(err = partialDerivative(userdata, fitfun, derivativeMatrix, ii, varparams[ii], coefs, numcoefs, xdata, datapoints, numDataDims))
			return err;
	} 
	return err;
}


static int calculateAlphaElement(int row, int col, double **alpha, double **derivativeMatrix, const double *edata, long datapoints, double lambda) {
	int err = 0;
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
	return err;
}


/** packs the upper right elements of the alpha matrix, because the alpha matrix should be symmetrical*/
static int packAlphaSymmetric(double** alpha, unsigned int numvarparams){  
	int err = 0,ii,jj;
	
	for(ii = 0 ; ii < numvarparams ; ii++)
		for(jj = numvarparams - 1 ; jj > ii ; jj--)
			alpha[ii][jj] = alpha[jj][ii];
	
	return err;
}

/** Calculates the lower left elements for <code>this.alpha</code>. */
static int updateAlpha(double **alpha, double **derivativeMatrix,  unsigned int numvarparams, const double *edata, long datapoints, double lambda) {
	int err = 0, ii, jj;
	for (ii = 0; ii < numvarparams; ii++) {
		for (jj = 0; jj < ii+1 ; jj++) {
			if(err = calculateAlphaElement(ii, jj, alpha, derivativeMatrix, edata, datapoints, lambda))
				return err;
		}
	}
	if(err = packAlphaSymmetric(alpha, numvarparams))
		return err;
	return err;
}

/** 
 * @return An calculated element for the beta-matrix.
 * NOTE: Does not change the value of beta-matrix.
 */
double calculateBetaElement(int row, const double *edata, long datapoints) {
	int ii;
	double result = 0;
	for (ii = 0 ; ii < datapoints ; ii++) {
		result +=  edata[ii]; //* 
//		(dataPoints[1][i] - function.getY(dataPoints[0][i], parameters)) *
//		function.getPartialDerivate(dataPoints[0][i], parameters, row);
	}
	return result;
}

/** Calculates all elements for <code>this.beta</code>. */
void updateBeta(double *b, int numvarparams, const double *edata, long datapoints) {
	int ii;
	for (ii = 0; ii < numvarparams; ii++)
		b[ii] = calculateBetaElement(ii, edata, datapoints);
}


int getCovarianceMatrix(double **covarianceMatrix,
						void *userdata,
						fitfunction fitfun,
						double cost,
						double *coefs,
						int numcoefs,
						unsigned int *holdvector,
						double *ydata,
						double *edata,
						long datapoints,
						double **xdata,
						int numDataDims,
						int unitSD){
	int err;
	double **derivativeMatrix = NULL;
	double **reducedCovarianceMatrix = NULL;
	double hessianDeterminant = 0;
	double val = 0;
	unsigned int *varparams = NULL;
	int ii,jj, numvarparams;
	err = 0;
	
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
	
	if(err = updatePartialDerivative(userdata, fitfun, derivativeMatrix, coefs, numcoefs, varparams, numvarparams, xdata, datapoints, numDataDims))
	   goto done;
	   	
	if(err = updateAlpha(reducedCovarianceMatrix, derivativeMatrix, numvarparams, edata, datapoints, 0))
		goto done;

	if(err = matrixInversion(reducedCovarianceMatrix, numvarparams, &hessianDeterminant)) goto done;
		
	if(unitSD)
		for(ii = 0; ii < numvarparams ; ii++)
			for(jj = 0 ; jj < numvarparams ; jj += 1)
				reducedCovarianceMatrix[ii][jj] *= cost/(datapoints - numvarparams);
	
	covarianceMatrix = (double**) malloc2d(numcoefs, numcoefs, sizeof(double));
	if(!covarianceMatrix){
		err = NO_MEMORY;
		goto done;
	}
	
	for(ii = 0 ; ii < numcoefs ; ii++)
		for(jj = 0 ; jj < numcoefs ; jj++)
			covarianceMatrix[ii][jj] = sqrt(-1);
	
	for (ii = 0; ii < numvarparams; ii++) {
		for(jj = 0 ; jj < numvarparams ; jj++){
			val = reducedCovarianceMatrix[ii][jj];
			covarianceMatrix[varparams[ii]][varparams[jj]];
		}
	}
		
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
static int
insertVaryingParams(double *coefs, const unsigned int* varparams, unsigned int numvarparams, double *vector){
	int err=0, ii;
	
	for(ii = 0 ; ii < numvarparams ; ii += 1)
		*(coefs + *(varparams + ii)) =  *(vector + ii);
	
	return err;
}

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
	int ii, jj, numvarcoefs = 0, iterations;
	unsigned int *varcoefs = NULL;
	double cost = -1, incrementedCost = -1, lambda;
	double **derivativeMatrix = NULL;
	double **alpha = NULL;
	double *reducedParameters = NULL;
	double *da = NULL;
	double *b = NULL;
	double *model = NULL;
	double *temp_coefs = NULL;
	double *incrementedParameters = NULL;
	gencurvefitOptions lgco;
	costfunction mycostfun;
	
	
	if(gco == NULL){
		lgco.iterations = 100;
		lgco.tolerance = 0.001;
	} else {
		lgco.iterations = gco->iterations;
		lgco.tolerance = gco->iterations;
	}
	
	if(costfun == NULL)
		mycostfun = &chisquared;
	else
		mycostfun = costfun;

		
	for(ii = 0; ii < numcoefs ; ii++)
		if(holdvector[ii] == 0) numvarcoefs ++;

	varcoefs = (unsigned int*) malloc(sizeof(unsigned int) * numvarcoefs);
	if(!varcoefs){
		err = NO_MEMORY;
		goto done;
	}
	
	derivativeMatrix = (double**) malloc2d(numvarcoefs, numvarcoefs, sizeof(double));
	if(!derivativeMatrix){
		err = NO_MEMORY;
		goto done;
	}
	
	alpha = (double**) malloc2d(numvarcoefs, numvarcoefs, sizeof(double));
	if(!alpha){
		err = NO_MEMORY;
		goto done;
	}
	
	reducedParameters = (double*) malloc(sizeof(double) * numvarcoefs);
	if(!reducedParameters){
		err = NO_MEMORY;
		goto done;
	}
	
	da = (double*) malloc(sizeof(double) * numvarcoefs);
	if(!da){
		err = NO_MEMORY;
		goto done;
	}
	
	b = (double*) malloc(sizeof(double) * numvarcoefs);
	if(!b){
		err = NO_MEMORY;
		goto done;
	}
	
	incrementedParameters = (double*) malloc(sizeof(double) * numvarcoefs);
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
	for(ii = 0 ; ii < numvarcoefs ; ii++){
		if(holdvector[ii] == 0){
			varcoefs[jj] = ii;
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
		insertVaryingParams(temp_coefs, varcoefs, numvarcoefs, reducedParameters);
		
		if(err = fitfun(userdata, temp_coefs, numcoefs, model, xdata, datapoints, numDataDims))
			goto done;
		
		cost = mycostfun(userdata, temp_coefs, numcoefs, model, ydata, edata, datapoints);
		
		if(err = updateAlpha(alpha, derivativeMatrix, numvarcoefs, edata, datapoints, 0))
		   goto done;
		   
		updateBeta(b, numvarcoefs, edata, datapoints);
		
//		solveIncrements();
		
		insertVaryingParams(temp_coefs, varcoefs, numvarcoefs, incrementedParameters);
		incrementedCost = mycostfun(userdata, temp_coefs, numcoefs, model, ydata, edata, datapoints);
		
		// The guess results to worse chi2 - make the step smaller
		if (incrementedCost >= cost) {
			lambda *= 10;
		}
		// The guess results to better chi2 - move and make the step larger
		else {
			lambda /= 10;
			memcpy(reducedParameters, incrementedParameters, sizeof(double) * numvarcoefs);
			insertVaryingParams(coefs, varcoefs, numvarcoefs, reducedParameters);
		}
		iterations++;
	} while ( iterations < lgco.iterations || lgco.tolerance > fabs(cost - incrementedCost));
	
	
done:
	if(temp_coefs)
		free(temp_coefs);
	if(alpha)
		free(alpha);
	if(model)
		free(model);
	if(derivativeMatrix)
		free(derivativeMatrix);
	if(da)
		free(da);
	if(b)
		free(b);
	if(reducedParameters)
		free(reducedParameters);
	if(incrementedParameters)
		free(incrementedParameters);
	
	if(varcoefs)
		free(varcoefs);
	
	return 0;
}

