/*
 *  myfitfunctions.h
 *  motoMC
 *
 *  Created by andrew on 29/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#define PI 3.14159265358979323846
#include <vector>
#define LINKAGE_MATRIX_WRONG 1

using namespace std;

//fitfunctions
int Abeles(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims);
int line(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims);
int gaussian(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims);
int AbelesModelWrapper(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims);
int globalFitWrapper(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims);


//costfunctions
double smoother(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts);
double log10ChiSquared(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts);

