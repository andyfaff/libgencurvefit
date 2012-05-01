/*
 *  untitled.cpp
 *  codevalidation
 *
 *  Created by Andrew Nelson on 1/05/12.
 *  Copyright 2012 __MyCompanyName__. All rights reserved.
 *
 */

#include "gencurvefit.h"
#include <vector>
#include <iostream>
#include "myfitfunctions.h"
#include "dataset.h"

using namespace std;

int main (int argc, char *argv[]) {
	int err = 0;
	int ii;
	vector<double> yy;
	vector<double> coefficients;
	double *qp;
	double chi2;
	
	dataset ds;
	ds.readDataFile(argv[1]);
	
	yy.resize(ds.xx.size());
	
	coefficients.push_back(2);
	coefficients.push_back(1);	
	coefficients.push_back(2.07);
	coefficients.push_back(6.36);	
	coefficients.push_back(1.0e-5);
	coefficients.push_back(2);	
	coefficients.push_back(10);
	coefficients.push_back(3.47);	
	coefficients.push_back(0);
	coefficients.push_back(2);	
	coefficients.push_back(211);	
	coefficients.push_back(0.5);
	coefficients.push_back(0);
	coefficients.push_back(2);	
	
	qp = &(ds.xx[0]);
	
	abeles(NULL,
		   (const double*) &coefficients[0],
		   coefficients.size(),
		   &yy[0],
		   (const double **) &qp,
		   ds.xx.size(),
		   1);

	chi2 = log10chisquared(NULL, NULL, 0, (const double*) &(ds.yy[0]), (const double*) &(yy[0]), (const double*) &(ds.dy[0]), ds.xx.size());
	
	
	for(ii = 0 ; ii < ds.xx.size() ; ii++)
		cout << scientific << ds.xx[ii] << " " << yy[ii] << endl;
	
    cout << chi2 << endl;
	return err;
}

