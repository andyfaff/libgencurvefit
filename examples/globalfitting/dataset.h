/*
 *  dataset.h
 *  motoMC
 *
 *  Created by Andrew Nelson on 1/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include <vector>
#include <algorithm>
#include <fstream>
#include <string>
#include "gencurvefit.h"

using namespace std;

class dataset {
public:
	vector<double> xx;
	vector<double> yy;
	vector<double> dy;
	vector<double> dx;
	long datapoints;
	unsigned numDataDims;
	const char* filename;
	
	dataset();
	~dataset();
	
	int readDataFile(const char* filename);	
	
};

void Tokenize(const string &str, vector<string> &tokens, const char* DELIMITERS, int szDELIMITERS);

int readCoefficientFile(const char* filename,
						vector <double> *coefficients,
						vector <unsigned int> *holdvector,
						vector <double> *lowlim,
						vector <double> *hilim,
						string *fitfunctionStr,
						string *costfunctionStr);


