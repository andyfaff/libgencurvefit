/*
 *  dataset.cpp
 *  motoMC
 *
 *  Created by Andrew Nelson on 1/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */

#include "dataset.h"
#include <vector>
#include <fstream>
#include <math.h>
#include "myfitfunctions.h"
#include "gencurvefit.h"
#include "iostream"

void Tokenize(const string &str, vector<string> &tokens, const char* DELIMITERS, int szDELIMITERS){
    // Skip delimiters at beginning.
	string delimiters(DELIMITERS,szDELIMITERS);		//WARNING, DELIMITERS ISNT NULL TERMINATED
	
	string::size_type lastPos = str.find_first_not_of(delimiters, 0);
	// Find first "non-delimiter".
	string::size_type pos     = str.find_first_of(delimiters, lastPos);
	
	while (string::npos != pos || string::npos != lastPos)
	{
		// Found a token, add it to the vector.
		tokens.push_back(str.substr(lastPos, pos - lastPos));
		// Skip delimiters.  Note the "not_of"
		lastPos = str.find_first_not_of(delimiters, pos);
		// Find next "non-delimiter"
		pos = str.find_first_of(delimiters, lastPos);
	}
}

dataset::dataset(){
	xx.clear();
	yy.clear();
	dy.clear();
	dx.clear();
	filename = NULL;
	datapoints = 0;
	numDataDims = 0;
};

dataset::~dataset(){
}

int dataset::readDataFile(const char* filename){
		int err = 0;
		ifstream file_to_read;
		vector<string> columndata;
		string linein;
		
		xx.clear();
		yy.clear();
		dy.clear();
		dx.clear();
		datapoints = 0;
		numDataDims = 1;
		filename = filename;
		
		//read the coefficient file, 1st two lines are headers.
		file_to_read.open(filename, ios::in);
		
		if(!file_to_read)
			return 1;
		
		columndata.clear();
		
		//read the data
		while(getline(file_to_read, linein, '\n')){
			Tokenize(linein, columndata, "\t", sizeof(char));			
			xx.push_back(strtod(columndata[0].c_str(), NULL));
			yy.push_back(strtod(columndata[1].c_str(), NULL));
			dy.push_back(fabs(strtod(columndata[2].c_str(), NULL)));
			columndata.clear();
		}
		datapoints = yy.size();

done:
		
		if(file_to_read.is_open())
			file_to_read.close();
		
		return err;
		
	};

/*
 read a coefficient file, parse the coefficients and holdvector and limits.
 return 0 if no error. return a number if an error.
 
 filename			- a c string containing the filename to be read
 
 coefficients		- a vector holding the coefficients to be fitted
 
 holdvector			- a vector that indicates whether each parameter will be fitted (0) or held (1)
 
 lowlim				- a vector containing the lower limits
 
 hilim				- a vector containing the upper limits
 
 The file format is 2 lines followed by lines containing the coefficients:
 Each of the coefficient lines is:
 number toHold lowerLim upperLim
 where number is a coefficient, toHold is either 0 or 1.  (see above) and lowerlim and upperlim
 are the limits for that parameter.  All numbers are SPACE delimited.  The end of line character is \n.
 
 e.g. 
 
 stuff1\n
 stuff2\n
 1.023 1 0.9 1.1\n
 2.07 0 2 3\n
 */

int readCoefficientFile(const char* filename, vector <double> &coefficients, vector <unsigned int> &holdvector, vector <double> &lowlim, vector <double> &hilim){
	int err = 0;
	
	ifstream file_to_read;
	vector<string> columndata;
	string linein;
	
	//read the coefficient file, 1st two lines are headers.
	file_to_read.open(filename, ios::in);
	
	if(!file_to_read)
		return 1;
	
	columndata.clear();
	
	getline(file_to_read, linein, '\n');		//chi2value header
	getline(file_to_read, linein, '\n');		//header describing columns
	
	//now read each line
	while(getline(file_to_read, linein, '\n')){
		Tokenize(linein, columndata, " ", sizeof(char));
		//		std::cout << columndata[0] << "\t" << columndata[1] << "\t" << columndata[2] << "\r";
		
		coefficients.push_back(strtod(columndata[0].c_str(), NULL));
		holdvector.push_back(strtol(columndata[1].c_str(), NULL, 10));			
		lowlim.push_back(strtod(columndata[2].c_str(), NULL));
		hilim.push_back(strtod(columndata[3].c_str(), NULL));
		
		columndata.clear();
	}
	
done:
	if(file_to_read.is_open())
		file_to_read.close();
	
	return err;
	
}



