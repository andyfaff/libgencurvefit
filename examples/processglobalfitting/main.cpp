
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <iomanip>
#include "math.h"

#include "dataset.h"
#include "globalfitfunction.h"

using namespace std;

std::string to_a_string(double *nums, long numthings)
{
	std::ostringstream oss;
	long ii;
	for(ii=0 ; ii< numthings ; ii+=1)
		oss << nums[ii] << " ";
	
	return oss.str();
}	

int expandGlobalLinkageTable(globalFit &gFS, string &coefsToExpand, vector<double> *expandedCoefs){
	int err = 0;
	int ii = 0, jj = 0;
	vector<string> theNumbers;
	int valpos;
	double val;
	for(ii = 0 ; ii< gFS.numDataSets ; ii++)
		expandedCoefs[ii].clear();
	
	Tokenize(coefsToExpand, theNumbers, " ", sizeof(char));
	
	//first element is chi2
	theNumbers.erase(theNumbers.begin());
	
	for(ii = 0 ; ii <gFS.numDataSets ; ii++){
		for(jj = 0 ; jj < gFS.globalFitIndividualArray[ii].numcoefs ; jj++){
			valpos = gFS.parameterLinkageTable[ii].at(jj);
			val = strtod(theNumbers.at(valpos).c_str(), NULL) ;
			expandedCoefs[ii].push_back(val);
		}
	}
	return err;
}

int main (int argc, char *argv[]) {
	int err = 0;
	
	globalFit gFS;
	dataset theDataSet;
	ifstream file_to_read;

	vector<double> coefs;
	vector<double> lowlim;
	vector<double> hilim;
	double val;
	
	vector<unsigned int> bs;
	string coefsToExpand;
	string poop;
	int ii, jj, numiterations = 0;
	vector<double> *expandedCoefs = NULL;
	vector<double> *mean = NULL;
	vector<double> *standard_deviation = NULL;	
	
	ofstream *outputfiles = NULL;
	
	if(argc != 3){
		cout << "Useage:\n ./motoMC globalpilot iterations\n";
		err = 1;
	}

	//set up global fit structure
	if(err = parseGlobalPilotFile(argv[1], gFS, theDataSet, coefs, bs, lowlim, hilim))
		goto done;
	
	file_to_read.open(argv[2], ios::in);
	
	expandedCoefs = new (nothrow) vector<double>[gFS.numDataSets];
	if(!expandedCoefs)
		return 1;
	
	mean = new (nothrow) vector<double>[gFS.numDataSets];
	if(!mean)
		return 1;
	
	standard_deviation = new (nothrow) vector<double>[gFS.numDataSets];
	if(!standard_deviation)
		return 1;
	
	outputfiles = new (nothrow) ofstream[gFS.numDataSets];
	if(!outputfiles)
		return 1;
	
	for(ii = 0 ; ii < gFS.numDataSets ; ii++)
		outputfiles[ii].open((gFS.globalFitIndividualArray[ii].datafilename + ".dat").c_str());
	
	
	while(getline(file_to_read, coefsToExpand, '\n')){
		if(err = expandGlobalLinkageTable(gFS, coefsToExpand, expandedCoefs))
			goto done;
				
		for(ii = 0 ; ii < gFS.numDataSets ; ii++){
			outputfiles[ii] << to_a_string(&((expandedCoefs[ii])[0]), expandedCoefs[ii].size()) << endl ;
			mean[ii].resize(expandedCoefs[ii].size());
			standard_deviation[ii].resize(expandedCoefs[ii].size());
			
			for(jj = 0 ; jj < expandedCoefs[ii].size() ; jj ++){
				mean[ii][jj] += expandedCoefs[ii][jj];
				standard_deviation[ii][jj] += pow(expandedCoefs[ii][jj], 2);
			}
		}
		numiterations ++;		
	}
	
	
	
	for(ii = 0 ; ii < gFS.numDataSets ; ii++){
		cout << "dataset" << ii << endl;
		cout << "mean" << endl;
		for(jj = 0 ; jj < mean[ii].size() ; jj++){
			mean[ii][jj]/=numiterations;
			cout << setw(13) << scientific << mean[ii][jj] << " ";
		}
		cout << endl;
		cout << "standard deviation" << endl;
		for(jj = 0 ; jj < mean[ii].size() ; jj++){
			val =  sqrt((standard_deviation[ii][jj]/numiterations) - pow(mean[ii][jj], 2));
			if(!isfinite(val))
				val = 0;
			cout << setw(13) << scientific << val << " ";
		}
		cout << endl;
		cout << endl;
	}
	
	for(ii = 0 ; ii<gFS.numDataSets ; ii++)
		outputfiles[ii].close();
	
done:
	if(expandedCoefs)
		delete [] expandedCoefs;
	if(mean)
		delete [] mean;
	if(standard_deviation)
		delete [] standard_deviation;
	
	if(outputfiles)
		delete [] outputfiles;
    return err;
}



