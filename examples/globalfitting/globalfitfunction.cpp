/*
 *  globalfitfunction.cpp
 *  motoMC
 *
 *  Created by Andrew Nelson on 1/09/10.
 *  Copyright 2010 __MyCompanyName__. All rights reserved.
 *
 */
#include "myfitfunctions.h"
#include "dataset.h"
#include "globalfitfunction.h"

#include <iostream>

using namespace std;

globalFit::globalFit(){
	 numDataSets = 0;
	parameterLinkageTable = NULL;	
	globalFitIndividualArray = NULL;
};

globalFit::~globalFit(){
	if(parameterLinkageTable)
		delete[] parameterLinkageTable;
	if(globalFitIndividualArray)
		delete[] globalFitIndividualArray;
};

int globalFit::makedatasets(int dataSetsToCreate){
	int err = 0;
	
	globalFitIndividualArray = new (nothrow) globalFitIndividual[dataSetsToCreate];
	if(!globalFitIndividualArray)
		return NO_MEMORY;
	
	parameterLinkageTable = new (nothrow) vector<int>[dataSetsToCreate];
	if(!parameterLinkageTable)
		return NO_MEMORY;
	numDataSets = dataSetsToCreate;
	
	return err;
}



globalFitIndividual::globalFitIndividual(){
	ffp = NULL;
	numcoefs = 0;
	datapoints = 0;
	numDataDims = 0;
	datafilename.clear();
	pilotfilename.clear();
}

/*
The following function checks the linkage array to check that it's the same size as the local coefficients array.
 Also that the linkage array is sane.
 Valid    Invalid
 0 5		0	8
 1 6		2	9
 2 7		3	10
 3			4
 3			5
 2			6
 4			7
 
 In other words, subsequent entries can be less than or equal to the highest parameter number so far PLUS one,
 but they can't be negative (the initial highest parameter number is -1).
 */
int sanityCheckOnLinkageValues(globalFit &gFS, vector<double> *localcoefficientsArray){
	int err = 0;
	int ii, jj, totalNumberOfParams = 0, highestParamSoFar = -1;
	
	
	for(ii = 0 ; ii < gFS.numDataSets ; ii += 1){
		if(gFS.parameterLinkageTable[ii].size() != localcoefficientsArray[ii].size()){
			cout << "Linkage array does not match coefficients array\n";
			return 1;
		}
		totalNumberOfParams += localcoefficientsArray[ii].size();
	}
	for(ii = 0 ; ii < gFS.numDataSets ; ii++){
		for(jj = 0 ; jj < gFS.parameterLinkageTable[ii].size() ; jj++){
			if(gFS.parameterLinkageTable[ii].at(jj) > highestParamSoFar + 1){
				cout << "Linkage Matrix is wrong\n";
				return 1;
			} else if (gFS.parameterLinkageTable[ii].at(jj) == highestParamSoFar + 1){
				highestParamSoFar += 1;
			} else if(gFS.parameterLinkageTable[ii].at(jj) < 0){
				cout << "Negative value in linkage array\n";
				return 1;
			}
		}
	}
	if(highestParamSoFar > totalNumberOfParams - 1){
		cout << "ERROR you had totalNumberOfParams being exceeded in the linkagematrix\n";
		return 1;
	}
	return err;
}

int parseGlobalPilotFile(const char* filename, 
						 globalFit &gFS,
						dataset &globalDataSet,
						 vector <double> &coefficients,
						 vector <unsigned int> &holdvector,
						 vector <double> &lowlim,
						 vector <double> &hilim){
	int err = 0;
	int ii, jj, highestParamSoFar, numdatasets;
	ifstream file_to_read;
	vector<string> columndata;
	
	vector<double> *localcoefficientsArray = NULL;
	vector<unsigned int> *localholdvectorArray = NULL;
	vector<double> *locallowlimArray = NULL;
	vector<double> *localhilimArray = NULL;

	string linein;
	
	//read the coefficient file, 1st two lines are headers.
	file_to_read.open(filename, ios::in);
	
	if(!file_to_read){
		cout << "ERROR global pilot file doesn't exist\n";
		return 1;
	}
	
	coefficients.clear();
	holdvector.clear();
	lowlim.clear();
	hilim.clear();

	//load each of the datafiles and put them in a global dataset
	if(!getline(file_to_read, linein)){
		cout << "error whilst reading datafiles from the global pilot file\n";
		err = 1;
		goto done;
	}
	
	Tokenize(linein, columndata, " ", sizeof(char));
	numdatasets = columndata.size();

	if(err = gFS.makedatasets(numdatasets))
		goto done;	

	globalDataSet.numDataDims = 2;
	for(ii = 0 ; ii < columndata.size() ; ii++){
		dataset tempdataset;
				
		if(err = tempdataset.readDataFile(columndata.at(ii).c_str()))
		   goto done;
		   
		globalDataSet.xx.insert( globalDataSet.xx.end(), tempdataset.xx.begin(), tempdataset.xx.end() );
		globalDataSet.yy.insert( globalDataSet.yy.end(), tempdataset.yy.begin(), tempdataset.yy.end() );
		globalDataSet.dy.insert( globalDataSet.dy.end(), tempdataset.dy.begin(), tempdataset.dy.end() );
		globalDataSet.dx.insert( globalDataSet.dx.end(), tempdataset.dx.begin(), tempdataset.dx.end() );
		globalDataSet.datapoints += tempdataset.datapoints;
		gFS.globalFitIndividualArray[ii].datafilename = columndata.at(ii);
		gFS.globalFitIndividualArray[ii].numDataDims = 2;
		gFS.globalFitIndividualArray[ii].datapoints = tempdataset.datapoints;
		//gFS.globalFitIndividualArray[ii].ffp = gaussian;
	}
	
	//now get the filenames for all the pilot files (coefficients/hold vectors, etc.
	if(!getline(file_to_read, linein)){
		cout << "error whilst reading individual pilot files from the global pilot file\n";
		err = 1;
		goto done;
	}
	columndata.clear();
	Tokenize(linein, columndata, " ", sizeof(char));
	
	//number of datasets has to be the same as the number of pilot files
	if(columndata.size() != gFS.numDataSets){
		cout << "error the number of datasets has to be equal to the number of pilot files in the global pilot file\n";
		err = 1;
		goto done;
	}
	
	for(ii = 0 ; ii < columndata.size() ; ii++)
		gFS.globalFitIndividualArray[ii].pilotfilename = columndata.at(ii);
	
	//lets read in the linkage matrix from the pilot file
	while(getline(file_to_read, linein)){
		columndata.clear();
		Tokenize(linein, columndata, " ", sizeof(char));

		for(ii = 0 ; ii < columndata.size() ; ii += 1){
			int val = strtol(columndata.at(ii).c_str(), NULL, 10);
			if(val >= 0)
				gFS.parameterLinkageTable[ii].push_back(val);
		}
	}
	
	//now read in the individual coefficient/pilot files.  They should have the same number of parameters as the linkage matrix.
	localcoefficientsArray = new (nothrow) vector<double>[gFS.numDataSets];
	locallowlimArray = new (nothrow) vector<double>[gFS.numDataSets];
	localhilimArray = new (nothrow) vector<double>[gFS.numDataSets];
	localholdvectorArray = new (nothrow) vector<unsigned int>[gFS.numDataSets];
	if(!localholdvectorArray || !locallowlimArray || !localhilimArray || !localholdvectorArray){
		err = NO_MEMORY;
		goto done;
	}

	for(ii = 0 ; ii < gFS.numDataSets ; ii += 1){
		if(err = readCoefficientFile(gFS.globalFitIndividualArray[ii].pilotfilename.c_str(),
									 localcoefficientsArray[ii],
									 localholdvectorArray[ii],
									 locallowlimArray[ii],
									 localhilimArray[ii])){
			cout << "Error whilst parsing one of the pilot files\n";
			goto done;
		}
		gFS.globalFitIndividualArray[ii].numcoefs = localcoefficientsArray[ii].size();
		if(gFS.globalFitIndividualArray[ii].numcoefs != gFS.parameterLinkageTable[ii].size()){
			err = 1;
			cout << "ERROR the number of parameters in the linkage table does not match the number of coeffcients in the pilot file\n";
			goto done;
		}
	}
	//now need to do sanity check on linkage arrays and coefficients
	if(err = sanityCheckOnLinkageValues(gFS, localcoefficientsArray))
		goto done;
	
	//create reduced parameter set for the global coefficients (i.e. remove linked coefficients).
	highestParamSoFar = -1;
	for(ii = 0 ; ii < gFS.numDataSets ; ii ++){
		for(jj = 0 ; jj < gFS.parameterLinkageTable[ii].size() ; jj++){
			if(gFS.parameterLinkageTable[ii].at(jj) == highestParamSoFar + 1){
				coefficients.push_back(localcoefficientsArray[ii].at(jj));
				holdvector.push_back(localholdvectorArray[ii].at(jj));
				lowlim.push_back(locallowlimArray[ii].at(jj));
				hilim.push_back(localhilimArray[ii].at(jj));
				highestParamSoFar += 1;
			}
		}
	}
	
done:
	if(file_to_read.is_open())
		file_to_read.close();
	if(localcoefficientsArray)
		delete[] localcoefficientsArray;
	if(localholdvectorArray)
		delete[] localholdvectorArray;
	if(locallowlimArray)
		delete[] locallowlimArray;
	if(localhilimArray)
		delete[] localhilimArray;
	
	
	return err;
};


//a fit function for making a global fit
int globalFitWrapper(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	globalFit *gFS = (globalFit*) userdata;	
	int val = 0, ii=0, jj=0;
	
	vector<double> individualCoefs;
	double *individualModel = model;
	const double **ourxdata = xdata;
	
	for(ii = 0 ; ii < gFS->numDataSets; ii++ ){
		
		globalFitIndividual *gFI = (gFS->globalFitIndividualArray + ii);
		individualCoefs.clear();
		
		//set up the coeffcients to be sent to the fitfunction
		for(jj = 0 ; jj < gFI->numcoefs ; jj ++ ){
			val = gFS->parameterLinkageTable[ii].at(jj);
			if(val < 0 || val > (numcoefs - 1)){
				err = LINKAGE_MATRIX_WRONG;
				goto done;
			}
			individualCoefs.push_back(coefs[val]);
		}
		if(individualCoefs.size() != gFI->numcoefs){
			cout << "ERROR, the expanded points did not match the number of fit parameters for one of the datasets\n";
			return 1;
		}
		
		if(err = (*(gFI->ffp))(NULL, (const double*) &individualCoefs[0], gFI->numcoefs, individualModel, (const double**)ourxdata, gFI->datapoints, gFI->numDataDims))
			goto done;
		
		//only increment the model pointer once you've called the function
		individualModel += gFI->datapoints;
		ourxdata += gFI->numDataDims;

	}
	
done:
	return err;
}

