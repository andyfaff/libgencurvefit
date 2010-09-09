
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <ctime>

//#include "gencurvefit.h"
#include "dataset.h"
#include "globalfitfunction.h"
#include "myfitfunctions.h"


#ifdef USE_MPI
#include <mpi.h>
#endif
#include <omp.h>

using namespace std;

int NUM_CPUS = 1;
#define FIT_FUNCTION globalFitWrapper
#define FIT_FUNCTION_INDIVIDUAL Abeles
#define COST_FUNCTION chisquared
#define GO_TOL 0.0005
#define GO_KM 0.7
#define GO_RECOMB 0.5
#define GO_POPSIZEMULTIPLIER 20
#define GO_ITERS 500
#define GO_STRATEGY 0


typedef struct{
	const double *coefP;
	unsigned int numcoefs;
	const double **limits;
	const double *yP;
	const double **xP;
	const double *eP;
	bool useErrors;
	const unsigned int *holdvector;
	
	long datapoints;
	int dimensions;
	
	fitfunction fitfun;
	costfunction costfun;
	long popsizeMultiplier;
	float k_m;
	float k_recomb;
	float tol;
	int gen_iters;
	
	double *coefResults;	//this is going to be where the coefficients for the fit are put
	double *chi2Results;
	void* userdata;
}  fitWorkerParm;

std::string to_a_string(double *nums, long numthings)
{
	std::ostringstream oss;
	long ii;
	for(ii=0 ; ii< numthings ; ii+=1)
		oss << nums[ii] << " ";
	
	return oss.str();
}	


void fitWorker(fitWorkerParm* p) { 
	int err = 0;
	
	int ii;
	gencurvefitOptions gco;
	double *yytemp = NULL;
	double *eetemp = NULL;
	
	string outputString;	
	double chi2;
	
	//copy the coefficients from the supplied coefficient vector into the results.
	memcpy(p->coefResults, p->coefP, sizeof(double) * p->numcoefs);
	
	yytemp = (double *)malloc(sizeof(double)*p->datapoints);
	if(!yytemp){
		err= NO_MEMORY;
		goto done;
	}
	memcpy(yytemp, p->yP, sizeof(double)*p->datapoints);
	
	eetemp = (double *)malloc(sizeof(double)*p->datapoints);
	if(!eetemp){
		err= NO_MEMORY;
		goto done;
	}
	memcpy(eetemp, p->eP, sizeof(double) * p->datapoints);
	
	if(p->useErrors)
		for(ii=0 ; ii< p->datapoints ; ii+=1)
			eetemp[ii] = log10((p->yP[ii] + p->eP[ii]) / p->yP[ii]);
	else{
		for(ii=0 ; ii< p->datapoints ; ii+=1)
			eetemp[ii] = 1.0;
	}
//
//	//add on the gaussian noise, we're going to be fitting on a log10 scale as well
//	//		for(ii=0 ; ii< p->datapoints ; ii+=1)
//	//			yytemp[ii] = log10(p->yP[ii] + gnoise(p->eP[ii]));
	for(ii=0 ; ii< p->datapoints ; ii+=1)
		yytemp[ii] = log10(p->yP[ii]);
	
	
	memset(&gco, 0, sizeof(gencurvefitOptions));
	gco.tolerance = p->tol;
	gco.k_m = p->k_m;
	gco.recomb = p->k_recomb;
	gco.popsizeMultiplier = p->popsizeMultiplier;
	gco.temp = -1;
	gco.iterations = p->gen_iters;
	gco.strategy = GO_STRATEGY;
			
	//at this point we have 3 columns of data and the coefficients
	//we can start doing the fit.
	//do a load of montecarlo iterations

	err = genetic_optimisation(p->fitfun,
							   p->costfun,
							   p->numcoefs,
							   p->coefResults,
							   p->holdvector,
							   p->limits,
							   p->datapoints,
							   yytemp,
							   p->xP,
							   eetemp,
							   1, 
							   &chi2,
							   &gco,
							   p->userdata);
	
	//output the results
	if(err){
		cout << err;
		goto done;
	}
	outputString.clear();
	
	outputString.append(to_a_string(&chi2, 1));
	outputString.append(" ");
//	outputString.append(to_a_string((double*)p->userdata, 1));
//	outputString.append(" ");
	outputString.append(to_a_string(p->coefResults, p->numcoefs));
	outputString.append("\r\n");
	cout << outputString;
	cout.flush();
	
done:
	if(yytemp)
		free(yytemp);
	if(eetemp)
		free(eetemp);
};


int main (int argc, char *argv[]) {
	int err = 0;
	long ii, highestNumberOfPoints = 0, pointOffset = 0;
	
	globalFit gFS;
	dataset theDataSet;
	
	vector<double> fityy;	

	vector<double> coefs;
	vector<double> lowlim;
	vector<double> hilim;
	double **limits = NULL;
	vector<unsigned int> bs;
	bool useErrors;
	string crap;
	int MCiters, myMCiters;
	double **xdata = NULL;
	double **fittedCoefs = NULL;
	double *fittedChi2 = NULL;
	fitWorkerParm *MC_arg = NULL;
	
	time_t time1, time2;
	
	int numprocs = 1;
	int myid = 0;	
#ifdef USE_MPI
	MPI_Init(&argc,&argv); /* all MPI programs start with MPI_Init; all 'N' processes exist thereafter */
	MPI_Comm_size(MPI_COMM_WORLD,&numprocs); /* find out how big the SPMD world is */
	MPI_Comm_rank(MPI_COMM_WORLD, &myid); /* and this processes' rank is */
#else
	numprocs = 1;
	myid = 0;
#endif
	
	if(argc != 4){
		cout << "Useage:\n ./motoMC globalpilot useerrors iterations\n";
		err = WRONG_NUMBER_OF_PARAMS;
		goto done;
	}
	
	useErrors = strtol(argv[2], NULL, 10);
	MCiters = (int) strtol(argv[3], NULL, 10);
	myMCiters = MCiters / numprocs;

	if(myid == numprocs - 1)
		myMCiters += MCiters - numprocs * myMCiters;
	
	time(&time1);
	
	//set up global fit structure
	if(err = parseGlobalPilotFile(argv[1], gFS, theDataSet, coefs, bs, lowlim, hilim))
		goto done;
	
	//setup the limits array
	limits = (double**)malloc2d(2, coefs.size(), sizeof(double));
	if(!limits){
		err = NO_MEMORY;
		goto done;
	}
	for(ii=0 ; ii<coefs.size(); ii+=1){
		limits[0][ii] = lowlim[ii];
		limits[1][ii] = hilim[ii];		
	}
	
	//say what fit function you want to use
	for(ii = 0 ; ii < gFS.numDataSets ; ii++)
		gFS.globalFitIndividualArray[ii].ffp = FIT_FUNCTION_INDIVIDUAL;
	
	//we have to put the xdata for the global fit wave in an array that the globalfitwrapper can understant.
	//currently they are in a vector.  We need them in a 2D array, where the rows are each dataset and the columns the datapoints
	//I can't be bothered making a ragged array, so lets just make it square, with the largest number of datapoints determining the
	//column size.
	for(ii = 0 ; ii < gFS.numDataSets ; ii++)
		if(gFS.globalFitIndividualArray[ii].datapoints > highestNumberOfPoints)
			highestNumberOfPoints = gFS.globalFitIndividualArray[ii].datapoints;

	xdata = (double**)malloc2d(gFS.numDataSets, highestNumberOfPoints, sizeof(double));
	if(!xdata){
		err = NO_MEMORY;
		goto done;
	}
	for(ii = 0 ; ii < gFS.numDataSets ; ii++){
		long datasetpoints = gFS.globalFitIndividualArray[ii].datapoints;
		memcpy(*(xdata + ii), &(theDataSet.xx[pointOffset]), datasetpoints * sizeof(double));
		pointOffset += datasetpoints;
	}
	
	//allocate memory for the fit parameter results.
	fittedCoefs = (double**) malloc2d(myMCiters, coefs.size(), sizeof(double));
	if(!fittedCoefs){
		err = NO_MEMORY;
		goto done;
	}
	
	//and allocate memory for the fitted chi2 value.
	fittedChi2 = (double*)malloc(myMCiters * sizeof(double));
	if(!fittedChi2){
		err = NO_MEMORY;
		goto done;
	}
	
	//allocate memory for the (openMP) thread arguments.
	MC_arg = (fitWorkerParm *) malloc(sizeof(fitWorkerParm) * myMCiters);
	if(!MC_arg){
		err = NO_MEMORY;
		goto done;
	}
	memset(MC_arg, 0, sizeof(fitWorkerParm));

#pragma omp parallel for shared(MC_arg) private(ii) 
	for (ii = 0 ; ii < myMCiters ; ii++){					
		MC_arg[ii].useErrors = useErrors;
		MC_arg[ii].fitfun = &FIT_FUNCTION;
		MC_arg[ii].costfun = &COST_FUNCTION;
		MC_arg[ii].holdvector = &bs[0];
		MC_arg[ii].chi2Results = fittedChi2 + ii;
		MC_arg[ii].coefResults = *(fittedCoefs + ii);
		MC_arg[ii].coefP = &coefs[0];
		MC_arg[ii].numcoefs = coefs.size();
		MC_arg[ii].dimensions = 1;
		MC_arg[ii].datapoints = theDataSet.datapoints;
		MC_arg[ii].yP = &(theDataSet.yy[0]);
		MC_arg[ii].xP = (const double **) xdata;

		MC_arg[ii].eP = &(theDataSet.dy[0]);
		MC_arg[ii].limits = (const double **) limits;
		
		MC_arg[ii].tol = GO_TOL;
		MC_arg[ii].k_m = GO_KM;
		MC_arg[ii].k_recomb = GO_RECOMB;
		MC_arg[ii].popsizeMultiplier = GO_POPSIZEMULTIPLIER;
		MC_arg[ii].gen_iters = GO_ITERS;
		MC_arg[ii].userdata = &gFS;
		
		fitWorker(MC_arg + ii);
	}
	
#ifdef USE_MPI	
	MPI_Finalize(); /* MPI Programs end with MPI Finalize; this is a weak synchronization point */
#endif
	time(&time2);
//	cout << difftime(time2, time1) << "\n";
done:
	
	if(limits)
		free(limits);
	if(MC_arg)
		free(MC_arg);
	if(fittedCoefs)
		free(fittedCoefs);
	if(fittedChi2)
		free(fittedChi2);
	
    return err;
}



