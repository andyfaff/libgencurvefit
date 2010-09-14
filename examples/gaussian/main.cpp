
#include <iostream>
#include <fstream>
#include <vector>
#include <sstream>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "gencurvefit.h"
#include "myfitfunctions.h"
#include "dataset.h"

#ifdef USE_MPI
#include <mpi.h>
#endif
#include <omp.h>

using namespace std;

int NUM_CPUS = 1;
#define FIT_FUNCTION gaussian
#define COST_FUNCTION chisquared
#define GO_TOL 0.0005
#define GO_KM 0.7
#define GO_RECOMB 0.5
#define GO_POPSIZEMULTIPLIER 20
#define GO_ITERS 500

typedef struct{
	const double *coefP;
	int numcoefs;
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
	
	double **coefResults;	//this is going to be successive fits, line after line.  i.e. have dimensions [numiters][coefnum];
	double *chi2Results;
	int iterationsToPerform;
	int iterationOffset;
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

void fitWorker(void* arg) { 
	int err = 0;
	int jj = 0;
	int ii = 0;
	
	fitWorkerParm *p = (fitWorkerParm *) arg;
	gencurvefitOptions gco;
	
	double *coefsTemp = NULL;
	double *yytemp = NULL;
	double *eetemp = NULL;
	
	string outputString;	
	double chi2;
	
	coefsTemp = (double *)malloc(sizeof(double)*p->numcoefs);
	if(!coefsTemp){
		err= NO_MEMORY;
		goto done;
	}
	memcpy(coefsTemp, p->coefP, sizeof(double) * p->numcoefs);
	
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
	
	if(p->useErrors)
		memcpy(eetemp, p->eP, sizeof(double) * p->datapoints);
	else{
		for(ii=0 ; ii< p->datapoints ; ii+=1)
			eetemp[ii] = 1.0;
	}
	
	gco.iterations = GO_ITERS;
	gco.popsizeMultiplier = GO_POPSIZEMULTIPLIER;
	gco.k_m = GO_KM;
	gco.recomb = GO_RECOMB;
	gco.tolerance = GO_TOL;
	gco.strategy = 0;
	gco.temp = 1;
	gco.updatefun = NULL;
	gco.updatefrequency = 0;
	gco.seed = -1;
	gco.useinitialguesses = 0;
	
	//at this point we have 3 columns of data and the coefficients
	//we can start doing the fit.
	//do a load of montecarlo iterations
	for(jj= 0 ; jj < p->iterationsToPerform ; jj+=1){
		//add on the gaussian noise, we're going to be fitting on a log10 scale as well
		//		for(ii=0 ; ii< p->datapoints ; ii+=1)
		//			yytemp[ii] = p->yP[ii] + gnoise(p->eP[ii]);
		
		//do the genetic optimisation
		err = levenberg_marquardt(p->fitfun,
								  p->costfun,
								  p->numcoefs,
								  coefsTemp,
								  p->holdvector,
								  p->datapoints,
								  yytemp,
								  p->xP,
								  eetemp,
								  1,
								  &chi2,
								  NULL,
								  p->userdata); 
/*		
		err = genetic_optimisation(p->fitfun,
								   p->costfun,
								   p->numcoefs,
								   coefsTemp,
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
*/	
		//output the results
		if(err){
			cout << err;
			goto done;
		}
		outputString.clear();
		
		outputString.append(to_a_string(&chi2, 1));
		outputString.append(" ");
		outputString.append(to_a_string(coefsTemp, p->numcoefs));
		outputString.append("\r\n");
		cout << outputString;
		cout.flush();
	}
	
done:
	if(yytemp)
		free(yytemp);
	if(eetemp)
		free(eetemp);
	if(coefsTemp)
		free(coefsTemp);
	
};


int main (int argc, char *argv[]) {
	int err = 0;
	long ii;
	
	vector<double> xx;
	vector<double> yy;
	
	vector<double> fityy;	
	vector<double> ee;
	vector<double> coefs;
	vector<double> lowlim;
	vector<double> hilim;
	double **limits = NULL;
	string linein;
	vector<string> columndata;
	vector<unsigned int> bs;
	ifstream file_to_read;
	bool useErrors;
	
	int MCiters, myMCiters;
	
	double **xdata = NULL;
	double **fittedCoefs = NULL;
	double *fittedChi2 = NULL;
	
	//create threads to do each fit.
	fitWorkerParm *MC_arg = NULL;
	
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
	
	if(argc != 5){
		cout << "Useage:\n ./motoMC data.txt pilot useErrors mciters\n";
		err = WRONG_NUMBER_OF_PARAMS;
		goto done;
	}
	
	useErrors = strtol(argv[3], NULL, 10);
	MCiters = (int) strtol(argv[4], NULL, 10);
	myMCiters = MCiters / numprocs;
	
	if(myid == numprocs - 1)
		myMCiters += MCiters - numprocs * myMCiters;
	
	ii=0;
	//read the data
	file_to_read.open(argv[1], ios::in);
	if(file_to_read.is_open()){		
		while(getline(file_to_read, linein, '\n')){
			Tokenize(linein, columndata, "\t", sizeof(char));			
		
			xx.push_back(strtod(columndata[0].c_str(), NULL));
			yy.push_back(strtod(columndata[1].c_str(), NULL));
			ee.push_back(fabs(strtod(columndata[2].c_str(), NULL)));
			columndata.clear();
			   ii+=1;
		}
		file_to_read.close();
	} else {
		goto done;
	}
	
	//read the coefficient file, 1st two lines are headers.
	file_to_read.open(argv[2], ios::in);
	columndata.clear();
	if(file_to_read.is_open()){
		getline(file_to_read, linein, '\n');		//chi2value header
		getline(file_to_read, linein, '\n');		//header describing columns
		
		while(getline(file_to_read, linein, '\n')){
			Tokenize(linein, columndata, " ", sizeof(char));
			
			coefs.push_back(strtod(columndata[0].c_str(), NULL));
			bs.push_back(strtol(columndata[1].c_str(), NULL, 10));			
			lowlim.push_back(strtod(columndata[2].c_str(), NULL));
			hilim.push_back(strtod(columndata[3].c_str(), NULL));
			columndata.clear();
		}
		file_to_read.close();
	} else {
		goto done;
	}
	
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
	
	//set up the xpoints
	xdata = (double**)malloc2d(1, yy.size(), sizeof(double));
	if(!xdata){
		err = NO_MEMORY;
		goto done;
	}
	for(ii=0 ; ii < yy.size(); ii+=1)
		xdata[0][ii] = xx[ii];
	
	//allocate memory for the fit parameter results.
	fittedCoefs = (double**)malloc2d(myMCiters, coefs.size(), sizeof(double));
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
	
	//allocate memory for the thread arguments.
	MC_arg = (fitWorkerParm *) malloc(sizeof(fitWorkerParm) * myMCiters);
	if(!MC_arg){
		err = NO_MEMORY;
		goto done;
	}
	memset(MC_arg, 0, sizeof(fitWorkerParm));
	
#pragma omp parallel for shared(MC_arg) private(ii) 
	for (ii = 0 ; ii < myMCiters ; ii++){					
		MC_arg[ii].iterationsToPerform = 1;
		MC_arg[ii].useErrors = useErrors;
		MC_arg[ii].iterationOffset = ii;
		MC_arg[ii].fitfun = &FIT_FUNCTION;
		MC_arg[ii].costfun = &COST_FUNCTION;
		MC_arg[ii].holdvector = &bs[0];
		MC_arg[ii].chi2Results = fittedChi2 + ii;
		MC_arg[ii].coefResults = *(fittedCoefs + ii); //&fittedcoefs[iterations_consumed][0]
		MC_arg[ii].coefP = &coefs[0];
		MC_arg[ii].numcoefs = coefs.size();
		MC_arg[ii].dimensions = 1;
		MC_arg[ii].datapoints = yy.size();
		MC_arg[ii].yP = &yy[0];
		MC_arg[ii].xP = (const double **) xdata;
		MC_arg[ii].eP = &ee[0];
		MC_arg[ii].limits = (const double **) limits;
		MC_arg[ii].userdata = NULL;
		
		fitWorker((void*)(MC_arg + ii));
	}
	
#ifdef USE_MPI	
	MPI_Finalize(); /* MPI Programs end with MPI Finalize; this is a weak synchronization point */
#endif
done:
	
	if(file_to_read.is_open())
		file_to_read.close();
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



