
#include <vector>
#include "gencurvefit.h"

#define LINKAGE_MATRIX_WRONG 1

using namespace std;

int globalFitWrapper(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims);
double globalCostWrapper(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts);


class globalFitIndividual{
public:
	fitfunction ffp;
	costfunction costfun;
	int numcoefs;
	long datapoints;
	int numDataDims;
	string datafilename;
	string pilotfilename;
	string fitfunctionname;
	string costfunctionname;
	globalFitIndividual();
};

class globalFit{
public:
	int numDataSets;
	/*
	 linkage table should have the list of parameters and how they are linked.
	 it will have dimensions [numdatasets][number of params].
	 The entry in the linkage table will refer to which parameter should be chosen from the master coefficient table.
	 Please note that this array may be ragged if the number of coefficients for each of the individual datasets is not constant.
	 */
	vector<int> *parameterLinkageTable;
	globalFitIndividual *globalFitIndividualArray;
	
	globalFit();
	~globalFit();
	int makedatasets(int numDataSets);
};


int parseGlobalPilotFile(const char* filename,
						 globalFit &gFS, 
						 dataset &globalDataSet,
						 vector <double> &coefficients,
						 vector <unsigned int> &holdvector,
						 vector <double> &lowlim,
						 vector <double> &hilim
);

