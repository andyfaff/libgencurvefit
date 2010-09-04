/*
 *  genetic.c
 *  GeneticOptimisation
 *
 *
 */
#include "gencurvefit.h"
#include "mt19937p.h"
#include <stdlib.h>
#include <time.h>
#include "math.h"
#include "string.h"


/*
 A structure to hold statistics of an array
 */
struct waveStats {
	double V_avg;
	double V_stdev;
	long V_maxloc;
	long V_minloc;
};
typedef struct waveStats waveStats;
typedef struct waveStats* waveStatsPtr;

struct genoptStruct {
	/*total number of datapoints*/
	long datapoints;
	
	/*the (original) fit coeffcients*/
	double *coefs;
	const double *ydata;
	/*xdata is 2D, can fit multidimensional data*/
	const double **xdata;
	const double *edata;
	
	//the number of dimensions for the fit
	unsigned int numDataDims;
	
	/*place to put model once it's been calculated with the fitfunction*/
	double *model;
	
	//short array to hold the parameters that are being held
    const unsigned int *holdvector;
	
	/*a place to put a temporary copy of the coefficients, used to calculate the model*/
	double *temp_coefs;
	
	/*the fitfunction for the data*/
	fitfunction fitfun;
	
	/*the costfunction (chi2?)*/
	costfunction costfun;
	
	/*an optional user defined update function to keep appraised of the current fit status.*/
	updatefunction updatefun;
	unsigned int updatefrequency;
	
	/*how many parameters are being varied*/
	unsigned int numvarparams;
	
	/*total number of parameters*/
	unsigned int numcoefs;
	
	/*parameters for the optimisation*/
	unsigned int popsizeMultiplier;
	double k_m;
	double recomb;
	double tolerance;
	unsigned int iterations;
	
	/*totalsize of the population = popsizeMultiplier * numvarparams*/
	long totalpopsize;
	
	/*which parameters are varying*/
	unsigned int *varparams;
	
	/*an array which holds all the different guesses for the fit.*/
	/*it has dimensions popsize*numvarparams, numvarparams*/
	double **gen_populationvector;
	/*an individual genetic guess.*/
	double *gen_trial;
	/*which genetic strategy do you want?*/
	int strategy;
	/*Monte Carlo tempering parameter, set to NaN if not required*/
	double MCtemp;
	
	/*an array which holds all the chi2 values for all the different guesses in the population vector.*/
	double *chi2Array;
	/*the current chi2*/
	double chi2;
	/*number of fititerations done*/
	long numfititers;
	
	/*a copy of the limits that are being used.*/
	const double **limits;
	
	//have a pointer to data that users can pass around.
	void *userdata;
	
	struct mt19937p myMT19937;
};
typedef struct genoptStruct genoptStruct;
typedef struct genoptStruct* genoptStructPtr;

/*
 * \brief Create a two-dimensional array in a single allocation
 *
 * The effect is the same as an array of "element p[ii][jj];
 * The equivalent declaration is "element** p;"
 * The array is created as an array of pointer to element, followed by an array of arrays of elements.
 * \param ii first array bound
 * \param jj second array bound
 * \param sz size in bytes of an element of the 2d array
 * \return NULL on error or pointer to array
 *
 * assign return value to (element**)
 */

/* to use this in practice one would write 
 
 double **pp = NULL;
 pp = (double**)malloc2d(5, 11, sizeof(double));
 if(pp==NULL)
 return NOMEM;
 
 <use pp as required>
 free(pp);
 
 Note you can access elements by
 *(*(p+i)+j) is equivalent to p[i][j]
 In addition *(p+i) points to a whole row.
 */

 void* malloc2d(int ii, int jj, int sz)
{
	void** p;
	int sz_ptr_array;
	int sz_elt_array;
	int sz_allocation;
	long i;
	
	sz_ptr_array = ii * sizeof(void*);
	sz_elt_array = jj * sz;
	sz_allocation = sz_ptr_array + ii * sz_elt_array;
	
	p = (void**) malloc(sz_allocation);
	if (p == NULL)
		return p;
	memset(p, 0, sz_allocation);
	for (i = 0; i < ii; ++i)
	{
		*(p+i) = (void*) ((long)p + sz_ptr_array + i * sz_elt_array);
	}
	return p;
}


/*
 randomInteger returns an integer between 0 and upper EXclusive
 i.e. you will never get upper returned.
 */
static int randomInteger (struct mt19937p *myMT19937, int upper){
	int val;
	while (upper <= (val = genrand_int(myMT19937) / (0x7fffffff / upper)));
//	while (upper <= (val = random() / (RAND_MAX/upper)));
	return val;
}

/*
 randomDouble returns a double value between lower <= x <= upper OR [lower,upper]
 */
static double randomDouble(struct mt19937p *myMT19937, double lower, double upper){
	return lower + (upper - lower) * genrand(myMT19937);
	//double val = lower + random()/(((double)RAND_MAX + 1)/(upper - lower));
	//return val;
}


//returns gaussian noise, which is from a distribution with sd = 2.
double gnoise(struct mt19937p *myMT19937, double sd){
	double en0, en1;
	do{
		en0 = randomDouble(myMT19937, 0, 1);
	} while(en0==1);
	do{
		en1 = randomDouble(myMT19937, 0, 1);
	} while(en1 == 1);
	return sd * sqrt(-2 * log(en0))*cos(2 * PI * en1);
}

void SelectSamples(struct mt19937p *myMT19937, long popsize, long candidate, long *r1, long *r2, long *r3, long *r4, long *r5){
	if (r1){
		do
			*r1 = randomInteger(myMT19937, popsize);	
		while (*r1 == candidate);
	}
	
	if (r2)	{
		do
			*r2 = randomInteger(myMT19937, popsize);
		while ((*r2 == candidate) || (*r2 == *r1));
	}
	
	if (r3){
		do
			*r3 = randomInteger(myMT19937, popsize);
		while ((*r3 == candidate) || (*r3 == *r2) || (*r3 == *r1));
	}
	
	if (r4){
		do
			*r4 = randomInteger(myMT19937, popsize);
		while ((*r4 == candidate) || (*r4 == *r3) || (*r4 == *r2) || (*r4 == *r1));
	}
	
	if (r5){
		do
			*r5 = randomInteger(myMT19937, popsize);
		while ((*r5 == candidate) || (*r5 == *r4) || (*r5 == *r3)
			   || (*r5 == *r2) || (*r5 == *r1));
	}
	
	return;
}


void Best1Bin(genoptStruct *p, long candidate){
	long r1, r2;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate, &r1, &r2, NULL, NULL, NULL);
	n = randomInteger(&(p->myMT19937),  p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0; i < p->numvarparams; i++) {
		if ((randomDouble(&(p->myMT19937), 0, 1) < p->recomb) || (i == (p->numvarparams - 1)))
			*(p->gen_trial + n) = p->gen_populationvector[0][n]
			+ p->k_m * (p->gen_populationvector[r1][n] - p->gen_populationvector[r2][n]);
		
		n = (n + 1) % p->numvarparams;
	}
	return;
}

void Best1Exp(genoptStruct *p, long candidate){
	long r1, r2;
	long n, i;
	
	SelectSamples(&(p->myMT19937),  p->totalpopsize, candidate, &r1, &r2, NULL, NULL, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0 ; (randomDouble(&(p->myMT19937),  0, 1) < p->recomb) && (i < p->numvarparams); i++){
		*(p->gen_trial + n) = p->gen_populationvector[0][n]
		+ p->k_m * (p->gen_populationvector[r1][n] - p->gen_populationvector[r2][n]);
		
		n = (n + 1) % p->numvarparams;
	}
	return;
}

void Rand1Exp(genoptStruct *p, long candidate){
	long r1, r2, r3;
	long n, i;
	
	SelectSamples(&(p->myMT19937),  p->totalpopsize, candidate,&r1,&r2,&r3, NULL, NULL);
	n = randomInteger(&(p->myMT19937),  p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0; (randomDouble(&(p->myMT19937), 0, 1) < p->recomb) && (i < p->numvarparams); i++) {
		*(p->gen_trial + n) = p->gen_populationvector[r1][n]
		+ p->k_m * (p->gen_populationvector[r2][n] - p->gen_populationvector[r3][n]);
		
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void RandToBest1Exp(genoptStruct *p, long candidate){
	long r1, r2;
	long n,  i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2, NULL, NULL, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0; (randomDouble(&(p->myMT19937), 0, 1) < p->recomb) && (i < p->numvarparams); i++) {
		*(p->gen_trial + n) += p->k_m * (p->gen_populationvector[0][n] - *(p->gen_trial + n))
		+ p->k_m * (p->gen_populationvector[r1][n]
					   - p->gen_populationvector[r2][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void Best2Exp(genoptStruct *p, long candidate){
	long r1, r2, r3, r4;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2,&r3,&r4, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0; (randomDouble(&(p->myMT19937), 0, 1) < p->recomb) && (i < p->numvarparams); i++) {
		*(p->gen_trial + n) = p->gen_populationvector[0][n] +
		p->k_m * (p->gen_populationvector[r1][n]
					 + p->gen_populationvector[r2][n]
					 - p->gen_populationvector[r3][n]
					 - p->gen_populationvector[r4][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void Rand2Exp(genoptStruct *p, long candidate){
	long r1, r2, r3, r4, r5;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2,&r3,&r4,&r5);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0; (randomDouble(&(p->myMT19937), 0, 1) < p->recomb) && (i < p->numvarparams); i++) {
		*(p->gen_trial + n) = p->gen_populationvector[r1][n]
		+ p->k_m * (p->gen_populationvector[r2][n]
					   + p->gen_populationvector[r3][n]
					   - p->gen_populationvector[r4][n]
					   - p->gen_populationvector[r5][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void RandToBest1Bin(genoptStruct *p, long candidate){
	long r1, r2;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate, &r1, &r2, NULL, NULL, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	for (i=0; i < p->numvarparams; i++) 
	{
		if ((randomDouble(&(p->myMT19937), 0, 1) < p->recomb) || (i  == (p->numvarparams - 1)))
			*(p->gen_trial + n) += p->k_m * (p->gen_populationvector[0][n] - *(p->gen_trial + n))
			+ p->k_m * (p->gen_populationvector[r1][n]
						   - p->gen_populationvector[r2][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void Best2Bin(genoptStruct *p, long candidate){
	long r1, r2, r3, r4;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2,&r3,&r4, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	for (i=0; i < p->numvarparams; i++) 
	{
		if ((randomDouble(&(p->myMT19937), 0, 1) < p->recomb) || (i  == (p->numvarparams - 1)))
			*(p->gen_trial + n) = p->gen_populationvector[0][n]
			+ p->k_m * (p->gen_populationvector[r1][n]
						   + p->gen_populationvector[r2][n]
						   - p->gen_populationvector[r3][n]
						   - p->gen_populationvector[r4][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void Rand2Bin(genoptStruct *p, long candidate){
	long r1, r2, r3, r4, r5;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2,&r3,&r4,&r5);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	for (i=0; i < p->numvarparams; i++) 
	{
		if ((randomDouble(&(p->myMT19937), 0, 1) < p->recomb) || (i  == (p->numvarparams - 1)))
			*(p->gen_trial + n) = p->gen_populationvector[r1][n]
			+ p->k_m * (p->gen_populationvector[r2][n]
						   + p->gen_populationvector[r3][n]
						   - p->gen_populationvector[r4][n]
						   - p->gen_populationvector[r5][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

void Rand1Bin(genoptStruct *p, long candidate){
	long r1, r2, r3;
	long n, i;
	
	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2,&r3,NULL, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);
	
	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));
	
	for (i=0; i < p->numvarparams; i++) {
		if ((randomDouble(&(p->myMT19937), 0, 1) < p->recomb) || (i  == (p->numvarparams - 1)))
			*(p->gen_trial + n) =  p->gen_populationvector[r1][n]
			+ p->k_m * (p->gen_populationvector[r2][n]
						   - p->gen_populationvector[r3][n]);
		n = (n + 1) % p->numvarparams;
	}
	
	return;
}

/*
 createTrialVector makes a mutated vector.  It fills the trialVector from the current pvector and from bPrime,
 in modulo.
 bPrime is created from two random population vectors and the best fit vector.
 */
void createTrialVector(genoptStruct *p, long currentpvector){
	void (*theStrategy)(genoptStruct*, long);
	
	switch(p->strategy){
		case 0:
			theStrategy = Best1Bin;
			break;
		case 1:
			theStrategy = Best1Exp;
			break;
		case 2:
			theStrategy = Rand1Exp;
		case 3:
			theStrategy = RandToBest1Exp;
			break;
		case 4:
			theStrategy = Best2Exp;
			break;
		case 5:
			theStrategy = Rand2Exp;
			break;
		case 6:
			theStrategy = RandToBest1Bin;
			break;
		case 7:
			theStrategy = Best2Bin;
			break;
		case 8:
			theStrategy = Rand2Bin;
			break;
		case 9:
			theStrategy = Rand1Bin;
			break;
		default:
			theStrategy = Best1Bin;
			break;
	}
	
	theStrategy(p, currentpvector);
}

static waveStats getWaveStats(double *sort, long length,int moment){
	long ii = 0;
	double minval = *sort, maxval = *sort;
	long minpos = 0, maxpos = 0;
	double nx2 = 0, nx = 0;
	struct waveStats retval;
	
	switch(moment){
		case 0:
			for(ii = 0 ; ii < length ; ii += 1){
				if(*(sort + ii) > maxval){
					maxval = *(sort + ii);
					maxpos = ii;
				}
				if(*(sort + ii) < minval){
					minval = *(sort + ii);
					minpos = ii;
				}
			}
			retval.V_maxloc = maxpos;
			retval.V_minloc = minpos;
			break;
		case 1:
			for(ii = 0 ; ii < length ; ii += 1){
				nx += (*(sort + ii));
				nx2 += *(sort + ii) * (*(sort + ii));
				if(*(sort + ii) > maxval){
					maxval = *(sort + ii);
					maxpos = ii;
				}
				if(*(sort + ii) < minval){
					minval = *(sort + ii);
					minpos = ii;
				}
			}
			retval.V_maxloc = maxpos;
			retval.V_minloc = minpos;
			retval.V_avg = nx/(double)length;
			retval.V_stdev = sqrt((nx2 / (double)length) - (retval.V_avg * retval.V_avg));
			break;
	}
	return retval;
}

/*
 ensureConstraints takes the current trial vector and makes sure that all the individual 
 parameters lie inbetween the upper and lower limits.
 returns void
 */
static void
ensureConstraints(genoptStruct *p){
	unsigned int ii;	
	for(ii = 0 ; ii < p->numvarparams ; ii+=1)
		if(*(p->gen_trial + ii) < p->limits[*(p->varparams + ii)][0] || *(p->gen_trial + ii) > (p->limits[*(p->varparams + ii)][1]))
			*(p->gen_trial + ii) = randomDouble(&(p->myMT19937), p->limits[*(p->varparams + ii)][0], p->limits[*(p->varparams + ii)][1]);
}

/*
 insertVaryingParams inserts the current pvector into an array copy of the coefficients,
 then into a temporary wave
 returns 0 if no error
 returns errorcode otherwise
 */
static int
insertVaryingParams(genoptStruct *p, double *vector, int vectorsize){
	int err=0,ii;
	
	for(ii = 0 ; ii < p->numvarparams ; ii += 1)
		*(p->temp_coefs + *(p->varparams + ii)) =  *(vector + ii);

	return err;
}

/*
 setPopVector sets the populationvector with index replace, with a double array
 returns 0 if no error
 returns errorcode otherwise
 */
static int 
setPopVector(genoptStruct *p, double* vector, int vectorsize, long replace){
	//p->gen_populationvector[replace][] = vector[q]
	memcpy(*(p->gen_populationvector + replace), vector, vectorsize * sizeof(double));
	return 0;
}

/*
 swapChi2values swaps two values (i,j) in the p->chi2array 
 */
static void
swapChi2values(genoptStruct *p, long i, long j){
	double temp = *(p->chi2Array+i);
	*(p->chi2Array + i) = *(p->chi2Array + j);
	*(p->chi2Array + j) = temp;
}

/*
 swapPopVector swaps the i vector with index j in the populationvector
 returns 0 if no error
 returns errorcode otherwise
 */
static int
swapPopVector(genoptStruct *p, long popsize, long i, long j){
	double *tempparams = NULL;
	//do swap with pointers
	tempparams = *(p->gen_populationvector + j);
	*(p->gen_populationvector + j) = *(p->gen_populationvector + i);
	*(p->gen_populationvector + i) = tempparams;
	return 0;
}


/*initialises the genetic optimisation for the optimiseloop, all memory allocation has occurred by now, just need to fill out the arrays*/
int initialiseFit(genoptStruct *p){
	int err = 0;
	
	long ii, jj;
	double bot, top, chi2;
	waveStats wavStats;
	
	//initialise population vector guesses, from between the limits
	for(jj = 0 ; jj < p->numvarparams ; jj += 1){
		bot = p->limits[*(p->varparams + jj)][0];
		top = p->limits[*(p->varparams + jj)][1];
		for(ii = 0 ; ii < p->totalpopsize ; ii += 1)
			p->gen_populationvector[ii][jj] = randomDouble(&(p->myMT19937), bot, top);
	}
	
	
	//initialise Chi2array, will require a bit of calculation of the model function for each of the initial guesses.
	for(ii = 0 ; ii < p->totalpopsize ; ii += 1){
		if(err = insertVaryingParams(p, *(p->gen_populationvector+ ii), p->numvarparams))
			goto done;
		
		//calculate the model
		if(err = (*(p->fitfun))(p->userdata, p->temp_coefs, p->numcoefs, p->model, p->xdata, p->datapoints, p->numDataDims))
			goto done;
		
		/*calculate the costfunction*/
		chi2 = (*(p->costfun))(p->userdata, p->temp_coefs, p->numcoefs, p->ydata, p->model, p->edata, p->datapoints);
		
		*(p->chi2Array + ii)= chi2;
	}
	//find best chi2 and put that into number 0 pos.
	wavStats = getWaveStats(p->chi2Array, p->totalpopsize, 0);
	
	swapChi2values(p, 0, wavStats.V_minloc);
	if(err = swapPopVector(p, p->totalpopsize, 0, wavStats.V_minloc))
		goto done;
	
	//put the best fit from the intialisation into the coeffcients to return		   
	if(err = insertVaryingParams(p,  *(p->gen_populationvector), p->numvarparams))
		goto done;
	
	if(p->updatefun && (4 & p->updatefrequency))
		if(err = (*(p->updatefun))(p->userdata, p->temp_coefs, p->numcoefs, 0, *(p->chi2Array)))
			goto done;
	
	memcpy(p->coefs, p->temp_coefs, p->numcoefs * sizeof(double));
	
	
done:
	return err;
}

/*
 optimiseloop performs the optimisation.  It takes the initial population and mutates it until we find the best fit solution
 returns 0 if no errors
 returns errorcode otherwise.
 */
int optimiseloop(genoptStruct *p){
	long ii, kk, currentpvector;
	int err = 0;
	double chi2pvector,chi2trial;
	int acceptMoveGrudgingly;
	waveStats wavStats;
			
	/* the user sets how many times through the entire population*/
	for(kk = 1; kk <= p->iterations ; kk += 1){
		p->numfititers = kk;
		
		if(p->updatefun && (8 & p->updatefrequency)){
			if(err = insertVaryingParams(p, *(p->gen_populationvector), p->numvarparams))
				goto done;
			
			if(err = (*(p->updatefun))(p->userdata, p->temp_coefs, p->numcoefs, kk, *(p->chi2Array)))
				goto done;
		}
		
		/*iterate over all the individual members of the population*/
		for(ii = 0 ; ii < p->totalpopsize ; ii += 1){			
			currentpvector = ii;
			/*now set up the trial vector using a wave from the populationvector and bprime
			//first set the pvector 
			// create a mutated trial vector from the best fit and two random population members*/
			createTrialVector(p, currentpvector);
			/*/ make sure the trial vector parameters lie between the user defined limits*/
			ensureConstraints(p);
			
			chi2pvector = *(p->chi2Array + ii);
			/*
			 find out the chi2 value of the trial vector		
			 */
			
			if(err = insertVaryingParams(p, p->gen_trial, p->numvarparams))
				goto done;
			
			if(err = (*(p->fitfun))(p->userdata, p->temp_coefs, p->numcoefs, p->model, p->xdata, p->datapoints, p->numDataDims))
				goto done;

			/*calculate the costfunction*/
			chi2trial = (*(p->costfun))(p->userdata, p->temp_coefs, p->numcoefs, p->ydata, p->model, p->edata, p->datapoints);

			acceptMoveGrudgingly = 0;
			if(isfinite(p->MCtemp) && p->MCtemp > 0 && (exp(-chi2trial / chi2pvector / p->MCtemp) < randomDouble(&(p->myMT19937), 0, 1)) ){
				acceptMoveGrudgingly = 1;				
				if(p->updatefun && (2 & p->updatefrequency))
					if(err = (*(p->updatefun))(p->userdata, p->temp_coefs, p->numcoefs, kk, chi2trial))
						goto done;
			}
			
			/*
			 if the chi2 of the trial vector is less than the current populationvector then replace it
			 */
			if(chi2trial < chi2pvector  || (acceptMoveGrudgingly && ii)){
				if(err = setPopVector(p, p->gen_trial, p->numvarparams, currentpvector))
					goto done;
				
				*(p->chi2Array + ii) = chi2trial;
				/*
				 if chi2 of the trial vector is less than that of the best fit, then replace the best fit vector
				 */
				if(chi2trial < *(p->chi2Array)){		/*if this trial vector is better than the current best then replace it*/
					if(err = setPopVector(p, p->gen_trial, p->numvarparams, 0))
						goto done;
					
					/*
					 a user defined update function that can be used to halt the fit early, and keep 
					  appraised of fit progress.
					*/
					if(p->updatefun && (1 & p->updatefrequency))
						if(err = (*(p->updatefun))(p->userdata, p->temp_coefs, p->numcoefs, kk, chi2trial))
							goto done;
	
					/*
					 if the fractional decrease in chi2 is less than the tolerance then abort the fit
					 */
					wavStats = getWaveStats(p->chi2Array, p->totalpopsize, 1);
					
					/*
					 update the best chi2 if you've just found a better fit (but not yet reached termination
					 */
					*(p->chi2Array) = chi2trial;
					
					if( wavStats.V_stdev/wavStats.V_avg < p->tolerance)	/*if the fractional decrease is less and 0.5% stop.*/
						goto done;
				}
			}
		}
	}
	
done:
	return err;
}

int genetic_optimisation(fitfunction fitfun,
							 costfunction costfun,
							 unsigned int numcoefs,
							 double* coefs,
							 const unsigned int *holdvector,
							 const double** limits,
							 long datapoints,
							 const double* ydata,
							 const double** xdata,
							 const double *edata,
							 unsigned int numDataDims,
							 double *chi2,
							 const gencurvefitOptions* gco,	 
							 void* userdata
							){
	int err = 0, popsizeMultiplier = 20;
	
	unsigned int ii,jj;
	
	/*the overall structure to contain the entire fit*/
	genoptStruct gos;
	memset(&gos, 0, sizeof(gos));
		
	//initialise the random number generators
//	srandom(clock());
	if(gco == NULL || gco->seed < 1)
		sgenrand(clock(), &(gos.myMT19937));
	else
		sgenrand(gco->seed, &gos.myMT19937);
	
	//setup the data.
	gos.xdata = xdata;
	gos.ydata = ydata;
	gos.edata = edata;
	gos.datapoints = datapoints;
	gos.numDataDims = numDataDims;
	
	//setup the function pointers
	gos.fitfun = fitfun;
	if(!costfun)
		gos.costfun = &chisquared;
	else
		gos.costfun = costfun;

	//check that the total number of parameters matches the number of parameters in the holdvector
	gos.numcoefs = numcoefs;
	gos.coefs = coefs;
	gos.holdvector = holdvector;
		
	//add in userdata
	gos.userdata = userdata;
		
	//work out which parameters are being held
	for(ii = 0 ; ii < numcoefs ; ii += 1)
		if(holdvector[ii] == 0)
		    gos.numvarparams += 1;
	
	if(gos.numvarparams < 1){
		err = NO_VARYING_PARAMS;
		goto done;
	}
	
	//put the parameter numbers that are being held into an array
	gos.varparams = (unsigned int*) malloc (gos.numvarparams * sizeof(int));
	if(gos.varparams == NULL){
		err = NO_MEMORY;
		goto done;
	}
	jj=0;
	for(ii = 0 ; ii < numcoefs ; ii += 1){
		if(holdvector[ii] == 0){
		    gos.varparams[jj] = ii;
			jj+=1;
		}
	}
	
	//these are the upper and lower bounds for the limits, check them
	for(ii = 0 ; ii < numcoefs ; ii += 1){
		if(holdvector[ii] == 0 && limits[ii][0] > limits[ii][1]){
			err = INCORRECT_LIMITS;
			goto done;
		}
	}
	gos.limits = limits;
	
	//optional parameters
	if(!gco){
		gos.updatefun = NULL;
		gos.MCtemp = NAN;
		popsizeMultiplier = 20;
		gos.k_m = 0.7;
		gos.recomb = 0.5;
		gos.popsizeMultiplier = 20;
		gos.iterations = 100;
		gos.tolerance = 0.001;
		gos.updatefrequency = 1;
	} else {
		gos.updatefun = gco->updatefun;
		if(gco->temp <= 0)
			gos.MCtemp = NAN;
		else
			gos.MCtemp = gco->temp;
		popsizeMultiplier = gco->popsizeMultiplier;
		gos.popsizeMultiplier = gco->popsizeMultiplier;
		gos.k_m = gco->k_m;
		gos.recomb = gco->recomb;
		gos.iterations = gco->iterations;
		gos.tolerance = gco->tolerance;
		gos.updatefrequency = gco->updatefrequency;
	}
	gos.totalpopsize = gos.numvarparams * popsizeMultiplier;

	
	//now we can allocate memory for the rest of the structure members in gos.
	//initialise population vector
	gos.gen_populationvector = (double**)malloc2d(gos.totalpopsize, gos.numvarparams, sizeof(double));
	if(gos.gen_populationvector == NULL){
		err = NO_MEMORY;
		goto done;
	}
	//initialise Chi2array
	gos.chi2Array = (double*) malloc (gos.totalpopsize * sizeof(double));
	if(gos.chi2Array == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memset(gos.chi2Array, 0, gos.totalpopsize);
	
	//initialise the trial vector
	gos.gen_trial = (double*)malloc(gos.numvarparams * sizeof(double));
	if(gos.gen_trial == NULL){
		err = NO_MEMORY;
		goto done;
	}

	//initialise space for a full array copy of the coefficients
	gos.temp_coefs = (double*)malloc(numcoefs * sizeof(double));
	if(gos.temp_coefs == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memcpy(gos.temp_coefs, coefs, numcoefs*sizeof(double));
	
	//initialise space for a full array copy of the coefficients
	gos.model = (double*)malloc(datapoints* sizeof(double));
	if(gos.model == NULL){
		err = NO_MEMORY;
		goto done;
	}	
	
	//at this stage we can initialise the fit, before looping.
	if(err = initialiseFit(&gos))
		goto done;
	
	//now iterate through, fitting the data.
	err = optimiseloop(&gos);

	//at this point we have the best fit.  Now return the results
	*chi2 = *(gos.chi2Array);
	
	/*
	 put the best fit into the coeffcient array that will be returned to the user
	 */
	insertVaryingParams(&gos,  *(gos.gen_populationvector), gos.numvarparams);
	memcpy(gos.coefs, gos.temp_coefs, gos.numcoefs * sizeof(double));
	
	
done:

	if(gos.chi2Array)
		free(gos.chi2Array);
	if(gos.gen_populationvector)
		free(gos.gen_populationvector);
	if(gos.gen_trial)
		free(gos.gen_trial);
	if(gos.model)
		free(gos.model);
	if(gos.temp_coefs)
		free(gos.temp_coefs);
	if(gos.varparams)
		free(gos.varparams);
	
	return err;
}

double chisquared(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts){
	long ii;
	double chi2 = 0;
	double val=0;
	for (ii=0; ii<numpnts; ii+=1){
		double temp1, temp2, temp3;
		temp1 = data[ii];
		temp2 = model[ii];
		temp3 = errors[ii];
		val = pow((fabs((data[ii] - model[ii])/errors[ii])),2);
		if(isfinite(val))
			chi2 += val;
	}
	
	return chi2;
}

double robust(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts){
 long ii;
 double chi = 0;
 double val=0;
 for (ii=0; ii<numpnts; ii+=1){
	 val = fabs((data[ii] - model[ii])/errors[ii]);
	 if(isfinite(val))
		 chi += val;
 }
 return chi;
 };
 
