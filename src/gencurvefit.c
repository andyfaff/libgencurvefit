/*
 *  genetic.c
 *  GeneticOptimisation
 *
 *
 */
#include "gencurvefit.h"
#include "mt19937p.h"
#include "stdlib.h"
#include <time.h>
#include "math.h"
#include "string.h"

#ifdef _WINDOWS
#ifndef __GNUC__
#include "float.h"

	int isfinite(double x){
		return _finite(x);
	}
#endif
#endif

#define MIN(x, y) (((x) < (y)) ? (x) : (y))

/**
 A structure to hold statistics of an array
 */
struct waveStats {
	double V_avg;
	double V_stdev;
	long V_maxloc;
	long V_minloc;
};
typedef struct waveStats waveStats;

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

	/*
	a place to put a temporary copy of the coefficients, used to calculate the
	model
	*/
	double *temp_coefs;

	/*the fitfunction for the data*/
	fitfunction fitfun;

	/*the costfunction (chi2?)*/
	costfunction costfun;

	/*
	 an optional user defined update function to keep appraised of the current
	 fit status.
	 */
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
	double dither[2];

	/*totalsize of the population = popsizeMultiplier * numvarparams*/
	unsigned int totalpopsize;

	/*did you want to use the intial guesses to initialise the fit*/
	int useinitialguesses;
    
    /*did you supply an initial population  to use*/
    const double *initial_population;
    unsigned int initial_population_rows;

	/*which parameters are varying*/
	unsigned int *varparams;

	/*
	 an array which holds all the different guesses for the fit.
	it has dimensions popsize*numvarparams, numvarparams
	it is scaled between 0 and 1
	 */
	double **gen_populationvector;

	/*an individual genetic guess.*/
	double *gen_trial;

    /*scale factors*/
    double **scale_factors;
    
	/*which genetic strategy do you want?*/
	int strategy;

	/*
	an array which holds all the chi2 values for all the different guesses in
	the population vector.
	*/
	double *chi2Array;

	/*the current chi2*/
	double chi2;

	/*number of fititerations done*/
	long numfititers;

	/*a copy of the limits that are being used.*/
	const double **limits;

	//have a pointer to data that users can pass around.
	void *userdata;

	mt19937p myMT19937;
};
typedef struct genoptStruct genoptStruct;

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
 void* malloc2d(int ii, int jj, int sz){
    void** p;
    size_t sz_ptr_array;
    size_t sz_elt_array;
    size_t sz_allocation;
    long i = 0;
    char *c = NULL;
    
    sz_ptr_array = ii * sizeof(void*);
    sz_elt_array = jj * sz;
    sz_allocation = sz_ptr_array + ii * sz_elt_array;
    
    p = (void**) malloc(sz_allocation);
    if (p == NULL)
        return p;
    memset(p, 0, sz_allocation);
    
    c = ((char*) p) + sz_ptr_array;
    for (i = 0; i < ii; ++i)
    {
        //*(p+i) = (void*) ((long)p + sz_ptr_array + i * sz_elt_array);
        p[i] = (void*) (c + i * sz_elt_array);
    }
    return p;
}


/*
 randomInteger returns an integer between 0 and upper EXclusive
 i.e. you will never get upper returned.
 */
unsigned long randomInteger (mt19937p *myMT19937, unsigned long upper){
	unsigned long val;
	while (upper <= (val = genrand_int(myMT19937) / (0xffffffff / upper)));
//	while (upper <= (val = random() / (RAND_MAX/upper)));
	return val;
}

/*
 randomDouble returns a double value between lower <= x <= upper OR [lower,upper]
 */
double randomDouble(mt19937p *myMT19937, double lower, double upper){
	return lower + (upper - lower) * genrand(myMT19937);
	//double val = lower + random()/(((double)RAND_MAX + 1)/(upper - lower));
	//return val;
}

//returns gaussian noise, which is from a distribution with sd = 2.
double gnoise(mt19937p *myMT19937, double sd){
	double en0, en1;
	do{
		en0 = randomDouble(myMT19937, 0, 1);
	} while(en0==1);
	do{
		en1 = randomDouble(myMT19937, 0, 1);
	} while(en1 == 1);
	return sd * sqrt(-2 * log(en0))*cos(2 * PI * en1);
}

static void SelectSamples(mt19937p *myMT19937,
                          long popsize,
                          long candidate,
                          long *r1,
                          long *r2,
                          long *r3,
                          long *r4,
                          long *r5){
	if (r1){
		do
			*r1 = randomInteger(myMT19937, popsize);
		while (*r1 == candidate);
	}

	if (r2)	{
		do
			*r2 = randomInteger(myMT19937, popsize);
		while ((*r2 == candidate)
		        || (*r2 == *r1));
	}

	if (r3){
		do
			*r3 = randomInteger(myMT19937, popsize);
		while ((*r3 == candidate)
		       || (*r3 == *r2)
		       || (*r3 == *r1));
	}

	if (r4){
		do
			*r4 = randomInteger(myMT19937, popsize);
		while ((*r4 == candidate)
		       || (*r4 == *r3)
		       || (*r4 == *r2)
		       || (*r4 == *r1));
	}

	if (r5){
		do
			*r5 = randomInteger(myMT19937, popsize);
		while ((*r5 == candidate)
		       || (*r5 == *r4)
		       || (*r5 == *r3)
			   || (*r5 == *r2)
			   || (*r5 == *r1));
	}
}


static void Best1Bin(genoptStruct *p, long candidate){
	long r1, r2;
	unsigned long n, i;

	SelectSamples(&(p->myMT19937),
	              p->totalpopsize,
	              candidate,
	              &r1,
	              &r2,
	              NULL,
	              NULL,
	              NULL);
	n = randomInteger(&(p->myMT19937),  p->numvarparams);

	memcpy(p->gen_trial,
	       *(p->gen_populationvector + candidate),
	       p->numvarparams * sizeof(double));

	for (i=0; i < p->numvarparams; i++) {
		if ((randomDouble(&(p->myMT19937), 0, 1) < p->recomb) || (i == (p->numvarparams - 1)))
			*(p->gen_trial + n) = p->gen_populationvector[0][n]
			                        + p->k_m * (p->gen_populationvector[r1][n]
			                        - p->gen_populationvector[r2][n]);

		n = (n + 1) % p->numvarparams;
	}
}

static void Best1Exp(genoptStruct *p, long candidate){
	long r1, r2;
	unsigned long n, i;

	SelectSamples(&(p->myMT19937),  p->totalpopsize, candidate, &r1, &r2, NULL, NULL, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);

	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));

	for (i=0 ; (randomDouble(&(p->myMT19937),  0, 1) < p->recomb) && (i < p->numvarparams); i++){
		*(p->gen_trial + n) = p->gen_populationvector[0][n]
		+ p->k_m * (p->gen_populationvector[r1][n] - p->gen_populationvector[r2][n]);

		n = (n + 1) % p->numvarparams;
	}
}

static void Rand1Exp(genoptStruct *p, long candidate){
	long r1, r2, r3;
	unsigned long n, i;

	SelectSamples(&(p->myMT19937),  p->totalpopsize, candidate,&r1,&r2,&r3, NULL, NULL);
	n = randomInteger(&(p->myMT19937),  p->numvarparams);

	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));

	for (i=0; (randomDouble(&(p->myMT19937), 0, 1) < p->recomb) && (i < p->numvarparams); i++) {
		*(p->gen_trial + n) = p->gen_populationvector[r1][n]
		+ p->k_m * (p->gen_populationvector[r2][n] - p->gen_populationvector[r3][n]);

		n = (n + 1) % p->numvarparams;
	}
}

static void RandToBest1Exp(genoptStruct *p, long candidate){
	long r1, r2;
	unsigned long n,  i;

	SelectSamples(&(p->myMT19937), p->totalpopsize, candidate,&r1,&r2, NULL, NULL, NULL);
	n = randomInteger(&(p->myMT19937), p->numvarparams);

	memcpy(p->gen_trial, *(p->gen_populationvector + candidate), p->numvarparams * sizeof(double));

	for (i=0; (randomDouble(&(p->myMT19937), 0, 1) < p->recomb) && (i < p->numvarparams); i++) {
		*(p->gen_trial + n) += p->k_m * (p->gen_populationvector[0][n] - *(p->gen_trial + n))
		+ p->k_m * (p->gen_populationvector[r1][n]
					   - p->gen_populationvector[r2][n]);
		n = (n + 1) % p->numvarparams;
	}
}

static void Best2Exp(genoptStruct *p, long candidate){
	long r1, r2, r3, r4;
	unsigned long n, i;

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
}

static void Rand2Exp(genoptStruct *p, long candidate){
	long r1, r2, r3, r4, r5;
	unsigned long n, i;

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
}

static void RandToBest1Bin(genoptStruct *p, long candidate){
	long r1, r2;
	unsigned long n, i;

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
}

static void Best2Bin(genoptStruct *p, long candidate){
	long r1, r2, r3, r4;
	unsigned long n, i;

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
}

static void Rand2Bin(genoptStruct *p, long candidate){
	long r1, r2, r3, r4, r5;
	unsigned long n, i;

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
}

static void Rand1Bin(genoptStruct *p, long candidate){
	long r1, r2, r3;
	unsigned long n, i;

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
}

/*
 createTrialVector makes a mutated vector.  It fills the trialVector from the current pvector and from bPrime,
 in modulo.
 bPrime is created from two random population vectors and the best fit vector.
 */
static void createTrialVector(genoptStruct *p, long currentpvector){
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
			break;
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
	waveStats retval;

	memset(&retval, 0, sizeof(retval));

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
 ensureConstraints takes the current trial vector and makes sure that all the
 individual parameters lie inbetween the upper and lower limits.
 returns void
 */
static void
ensureConstraints(genoptStruct *p){
	unsigned int ii;
	for(ii = 0 ; ii < p->numvarparams ; ii++)
		if(p->gen_trial[ii] < 0 || p->gen_trial[ii] > 1)
			p->gen_trial[ii] = randomDouble(&(p->myMT19937), 0, 1);
}

/*
 scale_parameters scales values from a range [0, 1] to a value that
 is between the upper and lower limits
 */
 void
scale_parameters(double *coefs,
                 const unsigned int* varparams,
                 unsigned int numvarparams,
                 const double *scaledVector,
                 const double **scale_factors){
	unsigned int ii;
	int ival;
	double dval,scale0, scale1, dval4;

	for(ii = 0 ; ii < numvarparams ; ii += 1){
		ival = varparams[ii];
		dval = scaledVector[ii];
        scale0 = scale_factors[0][ii];
        scale1 = scale_factors[1][ii];

        dval4 = coefs[varparams[ii]] = scale0 + (dval - 0.5) * scale1;
	}
}

/*
 unscale_parameters scales values from a value that is between the upper and
 lower limits to a value in the range [0, 1]
 */
 void
unscale_parameters(double *scaledVector,
                   const unsigned int* varparams,
                   unsigned int numvarparams,
                   const double *coefs,
                   const double **scale_factors){
	unsigned int ii;
	int ival;
	double dval,scale0, scale1, dval4;

	for(ii = 0 ; ii < numvarparams ; ii += 1){
		ival = varparams[ii];
		dval = coefs[ival];
        scale0 = scale_factors[0][ii];
        scale1 = scale_factors[1][ii];

        dval4 = scaledVector[ii] = (dval - scale0) / scale1 + 0.5;
	}
}

/*
 setPopVector sets the populationvector with index replace, with a double array
 */
static void
setPopVector(genoptStruct *p, double* vector, int vectorsize, long replace){
	memcpy(*(p->gen_populationvector + replace), vector, vectorsize * sizeof(double));
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
 overwrites gen_trial.
 */
static void
swapPopVector(genoptStruct *p, long popsize, long i, long j){
	memcpy(p->gen_trial,
	       *(p->gen_populationvector + j),
	       p->numvarparams * sizeof(double));

	memcpy(*(p->gen_populationvector + j),
	       *(p->gen_populationvector + i),
	       p->numvarparams * sizeof(double));

	memcpy(*(p->gen_populationvector + i),
	       p->gen_trial,
	       p->numvarparams * sizeof(double));
}


/*
initialises the genetic optimisation for the optimiseloop, all memory allocation
has occurred by now, just need to fill out the arrays
*/
static int initialiseFit(genoptStruct *p){
	int err = 0;

	unsigned int ii = 0, startIt = 0;
	double *val = NULL;
	double chi2 = 0.;
	waveStats wavStats;

    /*initialise scale_factors*/
    for(ii = 0 ; ii < p->numvarparams ; ii += 1){
        p->scale_factors[0][ii] = 0.5 * (p->limits[0][p->varparams[ii]] + p->limits[1][p->varparams[ii]]);
        p->scale_factors[1][ii] = fabs(p->limits[0][p->varparams[ii]] - p->limits[1][p->varparams[ii]]);
    }
    
	/*
	 initialise population vector guesses
	 if you want to seed with the initial guesses, you have to insert them into
	 the first row of the population vector.
	 The population vector is scaled between 0 and 1, so you have to rescale the
	 coefficient from a position between the bottom and upper limits.
	 If the limits are the same than the value will be NaN (divide by 0), so set
	 the value to one of the limits.
	 */
    val = *(p->gen_populationvector);
    
	if(p->useinitialguesses){
        /*
        if you want to use the initial guesses to start the fit, then they need
        to be within the limits.
        */

        unscale_parameters(p->gen_trial,
                           p->varparams,
                           p->numvarparams,
                           p->coefs,
                           (const double **) p->scale_factors);
        for(ii = 0 ; ii < p->numvarparams ; ii++){
            if(p->gen_trial[ii] > 1. || p->gen_trial[ii] < 0.){
                err = COEFS_MUST_BE_WITHIN_LIMITS;
                goto done;
            }
			if(!isfinite(p->gen_trial[ii]))
				p->gen_trial[ii] = p->limits[1][p->varparams[ii]];

			*val++ = p->gen_trial[ii];
        }
        startIt = p->numvarparams;
	} else {
		startIt = 0;
    }

	/*
	 initialise population vector with random numbers between 0 and 1.
	 The loop counter, ii, is initialised before here. This is because we may
	 want to use the initial guesses to seed the fit.
	 If so, then it is initialised to p->numvarparams, which should correspond
	 to the second row of the population vector
	 (p->gen_populationvector[p->totalpopsize][p->numvarparams])
	 */
    for(ii = startIt; ii < p->numvarparams * p->totalpopsize ; ii++)
        *val++ = randomDouble(&(p->myMT19937), 0, 1);
    
    if(p->initial_population) {
        /* an initial population was specified, go through and initialise
         the population with it.
         */
        val = (double*) p->initial_population;
        for(ii = 0 ; ii < MIN(p->initial_population_rows, p->totalpopsize) ; ii++){
            unscale_parameters(p->gen_trial,
                               p->varparams,
                               p->numvarparams,
                               val,
                               (const double **) p->scale_factors);
            val += p->numcoefs;
            ensureConstraints(p);
            setPopVector(p, p->gen_trial, p->numvarparams, ii);
        }
    }
    
	/*
	initialise Chi2array, will require a bit of calculation of the model
	function for each of the initial guesses.
	*/
	for(ii = 0 ; ii < p->totalpopsize ; ii += 1){
		scale_parameters(p->temp_coefs,
		                 (const unsigned int*) p->varparams,
		                 p->numvarparams,
		                 (const double*) *(p->gen_populationvector + ii),
		                 (const double**) p->scale_factors);

		//calculate the model
		if((err = (*(p->fitfun))(p->userdata,
		                         p->temp_coefs,
		                         p->numcoefs,
		                         p->model,
		                         p->xdata,
		                         p->datapoints,
		                         p->numDataDims)))
			goto done;

		/*calculate the costfunction*/
		chi2 = (*(p->costfun))(p->userdata,
		                       p->temp_coefs,
		                       p->numcoefs,
		                       p->ydata,
		                       p->model,
		                       p->edata,
		                       p->datapoints);

		*(p->chi2Array + ii)= chi2;
	}
	//find best chi2 and put that into number 0 pos.
	wavStats = getWaveStats(p->chi2Array, p->totalpopsize, 0);

	swapChi2values(p, 0, wavStats.V_minloc);
	swapPopVector(p, p->totalpopsize, 0, wavStats.V_minloc);


	//put the best fit from the intialisation into the coeffcients to return
	scale_parameters(p->temp_coefs,
	                 p->varparams,
	                 p->numvarparams,
	                 *(p->gen_populationvector),
	                 (const double **)p->scale_factors);

	if(p->updatefun && (4 & p->updatefrequency))
		if((err = (*(p->updatefun))(p->userdata,
									p->temp_coefs,
									p->numcoefs,
									0,
									*(p->chi2Array),
									4,
									-1,
									(const double**) p->gen_populationvector,
									p->varparams,
									p->numvarparams,
									p->popsizeMultiplier * p->numvarparams,
									(const double *)p->chi2Array)))
			goto done;

	memcpy(p->coefs, p->temp_coefs, p->numcoefs * sizeof(double));

done:
	return err;
}

/*
 optimiseloop performs the optimisation.  It takes the initial population and
 mutates it until we find the best fit solution.
 returns 0 if no errors
 returns errorcode otherwise.
 */
static int optimiseloop(genoptStruct *p){
	int err = 0;
	unsigned int ii, kk, currentpvector;
	double chi2pvector,chi2trial;
	waveStats wavStats;
	double convergenceNumber = -1;

	/* the user sets how many times through the entire population*/
	for(kk = 1; kk <= p->iterations ; kk += 1){
		p->numfititers = kk;

		if(p->dither[0] > 0. && p->dither[1] > 0.)
		    p->k_m = randomDouble(&(p->myMT19937), p->dither[0], p->dither[1]);

		/*
		 if the SD of the population divided by it's average is less than
		 tolerance stop.
		 */
		wavStats = getWaveStats(p->chi2Array, p->totalpopsize, 1);
		convergenceNumber = wavStats.V_avg * p->tolerance / wavStats.V_stdev;
		if(convergenceNumber > 1){
			if(p->updatefun && (16 & p->updatefrequency))
				/*
				 send the best fit to the update function
				 */
				if((err = (*(p->updatefun))(p->userdata,
                                    p->coefs,
                                    p->numcoefs,
                                    kk,
                                    *(p->chi2Array),
                                    16,
                                    convergenceNumber,
                                    (const double **) p->gen_populationvector,
                                    p->varparams,
                                    p->numvarparams,
                                    p->popsizeMultiplier * p->numvarparams,
                                    (const double*)p->chi2Array)))
					goto done;
			goto done;
		}

		if(p->updatefun && (8 & p->updatefrequency)){
			/*
			 send the best fit to the update function
			 */
			if((err = (*(p->updatefun))(p->userdata,
                                    p->coefs,
                                    p->numcoefs,
                                    kk,
                                    *(p->chi2Array),
                                    8,
                                    convergenceNumber,
                                    (const double **) p->gen_populationvector,
                                    p->varparams,
                                    p->numvarparams,
                                    p->popsizeMultiplier * p->numvarparams,
                                    (const double*)p->chi2Array)))
				goto done;
		}

		/*iterate over all the individual members of the population*/
		for(ii = 0 ; ii < p->totalpopsize ; ii += 1){
			currentpvector = ii;
			/*
			 now set up the trial vector using a wave from the populationvector
			 and bprime first set the pvector and create a mutated trial vector
			 from the best fit and two random population members
			 */
			createTrialVector(p, currentpvector);

			/*
			 make sure the trial vector parameters lie between the user defined
			 limits
			 */
			ensureConstraints(p);

			chi2pvector = *(p->chi2Array + ii);

			/*
			 find out the chi2 value of the trial vector
			 */
			scale_parameters(p->temp_coefs,
			                 p->varparams,
			                 p->numvarparams,
			                 (const double*) p->gen_trial,
			                 (const double**) p->scale_factors);

			if((err = (*(p->fitfun))(p->userdata,
									 p->temp_coefs,
									 p->numcoefs,
									 p->model,
									 p->xdata,
									 p->datapoints,
									 p->numDataDims)))
				goto done;

			/*calculate the costfunction*/
			chi2trial = (*(p->costfun))(p->userdata,
										p->temp_coefs,
										p->numcoefs,
										p->ydata,
										p->model,
										p->edata,
										p->datapoints);

			/*
			 if the chi2 of the trial vector is less than the current
			 populationvector then replace it
			 */
			if(chi2trial < chi2pvector){
				setPopVector(p, p->gen_trial, p->numvarparams, currentpvector);
				*(p->chi2Array + ii) = chi2trial;

				/*
				 if chi2 of the trial vector is less than that of the best fit,
				 then replace the best fit vector
				 */
				if(chi2trial < *(p->chi2Array)){
				/*
				if this trial vector is better than the current best then
				replace it
				*/
					setPopVector(p, p->gen_trial, p->numvarparams, 0);

					/*
					 update the best fit vector (p->temp_coefs should hold the
					 best vector still)
					 */
					memcpy(p->coefs,
					       p->temp_coefs,
					       p->numcoefs * sizeof(double));

					/*
					 if the fractional decrease in chi2 is less than the
					 tolerance then abort the fit
					 */
					wavStats = getWaveStats(p->chi2Array, p->totalpopsize, 1);

					/*
					 update the best chi2 if you've just found a better fit (but
					 not yet reached termination
					 */
					*(p->chi2Array) = chi2trial;

					/*
					 a user defined update function that can be used to halt the
					 fit early, and keep appraised of fit progress.
					 */
					if(p->updatefun && (1 & p->updatefrequency))
						if((err = (*(p->updatefun))(p->userdata,
                                    p->coefs,
                                    p->numcoefs,
                                    kk,
                                    *(p->chi2Array),
                                    1,
                                    convergenceNumber,
                                    (const double **) p->gen_populationvector,
                                    p->varparams,
                                    p->numvarparams,
                                    p->popsizeMultiplier * p->numvarparams,
                                    (const double*)	p->chi2Array)))
							goto done;

					/*
					 if the SD of the population divided by it's average is less
					 than tolerance stop.
					 */
					convergenceNumber = (wavStats.V_avg * p->tolerance
					                      / wavStats.V_stdev);
					if(convergenceNumber > 1){
						if(p->updatefun && (16 & p->updatefrequency))
							if((err = (*(p->updatefun))(p->userdata,
                                    p->coefs,
                                    p->numcoefs,
                                    kk,
                                    *(p->chi2Array),
                                    16,
                                    convergenceNumber,
                                    (const double **) p->gen_populationvector,
                                    p->varparams,
                                    p->numvarparams,
                                    p->popsizeMultiplier * p->numvarparams,
                                    (const double*)	p->chi2Array)))
								goto done;
						goto done;
					}
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

	long ii, jj;
	unsigned int kk;
	double *yyMC = NULL;
	/*the overall structure to contain the entire fit*/
	genoptStruct gos;

	//fit function must exist
	if(!fitfun)
		return NO_FIT_FUNCTION_SPECIFIED;

	//y, x, e arrays must exist
	if(!ydata)
		return NO_Y_ARRAY;

	if(!xdata)
		return NO_X_ARRAY;

	if(!edata)
		return NO_E_ARRAY;

	if(!coefs)
		return NO_COEFS_ARRAY;

	if(!limits)
		return NO_LIMITS_ARRAY;

    // the options structure is zeroed out to start with
	memset(&gos, 0, sizeof(gos));

	//initialise the random number generators
	if(gco == NULL || gco->seed < 1)
		sgenrand(clock(), &(gos.myMT19937));
	else
		sgenrand(gco->seed, &gos.myMT19937);

	//setup the data.
	//the user may wish to alter the data to do a monte carlo iteration
	if(gco && gco->monteCarlo){
		yyMC = (double*) malloc(sizeof(double) * datapoints);
		if(!yyMC){
			err = NO_MEMORY;
			goto done;
		}
		for( ii = 0 ; ii < datapoints ; ii++)
			yyMC[ii] = ydata[ii] + gnoise(&(gos.myMT19937), edata[ii]);

		gos.ydata = yyMC;
	} else {
		gos.ydata = ydata;
	}

	gos.xdata = xdata;
	gos.edata = edata;
	gos.datapoints = datapoints;
	gos.numDataDims = numDataDims;

	//setup the function pointers
	gos.fitfun = fitfun;
	if(!costfun)
		gos.costfun = &chisquared;
	else
		gos.costfun = costfun;

	//check that the total number of parameters matches the number of parameters
	// in the holdvector
	gos.numcoefs = numcoefs;
	gos.coefs = coefs;
	gos.holdvector = holdvector;

	//add in userdata
	gos.userdata = userdata;

	//work out which parameters are being held
	for(kk = 0 ; kk < numcoefs ; kk += 1)
		if(holdvector[kk] == 0)
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
	for(kk = 0 ; kk < numcoefs ; kk += 1){
		if(holdvector[kk] == 0){
		    gos.varparams[jj] = kk;
			jj+=1;
		}
	}

	//optional parameters
	if(!gco){
		gos.updatefun = NULL;
		popsizeMultiplier = 20;
		gos.k_m = 0.7;
		gos.recomb = 0.5;
		gos.popsizeMultiplier = 20;
		gos.iterations = 200;
		gos.tolerance = 0.02;
		gos.updatefrequency = 1;
		gos.useinitialguesses = 0;
		gos.dither[0] = -1.;
		gos.dither[1] = -1.;
	} else {
		gos.updatefun = gco->updatefun;
		popsizeMultiplier = gco->popsizeMultiplier;
		gos.popsizeMultiplier = gco->popsizeMultiplier;
		gos.k_m = gco->k_m;
		gos.recomb = gco->recomb;
		gos.iterations = gco->iterations;
		gos.tolerance = gco->tolerance;
		gos.updatefrequency = gco->updatefrequency;
		gos.useinitialguesses = gco->useinitialguesses;
		gos.dither[0] = gco->dither[0];
		gos.dither[1] = gco->dither[1];
        gos.initial_population = gco->initial_population;
        gos.initial_population_rows = gco->initial_population_rows;
	}
	gos.totalpopsize = gos.numvarparams * popsizeMultiplier;

	/*
	 now we can allocate memory for the rest of the structure members in gos.
	initialise population vector
	 */
	gos.gen_populationvector = (double**)malloc2d(gos.totalpopsize,
	                                              gos.numvarparams,
	                                              sizeof(double));
	if(gos.gen_populationvector == NULL){
		err = NO_MEMORY;
		goto done;
	}
    
    /*
	 initialise Chi2array
	 */
	gos.chi2Array = (double*) malloc (gos.totalpopsize * sizeof(double));
	if(gos.chi2Array == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memset(gos.chi2Array, 0, gos.totalpopsize);

	/*
	 initialise the trial vector
	 */
	gos.gen_trial = (double*)malloc(gos.numvarparams * sizeof(double));
	if(gos.gen_trial == NULL){
		err = NO_MEMORY;
		goto done;
	}
    
    /* initialise the scale factors */
    gos.scale_factors = (double**)malloc2d(2, gos.numvarparams,
                                           sizeof(double));
    if(gos.scale_factors == NULL){
        err = NO_MEMORY;
        goto done;
    }

	/*
	initialise space for a full array copy of the coefficients
	*/
	gos.temp_coefs = (double*)malloc(numcoefs * sizeof(double));
	if(gos.temp_coefs == NULL){
		err = NO_MEMORY;
		goto done;
	}
	memcpy(gos.temp_coefs, coefs, numcoefs * sizeof(double));

	/*
	 initialise space for a full array copy of the coefficients
	 */
	gos.model = (double*)malloc(datapoints* sizeof(double));
	if(gos.model == NULL){
		err = NO_MEMORY;
		goto done;
	}

	gos.limits = limits;

	/*
	 at this stage we can initialise the fit, before looping.
	 */
	if((err = initialiseFit(&gos)))
		goto done;

	/*
	 now iterate through, fitting the data.
	 */
	err = optimiseloop(&gos);

	if(!err && gco->polishUsingLM)
		err = levenberg_marquardt(fitfun,
								  costfun,
								  numcoefs,
								  gos.temp_coefs,
								  holdvector,
								  datapoints,
								  ydata,
								  xdata,
								  edata,
								  numDataDims,
								  chi2,
								  gco,
								  userdata);

	/*
	 at this point we have the best fit.  Now return the results
	 */
	if(chi2)
		*chi2 = *(gos.chi2Array);

	/*
	 put the best fit into the coeffcient array that will be returned to the user
	 */
	if(!err){
		scale_parameters(gos.temp_coefs,
		                 gos.varparams,
		                 gos.numvarparams,
		                 (const double*) *(gos.gen_populationvector),
		                 (const double**) gos.scale_factors);
		memcpy(gos.coefs, gos.temp_coefs, gos.numcoefs * sizeof(double));
	}

done:
	if(yyMC)
		free(yyMC);
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
    if(gos.scale_factors)
        free(gos.scale_factors);

	return err;
}

double chisquared(void *userdata,
                  const double *params,
                  unsigned int numcoefs,
                  const double *data,
                  const double *model,
                  const double *errors,
                  long datapoints){
    long ii;
    double chi2 = 0;
    double val=0;
    for (ii = 0; ii < datapoints ; ii += 1){
        val = pow((fabs((data[ii] - model[ii])/errors[ii])), 2);
        if(isfinite(val))
            chi2 += val;
    }
    return chi2;
}

double robust(void *userdata,
              const double *params,
              unsigned int numcoefs,
              const double *data,
              const double *model,
              const double *errors,
              long datapoints){
     long ii;
     double chi = 0;
     double val=0;
     for (ii = 0; ii < datapoints ; ii += 1){
         val = fabs((data[ii] - model[ii])/errors[ii]);
         if(isfinite(val))
             chi += val;
     }
     return chi;
 }

