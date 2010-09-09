/*
 *  myfitfunctions.cpp
 *  motoMC
 *
 *  Created by andrew on 29/05/09.
 *  Copyright 2009 __MyCompanyName__. All rights reserved.
 *
 */

#include "myfitfunctions.h"
#include "gencurvefit.h"
#include <math.h>
#include "MyComplex.h"
#include "pthread.h"
#include <vector>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>


using namespace std;
using namespace MyComplexNumber;

int line(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	long ii;
		
	for(ii=0 ; ii<numpnts ; ii+=1){
	    *(model + ii) = *(coefs) + (*xdata)[ii] * (*(coefs+1));
	}
	
	return err;
}

int gaussian(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	long ii;
	
	for(ii = 0 ; ii < numpnts ; ii += 1)
	    model[ii] = coefs[0] + coefs[1] * exp(-1 * pow((coefs[2] - (*xdata)[ii])/coefs[3], 2));

	return err;
}

MyComplex fres(MyComplex a,MyComplex b,double rough){
	return (compexp(-2*rough*rough*a*b))*(a-b)/(a+b);
}

int 
AbelesCalcAll(void *userdata, const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep){
	int err = 0;
	int j;
	
	int ii=0,jj=0,kk=0;
	
	double scale,bkg,subrough;
	double num=0,den=0, answer=0,qq;
	double anum,anum2;
	MyComplex temp,SLD,beta,rj;
	double numtemp=0;
	int offset=0;
	MyComplex  MRtotal[2][2];
	MyComplex subtotal[2][2];
	MyComplex MI[2][2];
	MyComplex temp2[2][2];
	MyComplex qq2;
	MyComplex oneC = MyComplex(1,0);
	MyComplex *pj_mul = NULL;
	MyComplex *pj = NULL;
	double *SLDmatrix = NULL;
	double *SLDmatrixREP = NULL;
	
	int nlayers = (int)coefP[0];
	
	try{
		pj = new MyComplex [nlayers+2];
		SLDmatrix = new double [nlayers+2];
	} catch(...){
		err = NO_MEMORY;
		goto done;
	}
	
	memset(pj, 0, sizeof(pj));
	memset(SLDmatrix, 0, sizeof(SLDmatrix));
	
	scale = coefP[1];
	bkg = fabs(coefP[4]);
	subrough = coefP[5];
	
	//offset tells us where the multilayers start.
	offset = 4 * nlayers + 6;
	
	//fillout all the SLD's for all the layers
	for(ii=1; ii<nlayers+1;ii+=1){
		numtemp = 1.e-6 * ((100. - coefP[4*ii+4])/100.) * coefP[4*ii+3]+ (coefP[4*ii+4]*coefP[3]*1.e-6)/100.;		//sld of the layer
		
		*(SLDmatrix+ii) = 4*PI*(numtemp  - (coefP[2]*1e-6));
	}
	*(SLDmatrix) = 0;
	*(SLDmatrix+nlayers+1) = 4*PI*((coefP[3]*1e-6) - (coefP[2]*1e-6));
	
	
	if(Vmullayers > 0 && Vmulrep > 0 && Vmulappend >= 0){
		//set up an array for wavevectors
		try{
			SLDmatrixREP = new double [Vmullayers];
			pj_mul = new MyComplex [Vmullayers];
		} catch(...){
			err = NO_MEMORY;
			goto done;
		}
		memset(pj_mul, 0, sizeof(pj_mul));
		for(ii=0; ii<Vmullayers;ii+=1){
			numtemp = (coefP[3]*1e-6*coefP[(4*ii)+offset+2]/100) +(1e-6 * ((100 - coefP[(4*ii)+offset+2])/100) * coefP[(4*ii)+offset+1]);		//sld of the layer
			*(SLDmatrixREP+ii) = 4*PI*(numtemp  - (coefP[2]*1e-6));
		}
	}
	
	for (j = 0; j < npoints; j++) {
		//intialise the matrices
		memset(MRtotal,0,sizeof(MRtotal));
		MRtotal[0][0] = oneC ; MRtotal[1][1] = oneC;
		
		qq = xP[j]*xP[j]/4;
		qq2=MyComplex(qq,0);
		
		for(ii=0; ii<nlayers+2 ; ii++){			//work out the wavevector in each of the layers
			pj[ii] = (*(SLDmatrix+ii)>qq) ? compsqrt(qq2-MyComplex(*(SLDmatrix+ii),0)): MyComplex(sqrt(qq-*(SLDmatrix+ii)),0);
		}
		
		//workout the wavevector in the toplayer of the multilayer, if it exists.
		if(Vmullayers>0 && Vmulrep > 0 && Vmulappend >=0){
			memset(subtotal,0,sizeof(subtotal));
			subtotal[0][0]=MyComplex(1,0);subtotal[1][1]=MyComplex(1,0);
			pj_mul[0] = (*(SLDmatrixREP)>qq) ? compsqrt(qq2-MyComplex(*SLDmatrixREP, 0)) : MyComplex(sqrt(qq-*SLDmatrixREP),0);
		}
		
		//now calculate reflectivities
		for(ii = 0 ; ii < nlayers+1 ; ii++){
			//work out the fresnel coefficients
			//this looks more complicated than it really is.
			//the reason it looks so convoluted is because if there is no complex part of the wavevector,
			//then it is faster to do the calc with real arithmetic then put it into a complex number.
			if(Vmullayers>0 && ii==Vmulappend && Vmulrep>0 ){
				rj=fres(pj[ii],pj_mul[0],coefP[offset+3]);
			} else {
				if((pj[ii]).im == 0 && (pj[ii+1]).im==0){
					anum = (pj[ii]).re;
					anum2 = (pj[ii+1]).re;
					rj.re = (ii==nlayers) ?
					((anum-anum2)/(anum+anum2))*exp(anum*anum2*-2*subrough*subrough)
					:
					((anum-anum2)/(anum+anum2))*exp(anum*anum2*-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5]);
					rj.im = 0;
				} else {
					rj = (ii==nlayers) ?
					((pj[ii]-pj[ii+1])/(pj[ii]+pj[ii+1]))*compexp(pj[ii]*pj[ii+1]*MyComplex(-2*subrough*subrough, 0))
					:
					rj = ((pj[ii]-pj[ii+1])/(pj[ii]+pj[ii+1]))*compexp(pj[ii]*pj[ii+1]*MyComplex(-2*coefP[4*(ii+1)+5]*coefP[4*(ii+1)+5], 0));	
				};
			}
			
			//work out the beta for the (non-multi)layer
			temp.im = fabs(coefP[4*ii+2]);
			temp.re = 0;
			beta = (ii==0)? oneC : compexp(pj[ii] * temp);
			
			//this is the characteristic matrix of a layer
			MI[0][0]=beta;
			MI[0][1]=rj*beta;
			MI[1][1]=oneC/beta;
			MI[1][0]=rj*MI[1][1];
			
			memcpy(temp2, MRtotal, sizeof(MRtotal));
			
			//multiply MR,MI to get the updated total matrix.			
			matmul(temp2,MI,MRtotal);
			
			if(Vmullayers > 0 && ii == Vmulappend && Vmulrep > 0){
				//workout the wavevectors in each of the layers
				for(jj=1 ; jj < Vmullayers; jj++){
					pj_mul[jj] = (*(SLDmatrixREP+jj)>qq) ? compsqrt(qq2-MyComplex(*(SLDmatrixREP+jj),0)): MyComplex(sqrt(qq-*(SLDmatrixREP+jj)),0);
				}
				
				//work out the fresnel coefficients
				for(jj = 0 ; jj < Vmullayers; jj++){
					
					rj = (jj == Vmullayers-1) ?
					//if you're in the last layer then the roughness is the roughness of the top
					((pj_mul[jj]-pj_mul[0])/(pj_mul[jj]+pj_mul[0]))*compexp((pj_mul[jj]*pj_mul[0])*MyComplex(-2*coefP[offset+3]*coefP[offset+3],0))
					:
					//otherwise it's the roughness of the layer below
					((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*MyComplex(-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3],0));
					
					//Beta's
					beta = compexp(MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]);
					
					MI[0][0]=beta;
					MI[0][1]=rj*beta;
					MI[1][1]=oneC/beta;
					MI[1][0]=rj*MI[1][1];
					
					memcpy(temp2, subtotal, sizeof(subtotal));
					//				temp2[0][0] = subtotal[0][0];
					//				temp2[0][1] = subtotal[0][1];
					//				temp2[1][0] = subtotal[1][0];
					//				temp2[1][1] = subtotal[1][1];
					
					matmul(temp2,MI,subtotal);
				};
				
				for(kk = 0; kk < Vmulrep; kk++){		//if you are in the last multilayer
					if(kk==Vmulrep-1){					//if you are in the last layer of the multilayer
						for(jj=0;jj<Vmullayers;jj++){
							beta = compexp((MyComplex(0,fabs(coefP[4*jj+offset]))*pj_mul[jj]));
							
							if(jj==Vmullayers-1){
								if(Vmulappend==nlayers){
									rj = ((pj_mul[Vmullayers-1]-pj[nlayers+1])/(pj_mul[Vmullayers-1]+pj[nlayers+1]))*compexp((pj_mul[Vmullayers-1]*pj[nlayers+1])*MyComplex(-2*subrough*subrough,0));
								} else {
									rj = ((pj_mul[Vmullayers-1]-pj[Vmulappend+1])/(pj_mul[Vmullayers-1]+pj[Vmulappend+1]))*compexp((pj_mul[Vmullayers-1]*pj[Vmulappend+1])*MyComplex(-2*coefP[4*(Vmulappend+1)+5]*coefP[4*(Vmulappend+1)+5],0));
								};
							} else {
								rj = ((pj_mul[jj]-pj_mul[jj+1])/(pj_mul[jj]+pj_mul[jj+1]))*compexp((pj_mul[jj]*pj_mul[jj+1])*MyComplex(-2*coefP[4*(jj+1)+offset+3]*coefP[4*(jj+1)+offset+3],0));
							}
							
							MI[0][0]=beta;
							MI[0][1]=(rj*beta);
							MI[1][1]=oneC/beta;
							MI[1][0]=(rj*MI[1][1]);
							
							memcpy(temp2, MRtotal, sizeof(MRtotal));							
							matmul(temp2,MI,MRtotal);
						}
					} else {
						memcpy(temp2, MRtotal, sizeof(MRtotal));
						matmul(temp2,subtotal,MRtotal);
					};
				};
			};
			
		}
		
		den=compnorm(MRtotal[0][0]);
		num=compnorm(MRtotal[1][0]);
		answer=log10(((num/den)*scale)+bkg);
		
		*yP++ = answer;
	}
	
done:
	if(pj != NULL)
		delete [] pj;
	if(pj_mul !=NULL)
		delete[] pj_mul;
	if(SLDmatrix != NULL)
		delete[] SLDmatrix;
	if(SLDmatrixREP != NULL)
		delete[] SLDmatrixREP;
	
	return err;
}

typedef struct{
	long npoints;
	int Vmullayers;
	int Vappendlayer;
	int Vmulrep;
	const double *coefP;
	double *yP;
	const double *xP;
	void *userdata;
}  refCalcParm;

void *AbelesThreadWorker(void *arg){
	int err = 0;
	refCalcParm *p = (refCalcParm *) arg;
	err = AbelesCalcAll(p->userdata, p->coefP, p->yP, p->xP, p->npoints, p->Vmullayers, p->Vappendlayer, p->Vmulrep);
	
	pthread_exit((void*)err);
	return NULL;
}

int Abeles(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	pthread_t *threads;
	extern int NUM_CPUS;
	refCalcParm *arg;
	int ii = 0;
	long pointsEachThread = 0, pointsRemaining = 0, pointsConsumed = 0;
	
	threads = (pthread_t *) malloc((NUM_CPUS-1) * sizeof(pthread_t));
	if(!threads && NUM_CPUS > 1){
		err = NO_MEMORY;
		goto done;
	}
	
	arg=(refCalcParm *)malloc(sizeof(refCalcParm)*(NUM_CPUS-1));
	if(!arg && NUM_CPUS > 1){
		err = NO_MEMORY;
		goto done;	
	}
	
	pointsEachThread = (long)floorl(numpnts / NUM_CPUS);
	pointsRemaining = numpnts;
	
	for (ii = 0; ii < NUM_CPUS-1; ii++){
		arg[ii].coefP = coefs;
		arg[ii].npoints = pointsEachThread;
		arg[ii].Vmullayers = 0;
		arg[ii].Vappendlayer = 0;
		arg[ii].Vmulrep = 0;
		arg[ii].xP = (*xdata) + pointsConsumed;
		arg[ii].yP = model + pointsConsumed;
		pthread_create(&threads[ii], NULL, AbelesThreadWorker, (void *)(arg+ii));
		pointsRemaining -= pointsEachThread;
		pointsConsumed += pointsEachThread;
	}
	
	err = AbelesCalcAll(userdata, coefs, model + pointsConsumed, (*xdata) + pointsConsumed, pointsRemaining, 0, 0, 0);
	
	for (ii = 0; ii < NUM_CPUS - 1; ii++)
		pthread_join(threads[ii], NULL);
	
done:
	if(threads)
		free(threads);
	if(arg)
		free(arg);
		
	return err;
};

//this function is meant as a wrapper to allow analytic profiles.
int AbelesModelWrapper(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	vector<double> coefstemp;
	int ii;
	
	coefstemp.push_back(coefs[0]);
	coefstemp.push_back(coefs[1]);
	coefstemp.push_back(coefs[2]);
	coefstemp.push_back(coefs[3]);
	coefstemp.push_back(coefs[4]);
	coefstemp.push_back(coefs[5]);

/*	for(ii=0 ; ii<(int)coefs[0] ; ii+=1){
		coefstemp.push_back(coefs[6]);
		coefstemp.push_back(coefs[ii+8]);
		coefstemp.push_back(0);
		coefstemp.push_back(coefs[7]);
	}
*/

//this code has a different roughness on the top
	for(ii=0 ; ii < (int)coefs[0] ; ii+=1){
		coefstemp.push_back(coefs[6]);
		coefstemp.push_back(coefs[ii+9]);
		coefstemp.push_back(0);
		if(ii)
		    coefstemp.push_back(coefs[7]);
		else
		    coefstemp.push_back(coefs[8]);
	}
	
	//the following code would add an extra layer on top of the substrate.
/*	for(ii=0 ; ii<(int)coefs[0]-1 ; ii+=1){
		coefstemp.push_back(coefs[6]);
		coefstemp.push_back(coefs[ii+12]);
		coefstemp.push_back(0);
		if(ii)
		    coefstemp.push_back(coefs[7]);
		else
		    coefstemp.push_back(coefs[11]);
	}

	coefstemp.push_back(coefs[8]);
	coefstemp.push_back(coefs[9]);
	coefstemp.push_back(0);
	coefstemp.push_back(coefs[10]);
*/
	err = Abeles(userdata, &coefstemp[0], numcoefs, model, xdata, numpnts, numDataDims);
	return err;
}

/*
 cost functions
 */
double smoother(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts){
	
	long ii;
	double chi2 = 0;
	double val=0;
	double lambda = *((double*)userdata);
	double beta = 0;
	int nlayers = (int) params[0];
	
	for (ii=0; ii<numpnts; ii+=1){
		val = pow((fabs((data[ii] - model[ii])/errors[ii])),2);
		if(isfinite(val))
			chi2 += val;
	}
	
	for(ii = 9 ; ii < 9 + nlayers - 1 ; ii += 1){
		beta += pow(params[ii] - params[ii + 1], 2);
	}
	beta += pow(params[2] - params[9], 2);
	beta += pow(params[3] - params[9 + nlayers - 1], 2);
	
	return chi2 + lambda * beta;
}



