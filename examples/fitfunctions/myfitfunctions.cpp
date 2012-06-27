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
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define NUM_CPUS 1


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
abelescalcall(void *userdata, const double *coefP, double *yP, const double *xP,long npoints, int Vmullayers, int Vmulappend, int Vmulrep){
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
		answer = ((num / den) * scale) + bkg;
		
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

int smearedabeles(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	int ii;
	double *dyP;
	double *ddxP = NULL;
	double *yP = model;
	const double *xP = *xdata;
	const double *dxP = *(xdata + 1);

	int RESPOINTS = 13;
	double INTLIMIT = 3.5;		//integration between -3.5 and 3 sigma
	double FWHM = 2 * sqrt(2 * log(2.0));
	double va, vb, sigma;
	
	double weights[13] = {0.0404840047653159,
		0.0921214998377285,
		0.1388735102197872,
		0.1781459807619457,
		0.2078160475368885,
		0.2262831802628972,
		0.2325515532308739,
		0.2262831802628972,
		0.2078160475368885,
		0.1781459807619457,
		0.1388735102197872,
		0.0921214998377285,
		0.0404840047653159};
	double abscissa[13] = {-0.9841830547185881,
		-0.9175983992229779,
		-0.8015780907333099,
		-0.6423493394403402,
		-0.4484927510364469,
		-0.2304583159551348,
		0.,
		0.2304583159551348,
		0.4484927510364469,
		0.6423493394403402,
		0.8015780907333099,
		0.9175983992229779,
		0.9841830547185881};
	
	
	dyP = (double*)malloc(numpnts * sizeof(double));
	if(!dyP)
		err = NO_MEMORY;

	ddxP = (double*)malloc(numpnts * RESPOINTS * sizeof(double));
	if(!ddxP)
		err = NO_MEMORY;
	
	for(ii = 0 ; ii < numpnts * RESPOINTS ; ii ++){
		sigma = dxP[ii/RESPOINTS] / FWHM;
		va = -INTLIMIT * sigma + xP[ii / RESPOINTS];
		vb = INTLIMIT * sigma + xP[ii / RESPOINTS];
		ddxP[ii] = (abscissa[ii % RESPOINTS] * (vb-va) + vb + va)/2;
	}	

	if(err = abeles(userdata, coefs, numcoefs, dyP, (const double**) &ddxP, numpnts * RESPOINTS, 1))
		goto done;
	
	for(ii = 0, yP = model ; ii < numpnts ; ii ++, yP++){
		//assumes 13 point gaussian quadrature, over +/- 3.5 sigma
		*yP = dyP[ii * RESPOINTS] * weights[0] * 0.001057642102668805;
		*yP += dyP[ii * RESPOINTS + 1] * weights[1] * 0.002297100003792314;
		*yP += dyP[ii * RESPOINTS + 2] * weights[2] *0.007793859679303332;
		*yP += dyP[ii * RESPOINTS + 3] * weights[3] *0.0318667809686739;
		*yP += dyP[ii * RESPOINTS + 4] * weights[4] *0.1163728244269813;
		*yP += dyP[ii * RESPOINTS + 5] * weights[5] *0.288158781825899;
		*yP += dyP[ii * RESPOINTS + 6] * weights[6] * 0.3989422804014327;
		*yP += dyP[ii * RESPOINTS + 7] * weights[7] * 0.288158781825899;
		*yP += dyP[ii * RESPOINTS + 8] * weights[8] * 0.1163728244269813;
		*yP += dyP[ii * RESPOINTS + 9] * weights[9] * 0.0318667809686739;
		*yP += dyP[ii * RESPOINTS + 10] * weights[10] * 0.007793859679303332;
		*yP += dyP[ii * RESPOINTS + 11] * weights[11] *0.002297100003792314;
		*yP += dyP[ii * RESPOINTS + 12] * weights[12] * 0.001057642102668805;
		
		*yP *= 3.5;
	}
		
done:
	if(dyP)
		free(dyP);
	if(ddxP)
		free(ddxP);
	
	return err;
};


void *AbelesThreadWorker(void *arg){
	int err = 0;
	refCalcParm *p = (refCalcParm *) arg;
	err = abelescalcall(p->userdata, p->coefP, p->yP, p->xP, p->npoints, p->Vmullayers, p->Vappendlayer, p->Vmulrep);
	
	pthread_exit((void*)err);
	return NULL;
}

int abeles(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	pthread_t *threads;
//	extern int NUM_CPUS;
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
	
	err = abelescalcall(userdata, coefs, model + pointsConsumed, (*xdata) + pointsConsumed, pointsRemaining, 0, 0, 0);
	
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
int abelesmodelwrapper(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
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
	err = abeles(userdata, &coefstemp[0], numcoefs, model, xdata, numpnts, numDataDims);
	return err;
}




/*
  
 SAH SAMfloat_monolayer function
 
	w[0] = scale
	w[1] = SLD fronting
	w[2] = SLD backing
	w[3] = bkg
	w[4]=backing rough

	w[5]=oxide thickness
	w[6]=oxide SLD
	w[7]=oxide solvent
	w[8]=Si rough

	w[9]=TiO2 thickness
	w[10]=TiO2 SLD
	w[11]=TiO2 solvent
	w[12]=oxide rough

	w[13] = SAM A per mol
	w[14] = HG thickness
	w[15]= TiO2 roughness
	w[16] = Tail thickness
	w[17] = roughness of HG

 */
int stephenssamfloat_monolayer(void *userdata, const double *coefs, unsigned int numcoefs, double *model, const double **xdata, long numpnts, unsigned int numDataDims){
	int err = 0;
	
	vector<double> W_forreflectivity;
	
	double Vhead = 60.807;
	double Vtail = 231.82;
	double btail = -0.0001331;
	double bhead = 0.0002835;

	//#layers
	W_forreflectivity.push_back(4);
	//scale
	W_forreflectivity.push_back(coefs[0]);
	//fronting
	W_forreflectivity.push_back(coefs[1]);
	//backing
	W_forreflectivity.push_back(coefs[2]);
	//Bgd
	W_forreflectivity.push_back(coefs[3]);
	//backing rough
	W_forreflectivity.push_back(coefs[4]);
	
	//SiO2 Layer
	W_forreflectivity.push_back(coefs[5]);
	W_forreflectivity.push_back(coefs[6]);
	W_forreflectivity.push_back(coefs[7]);
	//Si/SiO2 roughness
	W_forreflectivity.push_back(coefs[8]);
	
	//TiO2
	W_forreflectivity.push_back(coefs[9]);
	W_forreflectivity.push_back(coefs[10]);
	W_forreflectivity.push_back(coefs[11]);
	//TiO2/SiO2 roughness
	W_forreflectivity.push_back(coefs[12]);
	
	//SAM HG
	W_forreflectivity.push_back(coefs[14]);
	W_forreflectivity.push_back(bhead  / (coefs[14] * coefs[13]) + (1 - Vhead / (coefs[14] * coefs[13])) * coefs[2]);
	W_forreflectivity.push_back(0);
	//SAM/TiO2 roughness
	W_forreflectivity.push_back(coefs[15]);
	
	
	//SAM Tail
	W_forreflectivity.push_back(coefs[16]);
	W_forreflectivity.push_back(btail  / (coefs[16] * coefs[13]) + (1 - Vtail / (coefs[16] * coefs[13])) * coefs[2]);
	W_forreflectivity.push_back(0);
	//tail/HG rough
	W_forreflectivity.push_back(coefs[17]);

	err = smearedabeles(userdata, (const double*) &W_forreflectivity[0], W_forreflectivity.size(), model, xdata, numpnts, numDataDims);
	return err;
	
}


/**
 regularising cost function for reflectivity
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


/**
 a log10 cost function for reflectivity
 */
double log10chisquared(void *userdata, const double *params, unsigned int numparams, const double *data, const double *model, const double *errors, long numpnts){
	
	long ii;
	double chi2 = 0;
	double val=0;
	double ln10 = log(10);
	
	for (ii = 0 ; ii < numpnts ; ii += 1){
		val = log10(data[ii]) - log10(model[ii]);
		val /= fabs(errors[ii] / data[ii] / ln10);
		val = pow(val, 2);		
		if(isfinite(val))
			chi2 += val;
	}
	
	return chi2;
	
}
