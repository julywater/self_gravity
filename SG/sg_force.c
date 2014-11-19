#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<fftw3.h>
#include "../Headers/SG.h"
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/header.h"
//phi force
void cal_force(struct poisson* thePoisson,int DIRECTION){
	int i,j,k;
	int N_p=thePoisson->N_p;
	int N_r=thePoisson->N_r;
	int N_z=thePoisson->N_z;
	double * Force=thePoisson->buffer;
	double * Potential=thePoisson->density;
	double dphi=2*M_PI/N_p;
        double dr=thePoisson->dr;
	double dz=thePoisson->dz;	
	if(DIRECTION==PHI_DIR)
	for(k=0;k<N_z;k++)
		for(i=0;i<N_r;i++){
			for(j=0;j<N_p;j++){
				if(j!=0&&j!=N_p-1)
					Force[j]=-(Potential[j+1]-Potential[j-1])/(2*sg_get_r(thePoisson,i)*dphi);
				else if(j==0)
					Force[0]=-(Potential[1]-Potential[N_p-1])/(2*sg_get_r(thePoisson,i)*dphi);
				else 
					Force[N_p-1]=-(Potential[0]-Potential[N_p-2])/(2*sg_get_r(thePoisson,i)*dphi);
			}
			Force+=N_p;
			Potential+=N_p;
		}
	else if(DIRECTION==R_DIR)
	for(k=0;k<N_z;k++)
                for(i=0;i<N_r;i++){
			if(i!=0&&i!=N_r-1)
                        for(j=0;j<N_p;j++)
                                Force[j]=-(Potential[j+N_p]-Potential[j-N_p])/(2*dr);
                        else if(i==0)
			for(j=0;j<N_p;j++)
				Force[j]=-(Potential[j+N_p]-Potential[j])/dr;
			else 
			for(j=0;j<N_p;j++)
				Force[j]=-(Potential[j]-Potential[j-N_p])/dr;
                        Force+=N_p;
                        Potential+=N_p;
                }
	else if(DIRECTION==Z_DIR)
	        for(k=0;k<N_z;k++)
                for(i=0;i<N_r;i++){
                        if(k!=0&&k!=N_z-1)
                        for(j=0;j<N_p;j++)
                                Force[j]=-(Potential[j+N_p*N_r]-Potential[j-N_p*N_r])/(2*dz);
                        else if(k==0)
                        for(j=0;j<N_p;j++)
                                Force[j]=-(Potential[j+N_p*N_r]-Potential[j])/dz;
                        else
                        for(j=0;j<N_p;j++)
                                Force[j]=-(Potential[j]-Potential[j-N_p*N_r])/dz;
                        Force+=N_p;
                        Potential+=N_p;
                }
	else
		printf("wrong DIRECTION in cal_force()!\n");
}


	
