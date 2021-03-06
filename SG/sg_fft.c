#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<fftw3.h>
#include "../SG.h"

void density_fft(struct poisson* thePoisson){
	int i,j,k;
	int N_p=thePoisson->N_p;
	int N_r=thePoisson->N_r;
	int N_z=thePoisson->N_z;
	double *rho_p=thePoisson->density;
	const int n[1]={N_p};
	const fftw_r2r_kind kind[1]={FFTW_R2HC};
	fftw_plan p=fftw_plan_many_r2r(1,n,N_z*N_r,rho_p,n,1,N_p,rho_p,n,1,N_p,kind,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	fftw_cleanup();
	
}
void inverse_fft(struct poisson* thePoisson){
	int i,j,k;
        int N_p=thePoisson->N_p;
        int N_r=thePoisson->N_r;
        int N_z=thePoisson->N_z;
        double *V_p=thePoisson->density;
/*
	fftw_plan p;
        for(k=0;k<N_z;k++)
                for(i=0;i<N_r;i++){
			p=fftw_plan_r2r_1d(N_p,V_p,V_p,FFTW_HC2R, FFTW_ESTIMATE);
                        fftw_execute(p);
			int is=0;
			for(is=0;is<N_p;is++)
				V_p[is]=V_p[is]/N_p;
                        V_p+=N_p;
			fftw_destroy_plan(p);		
        }
*/
	const int n[1]={N_p};
        const fftw_r2r_kind kind[1]={FFTW_HC2R};
        fftw_plan p=fftw_plan_many_r2r(1,n,N_z*N_r,V_p,n,1,N_p,V_p,n,1,N_p,kind,FFTW_ESTIMATE);
        fftw_execute(p);
        fftw_destroy_plan(p);
	fftw_cleanup();
	int idx=0;
	for(idx=0;idx<N_p*N_r*N_z;idx++)
		V_p[idx]/=N_p;
}
//density array shape phi,r,z;
