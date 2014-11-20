#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../SG.h"
#include"../paul.h"
void cell_write_Force(struct cell*,double ,int);
void interp(double *x_disc,double *y_disc,int N,double *y_cylinder,int M,int direction){
	int i=1;
	if(x_disc[0]<0) x_disc[0]+=2*M_PI;
	for(i=1;i<N;i++){
		if(x_disc[i]<0) 
			x_disc[i]+=2*M_PI;
		if(x_disc[i]<x_disc[i-1])
			break;
	}
	int l,lm;
	if(i==N){
		i=0;
		l=i;
		lm=N-1;	
	}
	else{
		l=i;
		lm=l-1;
	}
	double dx=2*M_PI/M;
	int k=0;
	double *x_cylinder=(double *)malloc(sizeof(double)*(M+1));
	for(k=0;k<M+1;k++)
		x_cylinder[k]=k*dx;
//	
	if(direction==FORWARD){
	for(k=0;k<M;k++){
		if(x_cylinder[k]>=x_disc[lm]){
			y_cylinder[k]=(x_cylinder[k]-x_disc[lm])*(y_disc[l]-y_disc[lm])/(x_disc[l]+2*M_PI-x_disc[lm])+y_disc[lm];
			continue;
		}
			
		while(x_cylinder[k]>x_disc[i]){
			i++;
			if(i>=N)
				i-=N;
		}
		int im=i-1;
		if(im==-1)
			im=N-1;
		
		if(i!=l)
			y_cylinder[k]=(x_cylinder[k]-x_disc[im])*(y_disc[i]-y_disc[im])/(x_disc[i]-x_disc[im])+y_disc[im];
		else
			y_cylinder[k]=(x_cylinder[k]-(x_disc[im]-2*M_PI))*(y_disc[i]-y_disc[im])/(x_disc[i]-(x_disc[im]-2*M_PI))+y_disc[im];
	}
	}
	else if(direction==BACKWARD){
	int m=0;
	i=l;
	for(k=0;k<N;k++){
		while(x_disc[i]>=x_cylinder[m]){
			m++;
		}
	if(m!=M)
		y_disc[i]=(x_disc[i]-x_cylinder[m-1])*(y_cylinder[m]-y_cylinder[m-1])/(x_cylinder[m]-x_cylinder[m-1])+y_cylinder[m-1];
	else
		y_disc[i]=(x_disc[i]-x_cylinder[M-1])*(y_cylinder[0]-y_cylinder[M-1])/(2*M_PI-x_cylinder[M-1])+y_cylinder[M-1];
	
	i++;
	if(i>=N)
		i-=N;
	}
	}
	else
		printf("wrong parameter\n");
	free(x_cylinder);
}
void cylinder_interp(struct domain* theDomain,struct poisson *thePoisson){
	struct cell ** theCells = theDomain->theCells;
	int i,j,k;
	int N_p=thePoisson->N_p;
	int N_r=thePoisson->N_r;
	int N_z=thePoisson->N_z;
	double *rho_new;
	int rmax=thePoisson->rmax;
	int * Np_disco = theDomain->Np;
	double *rho_disco=malloc(sizeof(double)*Np_disco[rmax]);
	double *phi_disco=malloc(sizeof(double)*Np_disco[rmax]);

	for( k=0;k<N_z; ++k ){
		int k_disco=k+thePoisson->zmin;
		for( i=0; i<N_r; ++i ){
			int i_disco=i+thePoisson->rmin;
			int ik = i_disco+theDomain->Nr*k_disco;
			rho_new=thePoisson->density+(k*N_r+i)*N_p;
			for( j=0 ; j<Np_disco[ik]; ++j){
				struct cell * c = &(theCells[ik][j]);
				if(Np_disco[ik]>NP_disco[rmax]){
					realloc(phi_disco,Np_disco[ik]);
					realloc(rho_disco,Np_disco[ik]);
				}
		        if(cell_prim(c,RHO)<pow(10,-4))
					rho_disco[j]=0;
				else
					rho_disco[j]=4.0*M_PI*cell_prim(c,RHO)*pow(thePoisson->dr,2);
				phi_disco[j]=cell_tiph(c)-0.5*cell_dphi(c);
			}
		    interp(phi_disco,rho_disco,Np_disco[ik],rho_new,N_p,FORWARD);
        } 
    }
	free(phi_disco);
	free(rho_disco);   

}	
/*
void disco_interp(struct domain *theDomain,struct poisson *thePoisson){
	struct cell ** theCells = theDomain->theCells;
	int i,j,k;
    int N_p=thePoisson->N_p;
    int N_r=thePoisson->N_r;
    int N_z=thePoisson->N_z;
    double *V=thePoisson->density;
	int rmax=thePoisson->rmax;
    double *phi_disco=malloc(sizeof(double)*sim_N_p(theSim,rmax));
    double *V_disco=malloc(sizeof(double)*sim_N_p(theSim,rmax));
    int * Np_disco = theDomain->Np;
	
	for( k=0;k<N_z; ++k ){
		int k_disco=k+thePoisson->zmin;
                for( i=0;i<N_r; ++i ){
					int i_disco=i+thePoisson->rmin;
					int ik = i_disco+theDomain->Nr*k_disco;
                    V=thePoisson->density+(k*N_r+i)*N_p;
                    for( j=0;j<Np_disco[ik]; ++j){
                        struct Cell* c = &(theCells[ik][j]);
                        phi_disco[j]=cell_tiph(c)-0.5*cell_dphi(c);
                    }

                        interp(phi_disco,V_disco,Np_disco[ik],V,N_p,BACKWARD);
					for(j=0;j<Np_disco[ik]; ++j){
						struct Cell* c = &(theCells[ik][j]);
						cell_write_V(c,V_disco[j]);
					}
		}
    }
	free(phi_disco);
	free(V_disco);

}
*/

void disco_force_interp(struct domain* theDomain,struct poisson *thePoisson,int direction){
	struct cell ** theCells = theDomain->theCells;
    int i,j,k;
    int N_p=thePoisson->N_p;
    int N_r=thePoisson->N_r;
    int N_z=thePoisson->N_z;
    double *Force=thePoisson->buffer;
	int rmax=thePoisson->rmax;
	int * Np_disco = theDomain->Np;
	double *phi_disco=malloc(sizeof(double)*Np_disco[rmax]);
	double *force_disco=malloc(sizeof(double)*Np_disco[rmax]);
        for( k=0;k<N_z; ++k ){
                int k_disco=k+thePoisson->zmin;
                for( i=0;i<N_r; ++i ){
                        int i_disco=i+thePoisson->rmin;
						int ik = i_disco+theDomain->Nr*k_disco;
                        Force=thePoisson->buffer+(k*N_r+i)*N_p;
						if(Np_disco[ik]>NP_disco[rmax]){
							realloc(phi_disco,Np_disco[ik]);
							realloc(force_disco,Np_disco[ik]);
						}
//get phi distribution of the disco grid
                        for( j=0 ; j<Np_disco[ik] ; ++j){
                            struct cell * c = &(theCells[ik][j]);
                            phi_disco[j]=cell_tiph(c)-0.5*cell_dphi(c);
                        }
                        interp(phi_disco,force_disco,sim_N_p(theSim,i_disco),Force,N_p,BACKWARD);
                        for(j=0;j<Np_disco[ik]; ++j){
                            struct cell * c = &(theCells[ik][j]);
                            cell_write_Force(c,force_disco[j],direction);
						}
                }
        }
	free(phi_disco);
	free(force_disco);
}



