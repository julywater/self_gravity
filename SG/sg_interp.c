#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../SG.h"
#include"../paul.h"
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
void cylinder_interp(struct Sim *theSim,struct Cell ***theCells,struct poisson *thePoisson,struct MPIsetup *theMPIsetup){
	int i,j,k;
	int N_p=thePoisson->N_p;
	int N_r=thePoisson->N_r;
	int N_z=thePoisson->N_z;
	double *rho_new;
	double *rho_old=malloc(sizeof(double)*sim_N_p(theSim,thePoisson->rmax));
	double *phi_old=malloc(sizeof(double)*sim_N_p(theSim,thePoisson->rmax));
	int rank=mpisetup_MyProc(theMPIsetup);
	for( k=0;k<N_z; ++k ){
		int k_disco=k+thePoisson->zmin;
		for( i=0; i<N_r; ++i ){
			int i_disco=i+thePoisson->rmin;
			rho_new=thePoisson->density+(k*N_r+i)*N_p;
			for( j=0 ; j<sim_N_p(theSim,i_disco) ; ++j){
				struct Cell* c=cell_single(theCells,i_disco,j,k_disco);
		        	if(cell_prim(c,RHO)<pow(10,-4))
					rho_old[j]=0;
				else
					rho_old[j]=4.0*M_PI*cell_prim(c,RHO)*pow(thePoisson->dr,2);
				phi_old[j]=cell_tiph(c)-0.5*cell_dphi(c);
			}
		        interp(phi_old,rho_old,sim_N_p(theSim,i_disco),rho_new,N_p,FORWARD);
               } 
        }
	free(phi_old);
	free(rho_old);   

}	

void disco_interp(struct Sim *theSim,struct Cell ***theCells,struct poisson *thePoisson){
	int i,j,k;
        int N_p=thePoisson->N_p;
        int N_r=thePoisson->N_r;
        int N_z=thePoisson->N_z;
        double *V=thePoisson->density;
	int rmax=thePoisson->rmax;
        double *phi_disco=malloc(sizeof(double)*sim_N_p(theSim,rmax));
        double *V_disco=malloc(sizeof(double)*sim_N_p(theSim,rmax));
/*
	FILE *fp;
	fp=fopen("potential.dat","w");
	double dz=1.0/N_z;
        double dr=1.0/N_r;
	for( k=0;k<N_z; ++k ){
                double z=(k+0.5)*dz;
                for( i=0; i<N_r; ++i ){
                        double r=(i+0.5)*dr;
			V=thePoisson->density+(k*N_r+i)*N_p;
			fprintf(fp,"%f   %f   %f\n",z,r,*V);
	}
	}
	fclose(fp);
*/
	for( k=0;k<N_z; ++k ){
		int k_disco=k+thePoisson->zmin;
		double zm = sim_FacePos(theSim,k_disco-1,Z_DIR);
		double zp = sim_FacePos(theSim,k_disco,Z_DIR);
		double z=0.5*(zm+zp);
                for( i=0;i<N_r; ++i ){
			int i_disco=i+thePoisson->rmin;
		        double rm = sim_FacePos(theSim,i-1,R_DIR);
                        double rp = sim_FacePos(theSim,i,R_DIR);
			double r=0.5*(rm+rp);
                        V=thePoisson->density+(k*N_r+i)*N_p;
                        for( j=0 ; j<sim_N_p(theSim,i_disco) ; ++j){
                                struct Cell* c=cell_single(theCells,i_disco,j,k_disco);
                                phi_disco[j]=cell_tiph(c)-0.5*cell_dphi(c);
                        }

                        interp(phi_disco,V_disco,sim_N_p(theSim,i_disco),V,N_p,BACKWARD);
			for(j=0;j<sim_N_p(theSim,i_disco) ; ++j){
				struct Cell* c=cell_single(theCells,i_disco,j,k_disco);
				cell_write_V(c,V_disco[j]);
//				double phi=cell_tiph(c)-0.5*cell_dphi(c);
//				double x=r*cos(phi);
//				double y=r*sin(phi);
//				double R=sqrt(pow((x-0.0),2)+y*y+z*z);
//				if(z>0&&z<=fabs(thePoisson->dr))
//				fprintf(fp,"%f   %f\n",R,V_disco[j]);
			}
//			free(phi_disco);
//                        free(V_disco);
		}
		
        }
	free(phi_disco);
	free(V_disco);

}


void disco_force_interp(struct domain* theDomain,struct poisson *thePoisson,int direction){
    int i,j,k;
    int N_p=thePoisson->N_p;
    int N_r=thePoisson->N_r;
    int N_z=thePoisson->N_z;
    double *Force=thePoisson->buffer;
	int rmax=thePoisson->rmax;
	double *phi_disco=malloc(sizeof(double)*sim_N_p(theSim,rmax));
	double *force_disco=malloc(sizeof(double)*sim_N_p(theSim,rmax));
        for( k=0;k<N_z; ++k ){
                int k_disco=k+thePoisson->zmin;
                for( i=0;i<N_r; ++i ){
                        int i_disco=i+thePoisson->rmin;
                        Force=thePoisson->buffer+(k*N_r+i)*N_p;
                       
//get phi distribution of the disco grid
                        for( j=0 ; j<sim_N_p(theSim,i_disco) ; ++j){
                                struct Cell* c=cell_single(theCells,i_disco,j,k_disco);
                                phi_disco[j]=cell_tiph(c)-0.5*cell_dphi(c);
                        }

                        interp(phi_disco,force_disco,sim_N_p(theSim,i_disco),Force,N_p,BACKWARD);
                        for(j=0;j<sim_N_p(theSim,i_disco) ; ++j){
                                struct Cell* c=cell_single(theCells,i_disco,j,k_disco);
                                cell_write_Force(c,force_disco[j],direction);
						}


                }

        }
	free(phi_disco);
	free(force_disco);
}



