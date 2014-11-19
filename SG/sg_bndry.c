#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>

#include "../Headers/SG.h"
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/header.h"
/*
void cal_multipole(struct Sim* theSim,struct Cell*** theCells,struct poisson* thePoisson){
	int i,j,k;
	double mass_center[2][3]={{0,0,0},{0,0,0}};
	double total_mass[2]={0,0};
        const int N_Q=sim_NUM_Q(theSim);
	for(k=thePoisson->zmin;k<thePoisson->zmax;++k){
		double zm = sim_FacePos(theSim,k-1,Z_DIR);
                double zp = sim_FacePos(theSim,k,Z_DIR);
        	double z=0.5*(zm+zp);
		double dz=zp-zm;
		for(i=thePoisson->rmin;i<thePoisson->rmax;i++){
	                double rm = sim_FacePos(theSim,i-1,R_DIR);
	                double rp = sim_FacePos(theSim,i,R_DIR);
        	        double r=0.5*(rm+rp);
			for(j=0;j<sim_N_p(theSim,i) ; ++j ){
                		 struct Cell * c = cell_single(theCells,i,j,k);
		                 double phi= cell_tiph(c)-0.5*cell_dphi(c);
                		 double dphi = cell_dphi(c);
                		 double dV = dphi*0.5*(rp*rp-rm*rm)*dz;
                                 double x=r*cos(phi);
                                 double y=r*sin(phi);
				 double mass_element;
				 if(cell_prim(c,RHO)>pow(10,-3))
					mass_element=cell_prim(c,RHO)*dV;
				 else
					mass_element=0;
				 double q[2];
				 int is=0;
				 for(is=0;is<N_Q;is++){
                                         double M=mass_element*cell_prim(c,5+is);
                                	 mass_center[is][0]+=M*x;
                              		 mass_center[is][1]+=M*y;
                                 	 mass_center[is][2]+=M*z;
			         	 total_mass[is]+=M;
				}
			}
		}
	}
	MPI_Allreduce(total_mass,(thePoisson->total_mass),2,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(mass_center[0],thePoisson->mass_center[0],3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(mass_center[1],thePoisson->mass_center[1],3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	int idx;
	for(idx=0;idx<3;idx++){
		thePoisson->mass_center[0][idx]/=thePoisson->total_mass[0];
		thePoisson->mass_center[1][idx]/=thePoisson->total_mass[1];
	}
//	printf("%f  %f  %f  %f\n",thePoisson->total_mass,thePoisson->mass_center[0],thePoisson->mass_center[1],thePoisson->mass_center[2]);
}
	
	
                        

double getpotential(struct poisson* thePoisson,double r,double z,double phi){
	double R[2]={0};
	double x=r*cos(phi);
	double y=r*sin(phi);
	int i;
	double V=0;
	for(i=0;i<2;i++){
		R[i]=sqrt(pow((x-thePoisson->mass_center[i][0]),2)+pow((y-thePoisson->mass_center[i][1]),2)+pow((z-thePoisson->mass_center[i][2]),2));
		V+=-thePoisson->total_mass[i]/R[i];
	}	
	return V;
}
*/

void cal_multipole(struct Sim* theSim,struct Cell*** theCells,struct poisson* thePoisson){
        int i,j,k;
        double mass_center[3]={0,0,0};
        double total_mass=0;
        for(k=thePoisson->zmin;k<thePoisson->zmax;++k){
                double zm = sim_FacePos(theSim,k-1,Z_DIR);
                double zp = sim_FacePos(theSim,k,Z_DIR);
                double z=0.5*(zm+zp);
                double dz=zp-zm;
                for(i=thePoisson->rmin;i<thePoisson->rmax;i++){
                        double rm = sim_FacePos(theSim,i-1,R_DIR);
                        double rp = sim_FacePos(theSim,i,R_DIR);
                        double r=0.5*(rm+rp);
                        for(j=0;j<sim_N_p(theSim,i) ; ++j ){
                                 struct Cell * c = cell_single(theCells,i,j,k);
                                 double phi= cell_tiph(c)-0.5*cell_dphi(c);
                                 double dphi = cell_dphi(c);
                                 double dV = dphi*0.5*(rp*rp-rm*rm)*dz;
                                 double x=r*cos(phi);
                                 double y=r*sin(phi);
                                 double mass_element;
                                 if(cell_prim(c,RHO)>pow(10,-3))
                                        mass_element=cell_prim(c,RHO)*dV;
                                 else
                                        mass_element=0;
				 const double M=mass_element;

                                 mass_center[0]+=M*x;
                                 mass_center[1]+=M*y;
                                 mass_center[2]+=M*z;
                                 total_mass+=M;
                                }
                        }
                }
        MPI_Allreduce(&total_mass,&(thePoisson->total_mass),1,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	MPI_Allreduce(mass_center,thePoisson->mass_center,3,MPI_DOUBLE,MPI_SUM,MPI_COMM_WORLD);
	thePoisson->mass_center[0]/=thePoisson->total_mass;
	thePoisson->mass_center[1]/=thePoisson->total_mass;
	thePoisson->mass_center[2]/=thePoisson->total_mass;
	}
double getpotential(struct poisson* thePoisson,double r,double z,double phi){
        double x=r*cos(phi);
        double y=r*sin(phi);
        double R=sqrt(pow((x-thePoisson->mass_center[0]),2)+pow((y-thePoisson->mass_center[1]),2)+pow((z-thePoisson->mass_center[2]),2));
        double V=-thePoisson->total_mass/R; 
//	double V=-2.50662827463/sqrt(x*x+y*y+z*z);
        return V;
}


void cal_r_bndry(int rbndry,int rinner,struct Sim* theSim,struct poisson* thePoisson){
                int i,j,k;
                int N_p=thePoisson->N_p;
                int N_r=thePoisson->N_r;
		int N_z=thePoisson->N_z;
		double r_size=thePoisson->r_size;
		double * rho=thePoisson->density;
                double rm = sim_FacePos(theSim,rbndry-1,R_DIR);
                double rp = sim_FacePos(theSim,rbndry,R_DIR);
                double r=0.5*(rm+rp);
                for(k=0;k<N_z;k++){
			int k_disco=k+thePoisson->zmin;
                        double zm = sim_FacePos(theSim,k_disco-1,Z_DIR);
                        double zp = sim_FacePos(theSim,k_disco,Z_DIR);
                        double z=0.5*(zm+zp);
                        for(j=0;j<N_p; ++j){
                                double phi=j*2*M_PI/N_p;
                                double Vbndry=getpotential(thePoisson,r,z,phi);
				int index=(k*N_r+N_r-1)*N_p+j;
	         		rho[index]+=-(1.0+1.0/(2*thePoisson->N_r_glob+1.0))*Vbndry;
			}
                       		
                }
		
	
}
void cal_z_bndry(int zbndry,int zinner,struct Sim* theSim,struct poisson* thePoisson){
                int i,j,k;
		k=zinner-thePoisson->zmin;
                int N_p=thePoisson->N_p;
                int N_r=thePoisson->N_r;
                int N_z=thePoisson->N_z;
//		double dr=thePoisson->dr;
                double * rho=thePoisson->density;
                double zm = sim_FacePos(theSim,zbndry-1,Z_DIR);
                double zp = sim_FacePos(theSim,zbndry,Z_DIR);
                double z=0.5*(zm+zp);		
                for(i=0;i<N_r;i++){
			int i_disco=i+thePoisson->rmin;
                        double rm = sim_FacePos(theSim,i_disco-1,R_DIR);
                        double rp = sim_FacePos(theSim,i_disco,R_DIR);
                        double r=0.5*(rm+rp);
                        for(j=0;j<N_p; ++j){
                                double phi=j*2*M_PI/N_p;
                                double Vbndry=getpotential(thePoisson,r,z,phi);
				int index=(k*N_r+i)*N_p+j;
				rho[index]+=-pow((thePoisson->dr/thePoisson->dz),2)*Vbndry;
			}
                }
						
}
void set_bndry(struct Sim* theSim,struct Cell*** theCells,struct MPIsetup* theMPIsetup,struct poisson *thePoisson){
	cal_multipole(theSim,theCells,thePoisson);
	double zmin=thePoisson->zmin;
	double zmax=thePoisson->zmax;
	double rmin=thePoisson->rmin;
	double rmax=thePoisson->rmax;
	
	if(mpisetup_check_rout_bndry(theMPIsetup)){
		int rbndry=rmax;
		int rinner=rmax-1;
		cal_r_bndry(rbndry,rinner,theSim,thePoisson);
	}
	 if(mpisetup_check_ztop_bndry(theMPIsetup)){
		int zbndry=zmax;
		int zinner=zmax-1;	
		cal_z_bndry(zbndry,zinner,theSim,thePoisson);
	}
	if(mpisetup_check_zbot_bndry(theMPIsetup)){
		int zbndry=zmin-1;
		int zinner=zmin;
		cal_z_bndry(zbndry,zinner,theSim,thePoisson);
	}
}


  			                     

