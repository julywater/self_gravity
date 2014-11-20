//#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include"../paul.h"
#include "../SG.h"
double get_dV( double * , double * );
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

void cal_multipole(struct domain* theDomain,struct poisson* thePoisson){
        int i,j,k;
        double mass_center[3]={0,0,0};
        double total_mass=0;
		double * r_jph = theDomain->r_jph;
		double * z_kph = theDomain->z_kph;
		int * Np = theDomain->Np;
        for(k=thePoisson->zmin;k<thePoisson->zmax;k++){
                double zm = z_kph[k-1];
                double zp = z_kph[k]
                double z=0.5*(zm+zp);
                for(j=thePoisson->rmin;j<thePoisson->rmax;j++){
                        double rm = r_jph[j-1];
                        double rp = r_jph[j];
                        double r=0.5*(rm+rp);
						int jk=j+theDomain->Nr*k;
                        for(i=0;i<Np[jk]; i++){
                            struct cell * c = &(theCells[jk][i]);
							double phip = c->piph;
							double phim = c->piph - c->dphi;
							double xp[3] = {r_jph[j]  ,phip,z_kph[k]  };
							double xm[3] = {r_jph[j-1],phim,z_kph[k-1]};
							double dV = get_dV(xp,xm);
                            double x=r*cos(phi);
                            double y=r*sin(phi);
                            double mass_element;
                                 if(cell_prim(c,RHO)>pow(10,-4))
                                        mass_element=cell_prim(c,RHO)*dV;
                                 else
                                        mass_element=0;
                                 mass_center[0]+=mass_element*x;
                                 mass_center[1]+=mass_element*y;
                                 mass_center[2]+=mass_element*z;
                                 total_mass+=mass_element;
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


void cal_r_bndry(int rbndry,int rinner,struct domain* theDomain,struct poisson* thePoisson){
        int i,j,k;
        int N_p=thePoisson->N_p;
        int N_r=thePoisson->N_r;
		int N_z=thePoisson->N_z;
		double r_size=thePoisson->r_size;
		double * rho=thePoisson->density;
		double * r_jph = theDomain->r_jph;
		double * z_kph = theDomain->z_kph;
        double rm = r_jph[rbndry-1];
        double rp = r_jph[rbndry];
        double r=0.5*(rm+rp);
        for(k=0;k<N_z;k++){
			int k_disco=k+thePoisson->zmin;
            double zm = z_kph[k_disco-1];
            double zp = z_kph[k_disco];
            double z=0.5*(zm+zp);
            for(j=0;j<N_p; ++j){
			    double phi=j*2*M_PI/N_p;
				double Vbndry=getpotential(thePoisson,r,z,phi);
				int index=(k*N_r+N_r-1)*N_p+j;
				rho[index]+=-(1.0+1.0/(2*thePoisson->N_r_glob+1.0))*Vbndry;
			}       		
        }
		
	
}
void cal_z_bndry(int zbndry,int zinner,struct domain* theDomain,struct poisson* thePoisson){
        int i,j,k;
		k=zinner-thePoisson->zmin;
        int N_p=thePoisson->N_p;
        int N_r=thePoisson->N_r;
        int N_z=thePoisson->N_z;
//		double dr=thePoisson->dr;
        double * rho=thePoisson->density;
		double * r_jph = theDomain->r_jph;
		double * z_kph = theDomain->z_kph;
        double zm = z_kph[zbndry-1];
        double zp = z_kph[zbndry];
        double z=0.5*(zm+zp);		
        for(i=0;i<N_r;i++){
			int i_disco=i+thePoisson->rmin;
            double rm = r_jph[i_disco-1];
            double rp = r_jph[i_disco];
            double r=0.5*(rm+rp);
            for(j=0;j<N_p; ++j){
                double phi=j*2*M_PI/N_p;
                double Vbndry=getpotential(thePoisson,r,z,phi);
				int index=(k*N_r+i)*N_p+j;
				rho[index]+=-pow((thePoisson->dr/thePoisson->dz),2)*Vbndry;
			}
        }
						
}
void set_bndry(struct domain* theDomain,struct poisson *thePoisson){
	cal_multipole(theDomain,thePoisson);
	double zmin=thePoisson->zmin;
	double zmax=thePoisson->zmax;
	double rmin=thePoisson->rmin;
	double rmax=thePoisson->rmax;
	if(mpisetup_check_rout_bndry(theDomain)){
		int rbndry=rmax;
		int rinner=rmax-1;
		cal_r_bndry(rbndry,rinner,thDomain,thePoisson);
	}
	 if(mpisetup_check_ztop_bndry(theDomain)){
		int zbndry=zmax;
		int zinner=zmax-1;	
		cal_z_bndry(zbndry,zinner,theDomain,thePoisson);
	}
	if(mpisetup_check_zbot_bndry(theDomain)){
		int zbndry=zmin-1;
		int zinner=zmin;
		cal_z_bndry(zbndry,zinner,theDomain,thePoisson);
	}
}


  			                     

