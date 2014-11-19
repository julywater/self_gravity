#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<fftw3.h>
#include "../Headers/SG.h"
#include "../Headers/Cell.h"
#include "../Headers/Sim.h"
#include "../Headers/Face.h"
#include "../Headers/GravMass.h"
#include "../Headers/MPIsetup.h"
#include "../Headers/header.h"
int array_sum(int *a,int i){
	if(i==0)
		return 0;
	else{	
		int tot=0;
		int k;
		for(k=0;k<i;k++)
			tot+=a[k];
		return tot;
	}
}
struct poisson* poisson_create(struct Sim *theSim,struct MPIsetup *theMPIsetup){
	struct poisson* thePoisson=(struct poisson*)malloc(sizeof(struct poisson));
	int N_p=sim_N_global(theSim,R_DIR);
	thePoisson->N_p=N_p;
        thePoisson->zmin=sim_Nghost_min(theSim,Z_DIR);
	thePoisson->zmax=sim_N(theSim,Z_DIR)-sim_Nghost_max(theSim,Z_DIR);
	thePoisson->N_z=thePoisson->zmax-thePoisson->zmin;
	if(mpisetup_check_rin_bndry(theMPIsetup)){
		thePoisson->rmin=0;
		thePoisson->rmax=sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
	}
	else{
		thePoisson->rmin=sim_Nghost_min(theSim,R_DIR);
		thePoisson->rmax=sim_N(theSim,R_DIR)-sim_Nghost_max(theSim,R_DIR);
	}
	thePoisson->N_r=thePoisson->rmax-thePoisson->rmin;

//	printf("%d\n",thePoisson->N_r);

	//make sure N_r and  N_z are the same for all the procs.
	
        thePoisson->z_size=mpisetup_dim_NumProcs(theMPIsetup,Z_DIR);	
	thePoisson->r_size=mpisetup_dim_NumProcs(theMPIsetup,R_DIR);
	thePoisson->z_rank=mpisetup_dim_MyProc(theMPIsetup,Z_DIR);
	thePoisson->r_rank=mpisetup_dim_MyProc(theMPIsetup,R_DIR);
	thePoisson->N_r_glob=sim_N_global(theSim,R_DIR);	
	thePoisson->N_z_glob=sim_N_global(theSim,Z_DIR);


	int size=mpisetup_NumProcs(theMPIsetup);
	thePoisson->size=size;
	thePoisson->dr=sim_FacePos(theSim,1,R_DIR)-sim_FacePos(theSim,0,R_DIR);
        thePoisson->dz=sim_FacePos(theSim,1,Z_DIR)-sim_FacePos(theSim,0,Z_DIR);
 	//make sure dr is the same for all the procs in the first step	
	int N_k=thePoisson->N_p/thePoisson->size;//N_p can be divided by size (at current stage)
	thePoisson->N_k=N_k;
	int rank=mpisetup_MyProc(theMPIsetup);
	thePoisson->k_start=rank*N_k;
//	thePoisson->densityold=(double *)malloc(sizeof(double)*thePoisson->N_p);
	int block_size=thePoisson->N_p*thePoisson->N_r*thePoisson->N_z;
	int slice_size=thePoisson->N_k*thePoisson->N_r_glob*thePoisson->N_z_glob;
	int larger_one;
	if(block_size>=slice_size)
		larger_one=block_size;
	else
		larger_one=slice_size;
	thePoisson->density=(double *)malloc(sizeof(double)*larger_one);
	thePoisson->buffer=(double *)malloc(sizeof(double)*larger_one);
	if(thePoisson->N_r_glob>=thePoisson->N_z_glob)
		thePoisson->shortbuffer=(double *)malloc(sizeof(double)*thePoisson->N_r_glob);
	else
		thePoisson->shortbuffer=(double *)malloc(sizeof(double)*thePoisson->N_z_glob);

	thePoisson->freq=(double *)malloc(sizeof(double)*N_p);
	
	int i=0;
	for(i=0;i<N_p;i++){
		if(i<=N_p/2)		
			thePoisson->freq[i]=i;	
		else
			thePoisson->freq[i]=N_p-i;
	}

//MPI all to all parameters
	thePoisson->sendcnts=(int *)malloc(sizeof(int)*size);
        thePoisson->sdispls=(int *)malloc(sizeof(int)*size);
        thePoisson->recvcnts=(int *)malloc(sizeof(int)*size);
        thePoisson->rdispls=(int *)malloc(sizeof(int)*size);
        int send_elements=thePoisson->N_r*thePoisson->N_z*thePoisson->N_k;
	thePoisson->recv_elements=malloc(sizeof(int)*size);
	MPI_Allgather(&send_elements,1,MPI_INT,thePoisson->recv_elements,1,MPI_INT,MPI_COMM_WORLD);
	
	thePoisson->N_r_array=malloc(sizeof(int)*size);
	MPI_Allgather(&thePoisson->N_r,1,MPI_INT,thePoisson->N_r_array,1,MPI_INT,MPI_COMM_WORLD);
        for(i=0;i<size;i++){
                thePoisson->sendcnts[i]=send_elements;
                thePoisson->sdispls[i]=i*send_elements;
                thePoisson->recvcnts[i]=thePoisson->recv_elements[i];
                thePoisson->rdispls[i]=array_sum(thePoisson->recv_elements,i);
        }	
	thePoisson->N0_r=malloc(sizeof(int)*(thePoisson->r_size+1));
        for(i=0;i<thePoisson->r_size+1;i++)
        	thePoisson->N0_r[i]=sim_N0_r(theSim,i);

	
//tridiag linear system parameters
	int Nx=thePoisson->N_r_glob;
        int Ny=thePoisson->N_z_glob;
        thePoisson->eigenvalues=(double *)malloc(sizeof(double)*Ny);
        thePoisson->a=(double *)malloc(sizeof(double)*Nx);
        thePoisson->b=(double *)malloc(sizeof(double)*Nx);
        thePoisson->c=(double *)malloc(sizeof(double)*Nx);
	return thePoisson;
    
}  
