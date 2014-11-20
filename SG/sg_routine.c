//#define SG_PRIVATE_DEFS
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include<fftw3.h>
#include "../SG.h"
#include "../paul.h"
/*
struct Tridiagsys{
	double *eigenvalues;
//	double *lamba;
	double *a;
	double *b;
	double *c;
	};
*/
void solvetridiag(double *a,double *b,double *c,double *d,double *v,int N){
	int i;	
	b[0]=b[0]/a[0];
	d[0]=d[0]/a[0];
	for(i=1;i<N;i++){
		double m=1.0/(a[i]-c[i]*b[i-1]);
		b[i]=b[i]*m;
		d[i]=(d[i]-d[i-1]*c[i])*m;}
	v[N-1]=d[N-1];	
	for(i=N-2;i>=0;i--)
		v[i]=d[i]-b[i]*v[i+1];	
	}
void solveVp(int rank,struct poisson* thePoisson){
	int i,j,k;
	int N_r_glob=thePoisson->N_r_glob;
	int N_k=thePoisson->N_k;
	int N_z_glob=thePoisson->N_z_glob;
	double * D_buf=thePoisson->shortbuffer;
	double * data=thePoisson->density;
	double n;
	double *V_buf=(double *)malloc(sizeof(double)*N_r_glob);
	for(j=0;j<N_k;j++){
		n=thePoisson->freq[thePoisson->k_start+j];	
		for(k=0;k<N_z_glob;k++){
			for(i=0;i<N_r_glob;i++){
				double lamba=1.0/(2*(i+1)-1);
				thePoisson->a[i]=thePoisson->eigenvalues[k]-4.0*pow(n,2)*pow(lamba,2);
				thePoisson->b[i]=1.0+1.0/(2*(i+1)-1);
				thePoisson->c[i]=1.0-1.0/(2*(i+1)-1);
				D_buf[i]=data[in_lookup(thePoisson,k,i,j)];
			}
			solvetridiag(thePoisson->a,thePoisson->b,thePoisson->c,D_buf,V_buf,N_r_glob);
			for(i=0;i<N_r_glob;i++)
				data[in_lookup(thePoisson,k,i,j)]=V_buf[i];
		}
	}	
	free(V_buf);
}
int in_lookup(struct poisson *thePoisson,int k,int i,int j){
	if(i>=thePoisson->N_r_glob||k>=thePoisson->N_z_glob||j>=thePoisson->N_k){
		printf("in_lookup  index out of bound\n");
		return -1;
	}
	double *data=thePoisson->density;
//	int p=i/thePoisson->N_r;
	int q=k/thePoisson->N_z;
//	int ip=i-p*thePoisson->N_r;
	int kp=k-q*thePoisson->N_z;
	int p=0;
	
	while(thePoisson->N0_r[p+1]<=i)
			p++;
		
	int ip=i-thePoisson->N0_r[p];
	int rank=p*thePoisson->z_size+q;
	int index=array_sum(thePoisson->recv_elements,rank)+(kp*thePoisson->N_r_array[rank]+ip)*thePoisson->N_k+j;
	return index;
}

void sinefft(struct poisson* thePoisson){
	int i,j,k;
	double *in=thePoisson->density;
	const int N_element=thePoisson->N_z_glob;
	int N_k=thePoisson->N_k;
	int N_r_glob=thePoisson->N_r_glob;
	int N_z_glob=thePoisson->N_z_glob;
	const double uni=sqrt(0.5*(N_element+1));
/*
	double * sine_buf=thePoisson->shortbuffer;
	fftw_plan p;
	for(j=0;j<N_k;j++)
	for(i=0;i<N_r_glob;i++){
		for(k=0;k<N_z_glob;k++)
		sine_buf[k]=in[in_lookup(thePoisson,k,i,j)];
		p=fftw_plan_r2r_1d(N_element, sine_buf,sine_buf,FFTW_RODFT00,FFTW_ESTIMATE);
		fftw_execute(p);
		for(k=0;k<N_z_glob;k++)
			in[in_lookup(thePoisson,k,i,j)]=sine_buf[k]/2.0/uni;
		fftw_destroy_plan(p);
	}
*/
	
	double * sine_buf=thePoisson->buffer;
	int idx=0;
	for(j=0;j<N_k;j++)
		for(i=0;i<N_r_glob;i++)
			for(k=0;k<N_z_glob;k++)
				sine_buf[idx++]=in[in_lookup(thePoisson,k,i,j)];
	const int n[1]={N_element};
        const fftw_r2r_kind kind[1]={FFTW_RODFT00};
        fftw_plan p=fftw_plan_many_r2r(1,n,N_k*N_r_glob,sine_buf,n,1,N_element,sine_buf,n,1,N_element,kind,FFTW_ESTIMATE);
	fftw_execute(p);
	fftw_destroy_plan(p);
	idx=0;
	i=0;j=0;k=0;
	for(j=0;j<N_k;j++)
                for(i=0;i<N_z_glob;i++)
                        for(k=0;k<N_z_glob;k++)
                                in[in_lookup(thePoisson,k,i,j)]=sine_buf[idx++]/2.0/uni;
	fftw_cleanup();	
//nnote: in and out should have the same size;same size after sine transform
	}
void mpi_arrange(double *in,double * out,struct poisson *thePoisson){
	int i,j,k;
	for(j=0;j<thePoisson->size;j++)
	for(k=0;k<thePoisson->N_z;k++)
		for(i=0;i<thePoisson->N_r;i++){
			int in_index=(k*thePoisson->N_r+i)*thePoisson->N_p+j*thePoisson->N_k;
			int out_index=j*(thePoisson->N_r*thePoisson->N_z*thePoisson->N_k)+(k*thePoisson->N_r+i)*thePoisson->N_k;
			//or out+=N_k
			int is=0;
			for(is=0;is<thePoisson->N_k;is++)
				out[out_index+is]=in[in_index+is];
	}
}
void inverse_mpi_arrange(double *in,double* out,struct poisson *thePoisson){
	int i,j,k;
	for(j=0;j<thePoisson->size;j++)
	for(k=0;k<thePoisson->N_z;k++)
		for(i=0;i<thePoisson->N_r;i++){
			int out_index=(k*thePoisson->N_r+i)*thePoisson->N_p+j*thePoisson->N_k;
			int in_index=j*(thePoisson->N_r*thePoisson->N_z*thePoisson->N_k)+(k*thePoisson->N_r+i)*thePoisson->N_k;
			//or out+=N_k
			int is=0;
			for(is=0;is<thePoisson->N_k;is++)
				out[out_index+is]=in[in_index+is];
	}
}
void tridiag_setup(struct poisson* thePoisson){
	int Nx=thePoisson->N_r_glob;
        int Ny=thePoisson->N_z_glob;
	int i;
        for(i=0;i<Ny;i++)
                thePoisson->eigenvalues[i]=-4.0+2*cos((i+1)*M_PI/(Ny+1));
}
/*	
void destroy_tridiagsys(struct Tridiagsys* tri){
	free(tri->a);
	free(tri->b);
	free(tri->c);
	free(tri->eigenvalues);
	free(tri);
}
*/
void sg_route(struct Domain *theDomain,struct poisson *thePoisson){
//	struct Tridiagsys* thetridiag=alloc_tridiagsys(thePoisson);
	tridiag_setup(thePoisson);
	int N_r=thePoisson->N_r;
	int N_z=thePoisson->N_z;
	int N_p=thePoisson->N_p;
	int N_k=thePoisson->N_k;
	int i,j,k;
	int size=thePoisson->size;
	int rank=theDomain->rank;
	cylinder_interp(theSim,theCells,thePoisson,theMPIsetup);
	set_bndry(theSim,theCells,theMPIsetup,thePoisson);
	density_fft(thePoisson);
	mpi_arrange(thePoisson->density,thePoisson->buffer,thePoisson);
	double *buffersend=thePoisson->buffer;
	double *bufferstore=thePoisson->density;

	int *sendcnts=thePoisson->sendcnts;
    int *sdispls=thePoisson->sdispls;
    int *recvcnts=thePoisson->recvcnts;
    int *rdispls=thePoisson->rdispls;
	
    MPI_Alltoallv(buffersend,sendcnts,sdispls,MPI_DOUBLE,bufferstore,recvcnts,rdispls,MPI_DOUBLE,MPI_COMM_WORLD);
    sinefft(thePoisson);
    solveVp(rank,thePoisson);
    sinefft(thePoisson);    
/*
	if(rank==0){
	int i,j,k;
	i=0;j=0;k=0;
	FILE *f;
	f=fopen("poten.dat","w");
	for(i=0;i<thePoisson->N_r_glob;i++)
		for(k=0;k<thePoisson->N_z_glob;k++)
			fprintf(f,"%f   %d   %f\n",i+0.5,k,thePoisson->density[in_lookup(thePoisson,k,i,0)]);
	fclose(f);
	}
*/

	buffersend=thePoisson->density;
	bufferstore=thePoisson->buffer;
	//this is the inverse MPI comm.SO the sendcnts and recvcnts are exchanged 
	MPI_Alltoallv(buffersend,recvcnts,rdispls,MPI_DOUBLE,bufferstore,sendcnts,sdispls,MPI_DOUBLE,MPI_COMM_WORLD);
	inverse_mpi_arrange(thePoisson->buffer,thePoisson->density,thePoisson);
	inverse_fft(thePoisson);
	int direction=0;
	for(direction=0;direction<3;direction++){
		cal_force(thePoisson,direction);	
		disco_force_interp(theDomain,thePoisson,direction);
	}
//	disco_interp(theSim,theCells,thePoisson);
//	destroy_tridiagsys(thetridiag);

}
