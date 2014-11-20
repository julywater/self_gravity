#include"paul.h"
#include"SG.h"
void planet_RK_copy( struct planet * );
void onestep( struct domain * , struct poisson *,double , double , int , double );
void add_diagnostics( struct domain * , double );

void timestep( struct domain * theDomain ,struct poisson * thepoisson, double dt ){
   
   struct cell ** theCells = theDomain->theCells;
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
   int * Np = theDomain->Np;
   int Npl = theDomain->Npl;

   int i,jk,p;


   for( jk=0 ; jk<Nr*Nz ; ++jk ){
      for( i=0 ; i<Np[jk] ; ++i ){
         struct cell * c = &(theCells[jk][i]);
         memcpy( c->RKcons , c->cons , NUM_Q*sizeof(double) );
         c->RKpiph = c->piph; 
      }
   }
   for( p=0 ; p<Npl ; ++p ){
      planet_RK_copy( theDomain->thePlanets + p );
   }

   onestep( theDomain ,thepoisson, 0.0 ,     dt , 0 , dt );
   onestep( theDomain ,thepoisson, 0.5 , 0.5*dt , 1 , dt );

   add_diagnostics( theDomain , dt );
   theDomain->t += dt;   
   theDomain->count_steps += 1;

int rank=theDomain->rank;
  if(rank==0)
  	printf("time=%f\n",theDomain->t);
}
