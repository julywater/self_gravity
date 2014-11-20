
#include "paul.h"
#include "helm.h"
int mpiSetup( struct domain * , int , char *[] );
void setupGrid( struct domain * );
void timestep( struct domain * , double );
void setupCells( struct domain * );
void regrid( struct domain * );
void exchangeData( struct domain * , int );
double getmindt( struct domain * );

void read_par_file( struct domain * );

void setupDomain( struct domain * );
void freeDomain( struct domain * );
void check_dt( struct domain * , double * );
void possiblyOutput( struct domain * , int );

void start_clock( struct domain * );
void generate_log( struct domain * );

double * table;




int main( int argc , char * argv[] ){
 
   MPI_Init(&argc,&argv);
   struct domain theDomain = {0};
   start_clock( &theDomain ); 
   read_par_file( &theDomain );
   
   int error = mpiSetup(&theDomain,argc,argv);
   if( error==1 ) return(0);

   if(theDomain.rank==0) remove("abort");

//
   table=(double*)malloc(sizeof(double)*577013);
   read_helm_table_0_(table);
// read HELMHOLTZ EOS table

   setupGrid( &theDomain );   
   setupDomain( &theDomain );
/*
   setICparams( &theDomain );
   setHydroParams( &theDomain );
   setGeometryParams( &theDomain );
   setRiemannParams( &theDomain );
*/
   setupCells( &theDomain );
/*
   if( theDomain.theParList.Initial_Regrid && !(theDomain.theParList.restart_flag) ) regrid( &theDomain );
*/
   if( theDomain.Nr > 1 ) exchangeData( &theDomain , 0 );
   if( theDomain.Nz > 1 ) exchangeData( &theDomain , 1 );

   if( theDomain.rank==0 && !(theDomain.theParList.restart_flag) ){
      FILE * rFile = fopen("report.dat","w");
      fclose(rFile);
   }

   MPI_Barrier(theDomain.theComm);
   while( !(theDomain.final_step) ){

      double dt = getmindt( &theDomain );
      check_dt( &theDomain , &dt );
      possiblyOutput( &theDomain , 0 );
      timestep( &theDomain , dt );

   }

   possiblyOutput( &theDomain , 1 );
   generate_log( &theDomain );
   MPI_Barrier(theDomain.theComm);
   freeDomain( &theDomain );
   free(table);
   MPI_Finalize();

   return(0);

}

