#include "paul.h"
#include"SG.h"
void AMR( struct domain * ); 
void move_BCs( struct domain * , double );

void clean_pi( struct domain * );
void set_wcell( struct domain * );

void adjust_RK_cons( struct domain * , double );
void adjust_RK_planets( struct domain * , double );
void move_cells( struct domain * , double , double );
void calc_dp( struct domain * );
void calc_prim( struct domain * );

void setup_faces( struct domain * , struct face ** , int * , int );
void phi_flux( struct domain * , double dt );
void trans_flux( struct domain * , struct face * , int , double dt , int );
void add_source( struct domain * , double dt );

void movePlanets( struct planet * , double , double );
int planet_motion_analytic(void);

void boundary_r( struct domain * );
void boundary_trans( struct domain * , struct face * , int * , int );
void exchangeData( struct domain * , int );

int get_num_rzFaces( int , int , int );

//
extern void burning(struct domain *,double);
//


void onestep( struct domain * theDomain ,struct poisson * thePoisson, double RK , double dt , int last_step , double global_dt ){


  
   int Nr = theDomain->Nr;
   int Nz = theDomain->Nz;
 
   set_wcell( theDomain );
   adjust_RK_cons( theDomain , RK );

   struct face * theFaces_1 = NULL;
   struct face * theFaces_2 = NULL;

   int NRZ1 = get_num_rzFaces( Nr , Nz , 1 );
   int NRZ2 = get_num_rzFaces( Nr , Nz , 2 );
   int * nfr = (int *) malloc( (NRZ1+1)*sizeof(int) );
   int * nfz = (int *) malloc( (NRZ2+1)*sizeof(int) );;

   phi_flux( theDomain , dt );

   setup_faces( theDomain , &theFaces_1 , nfr , 1 );
   trans_flux( theDomain , theFaces_1 , nfr[NRZ1] , dt , 1 );

   if( Nz > 1 ){
      setup_faces( theDomain , &theFaces_2 , nfz , 2 );
      trans_flux( theDomain , theFaces_2 , nfz[NRZ2] , dt , 2 );
   }

   //self-gravity
   sg_route(theDomain,thePoisson);
   add_source( theDomain , dt );
   move_cells( theDomain , RK , dt );
   if( !planet_motion_analytic() || !last_step ){
      adjust_RK_planets( theDomain , RK );
      movePlanets( theDomain->thePlanets , theDomain->t , dt );
   }
   clean_pi( theDomain );
   calc_dp( theDomain );

//    if(last_step)
//       burning(theDomain,global_dt);

   calc_prim( theDomain );

   if( last_step ){
      AMR( theDomain );
   }

   boundary_trans( theDomain , theFaces_1 , nfr , 1 );
   exchangeData( theDomain , 0 );
   if( Nz > 1 ){
      int Periodic = theDomain->theParList.Z_Periodic;
      if( !Periodic ) boundary_trans( theDomain , theFaces_2 , nfz , 2 );
      exchangeData( theDomain , 1 );
   }

   free( nfr );
   free( nfz );

   if( theFaces_1 ) free( theFaces_1 );
   if( theFaces_2 ) free( theFaces_2 );

   set_wcell( theDomain );

}

