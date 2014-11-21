#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "../SG.h"
#include "../paul.h"

double sg_get_r(struct poisson* thePoisson,int i){
	int p=thePoisson->r_rank;
	return (thePoisson->N0_r[p]+i+0.5)*thePoisson->dr;
}

