#include <stdio.h>
#include <stdlib.h>
#include "preproc.h"

/* Return random number uniformly distributed on [0,1), from a generator
   selected in a towhee_input file from several choices:

   random_luxlevel   generator     
		  		 
        0-4	     RANLUX        
          5	     DX-1597-2-7   
          7	     KISS          
          8	     MRG32k3a      

   These are documented in their own files (which, due to file parameters and
   variables, should not be combined).  Initialization occurs when rnginit_ is
   called.
*/

typedef unsigned int UI;

/* Stub that should not get called. */

double twh_premature(void) {
  printf("Error: called twh_random prior to rnginit and rngrestart.\n");
  exit(2);
  return(911.0);
}
/* global rngpointer */
double (*twh_rngpointer) (void) = twh_premature;

/* random number generator call */
double twh_random_() {
  double value;
  do
    value = twh_rngpointer();
  while (value == 0.);
  return value;
}

/* Code to demonstrate correctness of initialization, save, and restart
   routines. */
/*
int main() {
  void rnginit_(int *, int *);
  void rngsave_(int *, UI *);
  void rngrestart_(int *, UI *);
  int nrng;
  UI irng[1600];
  int rngchoice, seed, i;
  FILE *fp;
  for (rngchoice=0; rngchoice<=8; rngchoice++) {
    seed =  17;
    rnginit_(&rngchoice, &seed);
    for (i=0; i<20; i++)
      printf("%2d %20.16f\n", i, twh_random_());
    rnginit_(&rngchoice, &seed);
    for (i=0; i<10; i++)
      printf("%2d %20.16f\n", i, twh_random_());
    irng[0] = rngchoice;
    rngsave_(&nrng, irng);
    fp=fopen("rngstore.dat","w");
    fprintf(fp, "%d ", nrng);
    for (i=0; i<nrng; i++)
      fprintf(fp, "%d ", irng[i]);
    fclose(fp);
    seed = 117;
    rnginit_(&rngchoice, &seed);
    for (i=0; i<1600; i++) irng[i]=0;
    fp=fopen("rngstore.dat","r");
    fscanf(fp, "%d", &nrng);
    for (i=0; i<nrng; i++)
      fscanf(fp, "%d", irng+i);
    fclose(fp);
    rngrestart_(&nrng, irng);
    for (i=0; i<10; i++)
      printf("%2d %20.16f\n", i+10, twh_random_());
  }
}
*/

/* When rnginit_ is called, twh_rngpointer is set to the chosen random number
   generator, and the generator is initialized. */

void rnginit_(int *rngchoice, int *seed) {
  void rluxgo_(int *,int *,int *,int *);
  void  dxini(int);
  void kissinit_(int *);
  void mrgini(int);
  void twh_random_luxlevel_(int *, int *);
  double dranlux_(void);
  double rngdebug_(void);
  double dx1597(void);
  double kiss99_(void);
  double mrg32k3a(void);
  int zero = 0;
  int get = GLB_GET;
  int luxlevel;
  switch (*rngchoice) {
  case RNG_RANLUX:
    twh_rngpointer = dranlux_;
    twh_random_luxlevel_(&get,&luxlevel);
    rluxgo_(&luxlevel,seed,&zero,&zero);
    return;
  case RNG_DEBUG:
    twh_rngpointer = rngdebug_;
    return;
  case RNG_DX_1597_2_7:
    twh_rngpointer = dx1597;
    dxini(*seed);
    return;
  case RNG_KISS99:
    twh_rngpointer = kiss99_;
    kissinit_(seed);
    return;
  case RNG_MRG32K3A:
    twh_rngpointer = mrg32k3a;
    mrgini(*seed);
    return;
  default:
    printf("In rnginit, unrecognized choice of RNG, rngchoice: %d\n", *rngchoice);
    exit(1);
  }
}

/* Save state for restart */

void rngsave_(int * rngchoice, int *nrng, UI *irng) {
  void rluxut_(UI *);
  void dxsave(int *,UI *);
  void kisssave_(int *, UI *);
  void mrgsave(int *,UI *);
  switch (*rngchoice) {
  case RNG_RANLUX:
    rluxut_(irng+1);
    *nrng = 26;
    return;
  case RNG_DEBUG:
    *nrng = 0;
    return;
  case RNG_DX_1597_2_7:
    dxsave(nrng, irng);
    return;
  case RNG_KISS99:
    kisssave_(nrng, irng);
    return;
  case RNG_MRG32K3A:
    mrgsave(nrng, irng);
    return;
  default:
    printf("In rngsave, unrecognized choice of RNG, rngchoice: %d\n", *rngchoice);
    exit(1);
  }
}

/* Restart with saved state */

void rngrestart_(int * rngchoice, int *nrng, UI *irng) {
  void rluxin_(UI *);
  void dxrestart(int *,UI *);
  void kissrestart_(int *, UI *);
  void mrgrestart(int *,UI *);
  double dranlux_(void);
  double dx1597(void);
  double kiss99_(void);
  double mrg32k3a(void);
  switch (*rngchoice) {
  case RNG_RANLUX:
    twh_rngpointer = dranlux_;
    *nrng = *nrng -1;
    rluxin_(irng+1);
    return;
  case RNG_DX_1597_2_7:
    twh_rngpointer = dx1597;
    dxrestart(nrng, irng);
    return;
  case RNG_KISS99:
    twh_rngpointer = kiss99_;
    kissrestart_(nrng, irng);
    return;
  case RNG_MRG32K3A:
    twh_rngpointer = mrg32k3a;
    mrgrestart(nrng, irng);
    return;
  default:
    printf("In rngrestart, unrecognized choice of RNG, rngchoice: %d\n", *rngchoice);
    exit(1);
  }
}
