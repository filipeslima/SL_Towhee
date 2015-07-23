#include <stdio.h>
#include <stdlib.h>

/* From "Efficient and Portable Multiple Recursive Generators of Large Order,"
   Lih-Yuan Deng, ACM Transactions on Modeling and Computer Simulation 15 (1),
   1-13 (2005).

   x(i) = (-2^r-2^w)[x(i-t)+x(i-k)] mod (2^31-1)

   last modified 08-02-2001 by M.G. Martin
*/

/* Initialization follows RNGs in Gnu Scientific Library */

#define LCG(x) ((69069 * (x) + 1) & 0xffffffffUL)

typedef unsigned long UL;

#define k 1597
#define t 7
#define r 25
#define w 7
#define p 0x7fffffff  /* 2^31-1 */
#define scale 4.6566128752457969241e-10 /* 1/p */
#define MODP(x) (((x)>>31)+((x)&p))
#define MULMODP(b,x) (((x)>>(31-b))+(((x)<<(b))&p))

static UL state[k];
static int i;

void dxini(UL seed) {
  int j;
  UL s, q = LCG(seed ? seed : 4357);
  for (j=0; j<k; j++) {
    s  = (q & 0xffff0000UL);        q = LCG(q);
    s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
    state[j] = MODP(s);
  }
  i=k-1;
}


/* Return double, uniformly distributed on [0,1), with period approx 2^49507.
   Initialized by calling dxini with nonzero seed (unsigned long). */

double dx1597() {
  UL s;
  int j;
  if (++i >= k) i=0;
  if ((j=i-t)<0) j+=k;
  s = state[i]+state[j];
  s = MODP(s);
  s = (p-MULMODP(r,s))+(p-MULMODP(w,s));
  state[i] = MODP(s);
  return(scale * (double) state[i]);
}


/* Output the RNG state, to enable later restarting */

void dxsave(int *nrng, unsigned int *irng) {
  int j;
  *nrng = k+1;
  for (j=0; j<k; j++) irng[j] = state[j];
  irng[k] = i;
}


/* Input the saved RNG state, to restart */

void dxrestart(int *nrng, unsigned int *irng) {
  int j;
  if (*nrng != k+1) {
    printf("In dxrestart, *nrng has wrong value: %d\n", *nrng);
    exit(1);
  }
  for (j=0; j<k; j++) state[j] = irng[j];
  i = irng[k];
}
  
