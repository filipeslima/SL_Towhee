#include <stdio.h>
#include <stdlib.h>

/* From "Good Parameters and Implementations for Combined Multiple Recursive
   Random Number Generators," Pierre L'Ecuyer, Operations Research 47 (1),
   159-164 (1999).
*/

/* Initialization follows RNGs in Gnu Scientific Library */

#define LCG(x) ((69069 * (x) + 1) & 0xffffffffUL)

typedef unsigned long UL;

#define norm 2.328306549295728e-10
#define m1   4294967087.0
#define m2   4294944443.0
#define a12     1403580.0
#define a13n     810728.0
#define a21      527612.0
#define a23n    1370589.0

double  s10, s11, s12, s20, s21, s22;

void mrgini(UL seed) {
  UL s, q = LCG(seed ? seed : 4357);
  s  = (q & 0xffff0000UL);        q = LCG(q);
  s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
  s10 = (double) (s<m1 ? s : s-m1);
  s  = (q & 0xffff0000UL);        q = LCG(q);
  s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
  s11 = (double) (s<m1 ? s : s-m1);
  s  = (q & 0xffff0000UL);        q = LCG(q);
  s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
  s12 = (double) (s<m1 ? s : s-m1);
  s  = (q & 0xffff0000UL);        q = LCG(q);
  s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
  s20 = (double) (s<m2 ? s : s-m2);
  s  = (q & 0xffff0000UL);        q = LCG(q);
  s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
  s21 = (double) (s<m2 ? s : s-m2);
  s  = (q & 0xffff0000UL);        q = LCG(q);
  s |= (q & 0xffff0000UL) >> 16;  q = LCG(q);
  s22 = (double) (s<m2 ? s : s-m2);
}


/* Return double, uniformly distributed on [0,1), with period approx 2^191.
   Initialized by calling mrgini with nonzero seed (unsigned long). */

double mrg32k3a () {
  long k;
  double p1, p2;
  /* Component 1 */
  p1 = a12 * s11 - a13n * s10;
  k = p1 / m1;   p1 -= k * m1;   if (p1 < 0.0) p1 += m1;
  s10 = s11;   s11 = s12;   s12 = p1;
  /* Component 2 */
  p2 = a21 * s22 - a23n * s20;
  k  = p2 / m2;  p2 -= k * m2;   if (p2 < 0.0) p2 += m2;
  s20 = s21;   s21 = s22;   s22 = p2;
  /* Combination */
  if (p1 <= p2) return ((p1 - p2 + m1) * norm);
  else return ((p1 - p2) * norm);
}


/* Output the RNG state, to enable later restarting */

void mrgsave(int *nrng, unsigned int *irng) {
  unsigned int *tmp;
  *nrng = 12;
  tmp = (unsigned int *) &s10; irng[0] = *tmp; irng[1] = *(tmp+1);
  tmp = (unsigned int *) &s11; irng[2] = *tmp; irng[3] = *(tmp+1);
  tmp = (unsigned int *) &s12; irng[4] = *tmp; irng[5] = *(tmp+1);
  tmp = (unsigned int *) &s20; irng[6] = *tmp; irng[7] = *(tmp+1);
  tmp = (unsigned int *) &s21; irng[8] = *tmp; irng[9] = *(tmp+1);
  tmp = (unsigned int *) &s22; irng[10] = *tmp; irng[11] = *(tmp+1);
}


/* Input the saved RNG state, to restart */

void mrgrestart(int *nrng, unsigned int *irng) {
  double tmp;
  if (*nrng != 12) {
    printf("In mrgrestart, *nrng has wrong value: %d\n", *nrng);
    exit(1);
  }
  *((unsigned int *) &tmp) = irng[0]; *(1+(unsigned int *) &tmp) = irng[1]; s10 = tmp;
  *((unsigned int *) &tmp) = irng[2]; *(1+(unsigned int *) &tmp) = irng[3]; s11 = tmp;
  *((unsigned int *) &tmp) = irng[4]; *(1+(unsigned int *) &tmp) = irng[5]; s12 = tmp;
  *((unsigned int *) &tmp) = irng[6]; *(1+(unsigned int *) &tmp) = irng[7]; s20 = tmp;
  *((unsigned int *) &tmp) = irng[8]; *(1+(unsigned int *) &tmp) = irng[9]; s21 = tmp;
  *((unsigned int *) &tmp) = irng[10]; *(1+(unsigned int *) &tmp) = irng[11]; s22 = tmp;
}
