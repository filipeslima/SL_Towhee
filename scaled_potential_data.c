/*     $Id:  */
/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program
* Copyright (C) 2003-2005 Marcus G. Martin
* see the file license.gpl for the full license information
*
* This program is free software; you can redistribute it and/or
* modify it under the terms of the GNU General Public License
* as published by the Free Software Foundation; either version 2
* of the License, or (at your option) any later version.
*
* This program is distributed in the hope that it will be useful,
* but WITHOUT ANY WARRANTY; without even the implied warranty of
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
* GNU General Public License for more details.
*
* You should have received a copy of the GNU General Public
* License along with this program; if not, write to the Free
* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
* MA  02111-1307, USA.
*/
/* 
   scaled_potential.c holds a dynamically allocated array of 
   foreign lambdas (see scaled_potential.F), as well as related
   functions.  As much for demonstration purposes of dynamic memory
   allocation in Fortran 77 via C as anything else.

   written 04-29-2006 by MAW
   last modified 12-03-2008 by M.G. Martin
*/

#include "globalc.h"

int number;
double **foreign_lambda;
double *testvec;

/* see, http://yolinux.com/TUTORIALS/LinuxTutorialMixingFortranAndC.html
   apparently, functions must end in an underscore, two if underscores exist elsewhere in
   name */
void initializec_(int *newval) {
   printf("howdy!\n");
   printf("memorizing %d \n", *newval);
   number = *newval + 3;
}

void get_number__(int *newval) {
   printf("I remember  %d\n", number);
   *newval = number; 
}

void initialize_foreign_lambda__(int *length) {
  testvec = (double*) twh_allocateVector(sizeof(double), *length+1);

  testvec[1] = 1.;
  testvec[*length] = -1.;
     
  /* pi[x] = foo */
  /* free(pi) */
} 

void get_foreign_lambda__(int *index, double *v) {
  *v = testvec[*index];
} 
