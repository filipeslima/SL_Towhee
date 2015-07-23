/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program
* Copyright (C) 2006-2009 Marcus G. Martin
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
  last modified 05-02-2013 by M.G. Martin
  handles all of the global variables for the c routines
*/
#include "preproc.h"
#include "globalc.h"

/* global variable declarations */
double time_array[MAXTIMERS];
/* global pointers */
int * integer_pointer[MAX_INT_POINTERS];
int integer_2dregkey[MAX_INT_POINTERS];
int integer_3dregkey[MAX_INT_POINTERS];
double * double_pointer[MAX_DOUBLE_POINTERS];
int double_2dregkey[MAX_DOUBLE_POINTERS];
#if SAFE_COMPARE
int integer_bounds[MAX_INT_POINTERS];
int double_bounds[MAX_DOUBLE_POINTERS];
#endif
void twh_time_array(int flag, int timer, double * dptr)
{
#ifdef USEMPI
  if ( flag == GLB_SET ) {
    time_array[timer] = *dptr;
  }
  else if ( flag == GLB_GET ) {
    *dptr = time_array[timer];
  }
  else if ( flag == GLB_INCR ) {
    time_array[timer] = time_array[timer] + *dptr;
  }
  else {
    printf("TIME_ARRAY: unknown flag %d \n",flag);
    exit(1);
  }
#endif
  return;
}

/* allocation assistance routines */

/*
 * allocates vector size [M] of given elements size.  If malloc fails,
 * exits.  Returns type *void which must be cast to proper type.
 * Example:
     pi = (int*) twh_allocateVector(sizeof(int), N+1);
     
     pi[x] = foo
     free(pi)
 */
void *twh_allocateVector(size_t size, int M) {
  void *x;

  /*It appears that sometimes memory is used before being set.
    twh_com(GLB_GET, ..., xcmi, ycmi, zcmi) in initconf.F appears
    to be one such case. Towhee expects that memory will be zeroed before use.
    calloc() zeroes memory.
  */
  if (!(x = calloc(M,size))) {
    fprintf(stderr, "Fatal memory allocation error\n");
    exit(1);
  }
  return x;
}

/*
 * allocates matrix size [M][N], of given elements size.  If malloc fails,
 * exits.  Returns type **void which must be cast to proper type.
 *
 * Example:
 *  int **ptr = (int**) twh_allocate2dMatrix(sizeof(int), M, N);
 *  double **lv = (double**) twh_allocate2dMatrix(sizeof(double), M, N);
 *  for (i = 0; i < N; i++) {
 *    for (l = 0; l < M; l++) {
 *     ptr[l][i] = something;
 *     lv[l][i] = something else;
 *  }}
 *
 * free all variables mallocked here.
 * for (i = 0; i < M; i++) {
 *   free(lv[i]);
 *   free(ptr[i]);
 * }
 * free(lv);
 * free(ptr);
 *
 */

void **twh_allocate2dMatrix(size_t size, int M, int N) {
  void **x;
  int i;

  if (!(x = (void**) malloc(sizeof(void*) * M))) {
    fprintf(stderr, "Fatal memory allocation error\n");
    exit(1);
  }


  /*Towhee expects that memory will be zeroed before use.
    calloc() zeroes memory
  */
  for (i = 0; i < M; i++)
    if (!(x[i] = (void*) calloc(N,size))) {
      fprintf(stderr, "Fatal memory allocation error\n");
      exit(1);
    }

  return x;
}

/* data template routines */
int twh_1ddouble(int ptr_index, int * flag, const int index, double * value)
/* the template for 1 dimensional double array storage
   originally written 12-04-2008 by M.G. Martin
   last modified 08-17-2011 by M.G. Martin
*/
{
  int count;
  double * ptr_local;
  /*  printf("twh_1ddouble ptr_index: %d flag: %d \n index: %d \n value: %d \n"
      ,ptr_index,*flag,index,*value);
  */
  switch ( *flag ) {
  case GLB_ALLOC:
    /* printf("allocating ptr_index: %d size: %d \n",ptr_index,index); */
    double_pointer[ptr_index] = (double*) twh_allocateVector(sizeof(double), index+1);
#if SAFE_COMPARE
    double_bounds[ptr_index] = index;
#endif
    break;
  case GLB_FREE:
    /* printf("freeing ptr_index: %d pointer value %p \n",ptr_index,double_pointer[ptr_index]); */
    if ( double_pointer[ptr_index] != NULL ) {
      free(double_pointer[ptr_index]);
      double_pointer[ptr_index] = NULL;
    }
    break;
  case GLB_SET:
    ptr_local = double_pointer[ptr_index];
    ptr_local[index] = *value;
#if SAFE_COMPARE
    if ( index > double_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds index: %d in 1ddouble SET ptr_index %d \n"
	     ,index,ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_GET:
    ptr_local = double_pointer[ptr_index];
    *value = ptr_local[index];
#if SAFE_COMPARE
    if ( index > double_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds in 1ddouble GET ptr_index %d flag %d index %d value %f"
	     ,ptr_index,*flag,index,*value);
      printf(" double_bounds[ptr_index] %d \n",double_bounds[ptr_index]);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_GET_TRIPLE:
    ptr_local = double_pointer[ptr_index];
    value[0] = ptr_local[index];
    value[1] = ptr_local[index+1];
    value[2] = ptr_local[index+2];
#if SAFE_COMPARE
    if ( index+2 > double_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds in 1ddouble GET_TRIPLE ptr_index %d \n"
	     ,ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_SET_TRIPLE:
    ptr_local = double_pointer[ptr_index];
    ptr_local[index] = value[0];
    ptr_local[index+1] = value[1];
    ptr_local[index+2] = value[2];
#if SAFE_COMPARE
    if ( index+2 > double_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds index %d in 1ddouble SET_TRIPLE ptr_index %d \n"
	     ,index,ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_INCR_TRIPLE:
    ptr_local = double_pointer[ptr_index];
    ptr_local[index] += value[0];
    ptr_local[index+1] += value[1];
    ptr_local[index+2] += value[2];
#if SAFE_COMPARE
    if ( index+2 > double_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds index %d in 1ddouble INCR_TRIPLE ptr_index %d \n"
	     ,index,ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_INIT:
    for ( count = 0; count <= index; count++) {
      ptr_local = double_pointer[ptr_index];
      ptr_local[count] = *value;
    }
#if SAFE_COMPARE
    if ( index > double_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds in 1ddouble GET ptr_index %d flag %d index %d value %f \n"
	     ,ptr_index,*flag,index,*value);
      printf("double_bounds[ptr_index] %d \n",double_bounds[ptr_index]);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_DECR:
    ptr_local = double_pointer[ptr_index];
    ptr_local[index] -= *value;
    break;
  case GLB_INCR:
    ptr_local = double_pointer[ptr_index];
    ptr_local[index] += *value;
    break;
  case GLB_SCALE:
    ptr_local = double_pointer[ptr_index];
    ptr_local[index] *= *value;
    break;
  default:
    printf("Unknown flag: %d ",*flag);
    return EXIT_FAILURE;
    break;
  }
  return EXIT_SUCCESS;
}

int twh_2dregdouble(int ptr_index, int * flag, int mindex, int nindex, double * value)
/* the template for 2 dimensional regular double array storage
   a regular array is m by n in size
   originally written 06-19-2009 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int key,newindex,errorcode;
  /* only need to set the index, and perhaps the 2dregkey, and then call the 1d data structure */
  switch ( *flag ) {
  case GLB_ALLOC:
    newindex = (mindex+1)*(nindex+1);
    /* store the size of the nindex in double_2dregkey for future use */
    double_2dregkey[ptr_index] = nindex+1;
    break;
  default:
    key = double_2dregkey[ptr_index];
    newindex = mindex*key + nindex;
    break;
  }
  /* call the 1d storage array using the new index */
  errorcode = twh_1ddouble(ptr_index, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {exit(errorcode);}
  /* if we made it this far then exit normally */
  return EXIT_SUCCESS;
}

int twh_1dinteger(int ptr_index, int * flag, const int index, int * value)
/* the template for 1 dimensional integer array storage
   originally written 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int count;
  int * ptr_local;
  /*    printf("twh_1dinteger ptr_index: %d flag: %d \n index: %d \n value: %d \n"
	 ,ptr_index,*flag,index,*value);
  */
  switch ( *flag ) {
  case GLB_ALLOC:
    /*    printf("allocating ptr_index: %d size: %d \n",ptr_index,index); */
    integer_pointer[ptr_index] = (int*) twh_allocateVector(sizeof(int), index+1);
#if SAFE_COMPARE
    integer_bounds[ptr_index] = index+1;
#endif
    break;
  case GLB_FREE:
    /* printf("freeing ptr_index: %d pointer value %p \n",ptr_index,integer_pointer[ptr_index]); */
    if ( integer_pointer[ptr_index] != NULL ) {
      free(integer_pointer[ptr_index]);
      integer_pointer[ptr_index] = NULL;
    }
    break;
  case GLB_SET:
    ptr_local = integer_pointer[ptr_index];
    ptr_local[index] = *value;
#if SAFE_COMPARE
    if ( index > integer_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds index: %d in 1dint SET ptr_index %d \n"
	     ,index,ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_GET:
    ptr_local = integer_pointer[ptr_index];
    *value = ptr_local[index];
#if SAFE_COMPARE
    if ( index > integer_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds in 1dinteger GET ptr_index %d \n",ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_INIT:
    for ( count = 0; count <= index; count++) {
      ptr_local = integer_pointer[ptr_index];
      ptr_local[count] = *value;
    }
#if SAFE_COMPARE
    if ( index > integer_bounds[ptr_index] || index < 0 ) {
      printf("out of bounds index: %d in 1dinteger INIT ptr_index %d \n",index,ptr_index);
      exit(EXIT_FAILURE);
    }
#endif
    break;
  case GLB_DECR:
    ptr_local = integer_pointer[ptr_index];
    ptr_local[index] -= *value;
    break;
  case GLB_INCR:
    ptr_local = integer_pointer[ptr_index];
    ptr_local[index] += *value;
    break;
  default:
    printf("Unknown flag: %d ",*flag);
    return EXIT_FAILURE;
    break;
  }
  return EXIT_SUCCESS;
}

int twh_2dreginteger(int ptr_index, int * flag, int mindex, int nindex, int * value)
/* the template for 2 dimensional regular integer array storage
   a regular array is m by n in size
   originally written 06-24-2009 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int key,newindex,errorcode;
  /* only need to set the index, and perhaps the 2dregkey, and then call the 1d data structure */
  switch ( *flag ) {
  case GLB_ALLOC:
    newindex = (mindex+1)*(nindex+1);
    /* store the size of the nindex in integer_2dregkey for future use */
    integer_2dregkey[ptr_index] = nindex+1;
    break;
  default:
    key = integer_2dregkey[ptr_index];
    newindex = mindex*key + nindex;
    break;
  }
  /* call the 1d storage array using the new index */
  errorcode = twh_1dinteger(ptr_index, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {exit(errorcode);}
  /* if we made it this far then exit normally */
  return EXIT_SUCCESS;
}

int twh_3dreginteger(int ptr_index, int * flag, int mindex, int nindex, int oindex, int * value)
/* the template for 3 dimensional regular integer array storage
   a regular array is m by n by o in size
   originally written 06-25-2009 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int keyn,keyo,newindex,errorcode;
  /* only need to set the index, and perhaps the 3dregkey, 2dregkey, and then call the
     1d data structure */
  switch ( *flag ) {
  case GLB_ALLOC:
    newindex = (mindex+1)*(nindex+1)*(oindex+1);
    /* store the size of the nindex in integer_2dregkey and 3dregkey for future use */
    integer_2dregkey[ptr_index] = nindex+1;
    integer_3dregkey[ptr_index] = oindex+1;
    break;
  default:
    keyn = integer_2dregkey[ptr_index];
    keyo = integer_3dregkey[ptr_index];
    newindex = keyo*(mindex*keyn + nindex) + oindex;
    break;
  }
  /* call the 1d storage array using the new index */
  errorcode = twh_1dinteger(ptr_index, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {exit(errorcode);}
  /* if we made it this far then exit normally */
  return EXIT_SUCCESS;
}

/* data storage */

#ifdef INTEL_VISUAL_FORTRAN
#define twh_coordfield_ TWH_COORDFIELD
#define twh_coordstorage_ TWH_COORDSTORAGE
#define twh_coordtemp_ TWH_COORDTEMP
#define twh_cubeletweight_ TWH_CUBELETWEIGHT
#define twh_gyration_ TWH_GYRATION
#define twh_eam_rho_real_ TWH_EAM_RHO_REAL
#define twh_eam_rho_temp_ TWH_EAM_RHO_TEMP
#define twh_rcmu_ TWH_RCMU
#define twh_tmmc_weight_ TWH_TMMC_WEIGHT
#define twh_v_semigrand_ TWH_V_SEMIGRAND
#define twh_acnvol_ TWH_ACNVOL
#define twh_acsvol_ TWH_ACSVOL
#define twh_ewald_kmax_ TWH_EWALD_KMAX
#define twh_glist_ TWH_GLIST
#define twh_globalpos_ TWH_GLOBALPOS
#define twh_growfrom_ TWH_GROWFROM
#define twh_grownum_ TWH_GROWNUM
#define twh_growprev_ TWH_GROWPREV
#define twh_growvalidcount_  TWH_GROWVALIDCOUNT
#define twh_logical_exist_ TWH_LOGICAL_EXIST
#define twh_logical_exsched_ TWH_LOGICAL_EXSCHED
#define twh_logical_moveme_ TWH_LOGICAL_MOVEME
#define twh_logical_periodic_ TWH_LOGICAL_PERIODIC
#define twh_moltyp_ TWH_MOLTYP
#define twh_nboxi_ TWH_NBOXI
#define twh_rmvol_ TWH_RMVOL
#define twh_arbcmofield_ TWH_ARBCMOFIELD
#define twh_comfield_ TWH_COMFIELD
#define twh_comtempfield_ TWH_COMTEMPFIELD
#endif

/* 2d regular double array storage */

void twh_acncell_(int * flag, int * imove, int * ivector, double * value)
/* total number of attempted cell volume moves for each box combination and each vector dimension
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imove -1;
  nindex = *ivector -1;
  errorcode = twh_2dregdouble(PNT_ACNCELL, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acncell \n");
    exit(errorcode);
  }
  return;
}

void twh_acscell_(int * flag, int * imove, int * ivector, double * value)
/* total number of accepted cell volume moves for each box combination and each vector dimension
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imove -1;
  nindex = *ivector -1;
  errorcode = twh_2dregdouble(PNT_ACSCELL, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acscell \n");
    exit(errorcode);
  }
  return;
}

void twh_acncomp_(int * flag, int * imolty, int * ibox, double * value)
/* accumulator for attempted composite moves
   rewritten from Fortran 06-19-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACNCOMP, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acncomp \n");
    exit(errorcode);
  }
  return;
}

void twh_acscomp_(int * flag, int * imolty, int * ibox, double * value)
/* total count of accepted composite moves
   originally written in Fortran 03-25-2006 by MAW
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACSCOMP, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acscomp \n");
    exit(errorcode);
  }
  return;
}

void twh_acnrot_(int * flag, int * imolty, int * ibox, double * value)
/* accumulator for the rotation move
   rewritten from Fortran 06-18-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACNROT, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acnrot \n");
    exit(errorcode);
  }
  return;
}

void twh_acsrot_(int * flag, int * imolty, int * ibox, double * value)
/* total count of accepted rotation about the COM moves
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACSROT, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acsrot \n");
    exit(errorcode);
  }
  return;
}

void twh_acnswitch_(int * flag, int * ipair, int * ibox, double * value)
/* total count of attempted switch moves
   originally written in Fortran 06-07-2008 by I.A. Hijazi
   rewritten from Fortran 06-19-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *ipair -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACNSWITCH, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acnswitch \n");
    exit(errorcode);
  }
  return;
}

void twh_acsswitch_(int * flag, int * ipair, int * ibox, double * value)
/* total count of accepted switch moves
   originally written in Fortran 06-07-2008 by I.A. Hijazi
   rewritten from Fortran 06-19-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *ipair -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACSSWITCH, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acnswitch \n");
    exit(errorcode);
  }
  return;
}

void twh_acntraa_(int * flag, int * imolty, int * ibox, double * value)
/* total count of attempted single-atom translation moves
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-20-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACNTRAA, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acntraa \n");
    exit(errorcode);
  }
  return;
}

void twh_acstraa_(int * flag, int * imolty, int * ibox, double * value)
/* total count of accepted COM translation moves
   originally written 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-20-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACSTRAA, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acstraa \n");
    exit(errorcode);
  }
  return;
}

void twh_acntrac_(int * flag, int * imolty, int * ibox, double * value)
/* total count of attempted COM translation moves
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-20-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACNTRAC, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acntrac \n");
    exit(errorcode);
  }
  return;
}

void twh_acstrac_(int * flag, int * imolty, int * ibox, double * value)
/* total count of accepted single-atom translation moves
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-20-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_ACSTRAC, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acstrac \n");
    exit(errorcode);
  }
  return;
}

void twh_bacell_(int * flag, int * imove, int * ivector, double * value)
/* running count of accepted cell moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imove -1;
  nindex = *ivector -1;
  errorcode = twh_2dregdouble(PNT_BACELL, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_bacell \n");
    exit(errorcode);
  }
  return;
}

void twh_bncell_(int * flag, int * imove, int * ivector, double * value)
/* running count of attempted cell moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imove -1;
  nindex = *ivector -1;
  errorcode = twh_2dregdouble(PNT_BNCELL, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_bncell \n");
    exit(errorcode);
  }
  return;
}

void twh_barot_(int * flag, int * imolty, int * ibox, double * value)
/* running count of accepted rotation about COM moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-23-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_BAROT, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_barot \n");
    exit(errorcode);
  }
  return;
}

void twh_bnrot_(int * flag, int * imolty, int * ibox, double * value)
/* running count of attempted rotation about COM moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-23-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_BNROT, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_bnrot \n");
    exit(errorcode);
  }
  return;
}

void twh_batraa_(int * flag, int * imolty, int * ibox, double * value)
/* running count of accepted single-atom translation moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_BATRAA, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_batraa \n");
    exit(errorcode);
  }
  return;
}

void twh_bntraa_(int * flag, int * imolty, int * ibox, double * value)
/* running count of attempted single-atom translation moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_BNTRAA, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_bntraa \n");
    exit(errorcode);
  }
  return;
}

void twh_batrac_(int * flag, int * imolty, int * ibox, double * value)
/* running count of accepted COM translation moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_BATRAC, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_batrac \n");
    exit(errorcode);
  }
  return;
}

void twh_bntrac_(int * flag, int * imolty, int * ibox, double * value)
/* running count of attempted COM translation moves since the last update
   originally written in Fortran 03-21-2006 by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_BNTRAC, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_bntrac \n");
    exit(errorcode);
  }
  return;
}

void twh_blockvalue_(int * flag, int * iprop, int * iblock, double * value)
/* averaging invformation for various properties used to compute block averages
   rewritten from Fortran 02-03-2011 by M.G. Martin
   last modified 02-03-2011 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *iprop -1;
  nindex = *iblock -1;
  errorcode = twh_2dregdouble(PNT_BLOCKVALUE, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_blockvalue \n");
    exit(errorcode);
  }
  return;
}

void twh_c_matrix_(int * flag, int * ichain, int * iparam, double * value)
/* TMMC collection matrix
   originally written in Fortran 10-13-2008 by M.G. Martin
   rewritten from Fortran 06-23-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  /* ichain ranges from 0 to the maximum so no shift */
  mindex = *ichain;
  /* iparam ranges from -1 to 1 so shift up 1 */
  nindex = *iparam +1;
  errorcode = twh_2dregdouble(PNT_CMATRIX, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_c_matrix \n");
    exit(errorcode);
  }
  return;
}

void twh_rmcomrot_(int * flag, int * imolty, int * ibox, double * value)
/* the maximum rotational displacement for usein the composite move for each molecule type
   in each box
   originally written 09-24-2008 in Fortran by M.G. Martin
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_RMCOMROT, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rmcomrot \n");
    exit(errorcode);
  }
  return;
}

void twh_rmcomtra_(int * flag, int * imolty, int * ibox, double * value)
/* the maximum translations displacement for usein the composite move for each molecule
   type in each box
   originally written 09-24-2008 in Fortran by M.G. Martin
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_RMCOMTRA, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rmcomtra \n");
    exit(errorcode);
  }
  return;
}

void twh_rmtraa_(int * flag, int * imolty, int * ibox, double * value)
/* the maximum single-atom translational displacement for each molecule type in each box
   originally written 09-23-2008 in Fortran by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_RMTRAA, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rmtraa \n");
    exit(errorcode);
  }
  return;
}

void twh_rmtrac_(int * flag, int * imolty, int * ibox, double * value)
/* the maximum COM translational displacement for each molecule type in each box
   originally written 09-24-2008 in Fortran by M.G. Martin
   rewritten from Fortran 06-21-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_RMTRAC, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rmtrac \n");
    exit(errorcode);
  }
  return;
}

void twh_rmrot_(int * flag, int * imolty, int * ibox, double * value)
/* the maximum rotational displacement for each molecule type in each box
   originally written in Fortran 09-24-2008 by M.G. Martin
   rewritten from Fortran 06-22-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{ 
  int errorcode,mindex,nindex;
  mindex = *imolty -1;
  nindex = *ibox -1;
  errorcode = twh_2dregdouble(PNT_RMROT, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rmrot \n");
    exit(errorcode);
  }
  return;
}

/* 1d double array storage */
void twh_arbcmofield_(int * flag, int * index, double * value)
/* field molecule center of mass coordinate storage of triples
   originally written 06-23-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* no shifting required */
  newindex = *index;
  errorcode = twh_1ddouble(PNT_ARBCMOFIELD, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_arbcmofield \n");
    exit(errorcode);
  }
  return;
}

void twh_comfield_(int * flag, int * index, double * value)
/* center of mass coordinate storage of triples
   originally written 06-24-2009 by M.G. Martin
   last modified 06-24-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* no shifting required */
  newindex = *index;
  errorcode = twh_1ddouble(PNT_COMFIELD, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_comfield \n");
    exit(errorcode);
  }
  return;
}

void twh_comtempfield_(int * flag, int * index, double * value)
/* center of mass temporary coordinate storage of triples
   originally written 06-24-2009 by M.G. Martin
   last modified 06-24-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* no shifting required */
  newindex = *index;
  errorcode = twh_1ddouble(PNT_COMTEMPFIELD, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_comtempfield \n");
    exit(errorcode);
  }
  return;
}

void twh_coordfield_(int * flag, int * index, double * value)
/* field coordinate storage of triples
   rewritten from Fortran 12-06-2008 by M.G. Martin
   last modified 12-06-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_COORDFIELD, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_coordfield \n");
    exit(errorcode);
  }
  return;
}

void twh_coordstorage_(int * flag, int * index, double * value)
/* real coordinate storage of triples
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_COORDSTORAGE, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_coordstorage \n");
    exit(errorcode);
  }
  return;
}

void twh_coordtemp_(int * flag, int * index, double * value)
/* temporary coordinate storage of triples
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_COORDTEMP, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_coordtemp \n");
    exit(errorcode);
  }
  return;
}

void twh_cubeletweight_(int * flag, int * index, double * value)
/* the normalized biasing weight for each cubelet
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_CUBELETWEIGHT, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_cubeletweight \n");
    exit(errorcode);
  }
  return;
}

void twh_gyration_(int * flag, int * index, double * value)
/* the current radius of gyration for each molecule
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 12-04-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_GYRATION, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_gyration \n");
    exit(errorcode);
  }
  return;
}

void twh_eam_rho_real_(int * flag, int * index, double * value)
/* the real current embedding density for each molecule in the simulation
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_EAM_RHO_REAL, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_eam_rho_real \n");
    exit(errorcode);
  }
  return;
}

void twh_eam_rho_temp_(int * flag, int * index, double * value)
/* the temporary embedding density for each molecule in the simulation
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index - 1;
  errorcode = twh_1ddouble(PNT_EAM_RHO_TEMP, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_eam_rho_temp \n");
    exit(errorcode);
  }
  return;
}

void twh_rcmu_(int * flag, int * index, double * value)
/* the maximum distance from the center of mass to any atom for
   each molecule in the system
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index - 1;
  errorcode = twh_1ddouble(PNT_RCMU, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rcmu \n");
    exit(errorcode);
  }
  return;
}

void twh_rmvol_(int * flag, int * index, double * value)
/* the maximum volume displacement for each box or box pair
   originally written in Fortran 06-12-2006 by M.G. Martin
   rewritten from Fortran 06-23-2009 by M.G. Martin
   last modified 06-23-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down by 1 */
  newindex = *index - 1;
  errorcode = twh_1ddouble(PNT_RMVOL, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_rmvol \n");
    exit(errorcode);
  }
  return;
}

void twh_tmmc_weight_(int * flag, int * index, double * value)
/* TMMC biasing function
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* do not shift index */
  newindex = *index;
  errorcode = twh_1ddouble(PNT_TMMC_WEIGHT, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_tmmc_weight \n");
    exit(errorcode);
  }
  return;
}

void twh_v_semigrand_(int * flag, int * index, double * value)
/* TMMC semigrand average potential energy
   rewritten from Fortran 12-05-2008 by M.G. Martin
   last modified 12-05-2008 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* do not shift index */
  newindex = *index;
  errorcode = twh_1ddouble(PNT_V_SEMIGRAND, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_v_semigrand \n");
    exit(errorcode);
  }
  return;
}

void twh_wrap_foreign_energy_(int * flag, int * index, double * value)
/* The foreign energy
   originally written 07-02-2009 by M.G. Martin
   last modified 07-02-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_WRAP_FOREIGN_ENERGY, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_wrap_foreign_energy \n");
    exit(errorcode);
  }
  return;
}


void twh_wrap_foreign_lambda_lj_(int * flag, int * index, double * value)
/* The foreign lambda_lj (Lennard-Jones - generally the nonbonded noncoulombic)
   originally written 07-02-2009 by M.G. Martin
   last modified 07-02-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_WRAP_FOREIGN_LAMBDA_LJ, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_wrap_foreign_lambda_lj \n");
    exit(errorcode);
  }
  return;
}


void twh_wrap_foreign_lambda_c_(int * flag, int * index, double * value)
/* The foreign lambda_c (Coulombic parts)
   originally written 07-02-2009 by M.G. Martin
   last modified 07-02-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  /* shift index down 1 */
  newindex = *index -1;
  errorcode = twh_1ddouble(PNT_WRAP_FOREIGN_LAMBDA_C, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_wrap_foreign_lambda_c \n");
    exit(errorcode);
  }
  return;
}


/* 1d integer array storage */
void twh_acnvol_(int * flag, int * index, int * value)
/* total number of attempted volume moves for each box combination
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_ACNVOL, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acnvol \n");
    exit(errorcode);
  }
  return;
}

void twh_acsvol_(int * flag, int * index, int * value)
/* total number of successful volume moves for each box combination
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_ACSVOL, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_acnvol \n");
    exit(errorcode);
  }
  return;
}

void twh_ewald_kmax_(int * flag, int * index, int * value)
/* the maximum number of k-vectors in each box of the Ewald sum
   kmax controls the total number of reciprocal vectors
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_EWALD_KMAX, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_ewald_kmax \n");
    exit(errorcode);
  }
  return;
}

void twh_glist_(int * flag, int * index, int * value)
/* temp array to keep track of atoms in a molecule that still
   need something done to them
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_GLIST, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_glist \n");
    exit(errorcode);
  }
  return;
}

void twh_globalpos_(int * flag, int * index, int * value)
/* globalpos has the index used to located the first atom of each
   molecule in the data structures
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_GLOBALPOS, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_globalpos \n");
    exit(errorcode);
  }
  return;
}

void twh_growfrom_(int * flag, int * index, int * value)
/* contains the from atom for each step of a CBMC move
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_GROWFROM, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_growfrom \n");
    exit(errorcode);
  }
  return;
}

void twh_grownum_(int * flag, int * index, int * value)
/* contains the number of atoms to be grown for each step of a CBMC move
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_GROWNUM, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_grownum \n");
    exit(errorcode);
  }
  return;
}

void twh_growprev_(int * flag, int * index, int * value)
/* contains the prev atom for each step of a CBMC move
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_GROWPREV, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_growprev \n");
    exit(errorcode);
  }
  return;
}

void twh_growvalidcount_(int * flag, int * index, int * value)
/* contains the number of valid starting atoms for each molecule type 
   when using a CBMC move
   rewritten from Fortran 12-03-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_GROWVALIDCOUNT, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_growvalidcount \n");
    exit(errorcode);
  }
  return;
}

void twh_logical_exist_(int * flag, int * index, int * value)
/* temporary variable for keeping track of which units have been
   used for some purpose: in CBMC this is the beads that already exist
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_LOGICAL_EXIST, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_logical_exist \n");
    exit(errorcode);
  }
  return;
}

void twh_logical_exsched_(int * flag, int * index, int * value)
/* scheduling logical for use during CBMC growths
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_LOGICAL_EXSCHED, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_logical_exsched \n");
    exit(errorcode);
  }
  return;
}

void twh_logical_moveme_(int * flag, int * index, int * value)
/* scheduling logical for use during atomshift moves
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_LOGICAL_MOVEME, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_logical_moveme \n");
    exit(errorcode);
  }
  return;
}

void twh_logical_periodic_(int * flag, int * index, int * value)
/* logical as to whether a molecule currently wraps through the
   periodic boundaries in a manner that makes it infinite
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_LOGICAL_PERIODIC, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_logical_periodic \n");
    exit(errorcode);
  }
  return;
}

void twh_moltyp_(int * flag, int * index, int * value)
/* molecule type for each chain in the simulation
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_MOLTYP, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_moltyp \n");
    exit(errorcode);
  }
  return;
}

void twh_nboxi_(int * flag, int * index, int * value)
/* the current simulation box for each chain in the simulation
   rewritten from Fortran 12-04-2008 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,newindex;
  newindex = *index - 1;
  errorcode = twh_1dinteger(PNT_NBOXI, flag, newindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_nboxi \n");
    exit(errorcode);
  }
  return;
}

/* 2d integer array storage */
void twh_parall_(int * flag, int * itype, int * ipoint, int * value)
/* the list of all chains of a given type in any box
   originally written in Fortran 11-01-2007 by M.G. Martin
   rewritten from Fortran 06-24-2009 by M.G. Martin
   last modified 06-24-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex;
  mindex = *itype - 1;
  nindex = *ipoint -1;
  errorcode = twh_2dreginteger(PNT_PARALL, flag, mindex, nindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_parall \n");
    exit(errorcode);
  }
  return;
}

/* 3d integer array storage */
void twh_chainlist_(int * flag, int * ipoint, int * ibox, int * itype, int * value)
/* the list of chains of each molecule type in each simulation box
   originally written in Fortran 11-01-2007 by M.G. Martin
   rewritten from Fortran 06-25-2009 by M.G. Martin
   last modified 06-25-2009 by M.G. Martin
*/
{
  int errorcode,mindex,nindex,oindex;
  mindex = *ipoint - 1;
  /* ibox ranges from 0 to maximum */
  nindex = *ibox;
  oindex = *itype - 1;
  errorcode = twh_3dreginteger(PNT_CHAINLIST, flag, mindex, nindex, oindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_chainlist \n");
    exit(errorcode);
  }
  return;
}

void twh_torofcode_(int * flag, int * imolty, int * iunit, int * itor, int * value)
/* the list of one-four scaling codes for each torsion in each molecule type
   originally written 01-30-2011 by M.G. Martin
   last modified 01-30-2011 by M.G. Martin
*/
{
  int errorcode,mindex,nindex,oindex;
  mindex = *imolty - 1;
  nindex = *iunit - 1;
  oindex = *itor - 1;
  errorcode = twh_3dreginteger(PNT_TOROFCODE, flag, mindex, nindex, oindex, value);
  if ( errorcode == EXIT_FAILURE ) {
    printf("in twh_torofcode \n");
    exit(errorcode);
  }
  return;
}
