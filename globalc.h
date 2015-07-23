/*     $Id: globalc.h,v 1.19 2010/09/04 14:16:49 lperi Exp $ */
/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program     *
* Copyright (C) 2003-2005 Marcus G. Martin                       *
* see the file license.gpl for the full license information      *
*                                                                *
* This program is free software; you can redistribute it and/or  *
* modify it under the terms of the GNU General Public License    *
* as published by the Free Software Foundation; either version 2 *
* of the License, or (at your option) any later version.         *
*                                                                *
* This program is distributed in the hope that it will be useful,*
* but WITHOUT ANY WARRANTY; without even the implied warranty of *
* MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
* GNU General Public License for more details.                   *
*                                                                *
* You should have received a copy of the GNU General Public      *
* License along with this program; if not, write to the Free     *
* Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
* MA  02111-1307, USA.                                           
*/

/*
  globalc.h
  this file contains all of the information that is at the top 
  of every C routine in Towhee
  last modified 07-02-2009 by M.G. Martin
   generic c includes 
*/
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <errno.h>

/* Include general towhee definitions */
#include "preproc.h"

/* C specific preprocessor derectives */
/* maximum number of timers for keeping track of mpi_wtime */
#define MAXTIMERS 5

#define TRUE 1
#define FALSE 0

/* Structures */
/*
 * user_args contains information about command line arguments passed to
 * towhee.  Populated by parse_args, and passed around as a pointer to be
 * processed where needed.
 */
typedef struct {
  char rex_fn[MAXDIRLENGTH];
  char towhee_input_fn[MAXDIRLENGTH];
  int quiet_mode;
  int parstyle;
  int random_seed;
} user_args;

typedef struct {
  int num_nodes;
  int num_steps;
  int num_rounds;
  int num_swaps;
  float lambda_lj[MAX_FOREIGN_LAMBDA];
  float lambda_c[MAX_FOREIGN_LAMBDA];
} rex_params;

/* 
   template or function declarations 
*/
/* extern void calltest_(); */

#ifdef INTEL_VISUAL_FORTRAN
#define towheemainloop_ TOWHEEMAINLOOP
#define twhquantumloop_ TWHQUANTUMLOOP
#define twhendquantum_ TWHENDQUANTUM
#define swrap_ SWRAP
#endif

#ifdef MSVC
#define snprintf _snprintf
#endif

extern void towheemainloop_(const int *,const int *, const int *, const int *);
extern void twhquantumloop_();
extern void twhendquantum_();
void timekeeper(int, int, int);
void clearline (FILE *);
void tramonto_control(int, int);
double tramonto_();
void tramontoloop();
void jobfarm(FILE *, int, int);
int get_rex_directory(char*, const char*, int, int, int);
void initialize_wrapper_strings();
void set_towhee_input_file(const char*);
void set_output_dir(const char*);
int populate_directory_initial(char*, const char*, int);
int populate_directory_restart(char*, const char*, int, int);
void initialize_args( user_args*);
void parse_args(int, char**, user_args*);
int get_parstyle(int, int);
void print_usage(const char*);
int read_rex_params(const char*, rex_params*);

void open_towhee_parallel(const char*, FILE **);
void run_rex_server(user_args*, rex_params*);
void run_rex_client(user_args*, rex_params*, int);
int attempt_exchange(float, int, int, double**, int*, rex_params*, FILE*);
int write_lambda_file(const int*, const rex_params*, const char*, int, 
        const char*);
void timestamp(char*);
void **twh_allocate2dMatrix(size_t, int, int);
void *twh_allocateVector(size_t, int);

int initialize_MPI(int, char**);
void do_none(user_args*);
void do_tramonto(user_args*, int);
void do_jobfarm(user_args*, int);
void do_rex(user_args*, int);
void do_lcao(user_args*, int);

void lcao_control(int, int);

void twh_time_array(int, int, double*);
/* 
   USEMPI only  
*/
#ifdef USEMPI
#include <mpi.h>
#endif
/* 
   USETRAMONTO only 
*/
#ifdef USETRAMONTO
extern dftmain(double *);
#endif

/* data storage */
#ifdef INTEL_VISUAL_FORTRAN
#define twh_acncell_ TWH_ACNCELL
#define twh_acscell_ TWH_ACSCELL
#define twh_acncomp_ TWH_ACNCOMP
#define twh_acscomp_ TWH_ACSCOMP
#define twh_acnrot_ TWH_ACNROT
#define twh_acsrot_ TWH_ACSROT
#define twh_acnswitch_ TWH_ACNSWITCH
#define twh_acsswitch_ TWH_ACSSWITCH
#define twh_acntraa_ TWH_ACNTRAA
#define twh_acstraa_ TWH_ACSTRAA
#define twh_acntrac_ TWH_ACNTRAC
#define twh_acstrac_ TWH_ACSTRAC
#define twh_bacell_ TWH_BACELL
#define twh_bncell_ TWH_BNCELL
#define twh_barot_ TWH_BAROT
#define twh_bnrot_ TWH_BNROT
#define twh_batraa_ TWH_BATRAA
#define twh_bntraa_ TWH_BNTRAA
#define twh_batrac_ TWH_BATRAC
#define twh_bntrac_ TWH_BNTRAC
#define twh_c_matrix_ TWH_C_MATRIX
#define twh_rmcomrot_ TWH_RMCOMROT
#define twh_rmcomtra_ TWH_RMCOMTRA
#define twh_rmtraa_ TWH_RMTRAA
#define twh_rmtrac_ TWH_RMTRAC
#define twh_rmrot_ TWH_RMROT
#define twh_arbcmofield_ TWH_ARBCMOFIELD
#define twh_comfield_ TWH_COMFIELD
#define twh_comtempfield_ TWH_COMTEMPFIELD
#define twh_coordfield_ TWH_COORDFIELD
#define twh_coordstorage_ TWH_COORDSTORAGE
#define twh_coordtemp_ TWH_COORDTEMP
#define twh_cubeletweight_ TWH_CUBELETWEIGHT
#define twh_gyration_ TWH_GYRATION
#define twh_eam_rho_real_ TWH_EAM_RHO_REAL
#define twh_eam_rho_temp_ TWH_EAM_RHO_TEMP
#define twh_rcmu_ TWH_RCMU
#define twh_rmvol_ TWH_RMVOL
#define twh_tmmc_weight_ TWH_TMMC_WEIGHT
#define twh_v_semigrand_ TWH_V_SEMIGRAND
#define twh_wrap_foreign_energy_ TWH_WRAP_FOREIGN_ENERGY
#define twh_wrap_foreign_lambda_lj_ TWH_WRAP_FOREIGN_LAMBDA_LJ
#define twh_wrap_foreign_lambda_c_ TWH_WRAP_FOREIGN_LAMBDA_C
#define twh_acsvol_ TWH_ACSVOL
#define twh_ewald_kmax_ TWH_EWALD_KMAX
#define twh_glist_ TWH_GLIST
#define twh_globalpos_ TWH_GLOBALPOS
#define twh_growfrom_ TWH_GROWFROM 
#define twh_grownum_ TWH_GROWNUM
#define twh_growprev_ TWH_GROWPREV
#define twh_growvalidcount_ TWH_GROWVALIDCOUNT
#define twh_logical_exist_ TWH_LOGICAL_EXIST
#define twh_logical_exsched_ TWH_LOGICAL_EXSCHED
#define twh_logical_moveme_ TWH_LOGICAL_MOVEME
#define twh_logical_periodic_ TWH_LOGICAL_PERIODIC
#define twh_moltyp_ TWH_MOLTYP
#define twh_nboxi_ TWH_NBOXI
#define twh_parall_ TWH_PARALL
#define twh_chainlist_ TWH_CHAINLIST
#define twh_torofcode_ TWH_TOROFCODE
#endif
int twh_1ddouble(int, int *, const int, double *);
int twh_2dregdouble(int, int *, int, int, double *);
int twh_1dinteger(int, int *, const int, int *);
int twh_2dreginteger(int, int *, int, int, int *);
int twh_3dreginteger(int, int *, int, int, int, int *);
void twh_acncell_(int *, int *, int *, double *);
void twh_acscell_(int *, int *, int *, double *);
void twh_acncomp_(int *, int *, int *, double *);
void twh_acscomp_(int *, int *, int *, double *);
void twh_acnrot_(int *, int *, int *, double *);
void twh_acsrot_(int *, int *, int *, double *);
void twh_acnswitch_(int *, int *, int *, double *);
void twh_acsswitch_(int *, int *, int *, double *);
void twh_acntraa_(int *, int *, int *, double *);
void twh_acstraa_(int *, int *, int *, double *);
void twh_acntrac_(int *, int *, int *, double *);
void twh_acstrac_(int *, int *, int *, double *);
void twh_bacell_(int *, int *, int *, double *);
void twh_bncell_(int *, int *, int *, double *);
void twh_barot_(int *, int *, int *, double *);
void twh_bnrot_(int *, int *, int *, double *);
void twh_batraa_(int *, int *, int *, double *);
void twh_bntraa_(int *, int *, int *, double *);
void twh_batrac_(int *, int *, int *, double *);
void twh_bntrac_(int *, int *, int *, double *);
void twh_c_matrix_(int *, int *, int *, double *);
void twh_rmcomrot_(int *, int *, int *, double *);
void twh_rmcomtra_(int *, int *, int *, double *);
void twh_rmtraa_(int *, int *, int *, double *);
void twh_rmtrac_(int *, int *, int *, double *);
void twh_rmrot_(int *, int *, int *, double *);
void twh_arbcmofield_(int *, int *, double *);
void twh_comfield_(int *, int *, double *);
void twh_comtempfield_(int *, int *, double *);
void twh_coordfield_(int *, int *, double *);
void twh_coordstorage_(int *, int *, double *);
void twh_coordtemp_(int *, int *, double *);
void twh_cubeletweight_(int *, int *, double *);
void twh_gyration_(int *, int *, double *);
void twh_eam_rho_real_(int *, int *, double *);
void twh_eam_rho_temp_(int *, int *, double *);
void twh_rcmu_(int *, int *, double *);
void twh_rmvol_(int *, int *, double *);
void twh_tmmc_weight_(int *, int *, double *);
void twh_v_semigrand_(int *, int *, double *);
void twh_wrap_foreign_energy_(int *, int *, double *);
void twh_wrap_foreign_lambda_lj_(int *, int *, double *);
void twh_wrap_foreign_lambda_c_(int *, int *, double *);
void twh_acsvol_(int *, int *, int *);
void twh_ewald_kmax_(int *, int *, int *);
void twh_glist_(int *, int *, int *);
void twh_globalpos_(int *, int *, int *);
void twh_growfrom_(int *, int *, int *);
void twh_grownum_(int *, int *, int *);
void twh_growprev_(int *, int *, int *);
void twh_growvalidcount_(int *, int *, int *);
void twh_logical_exist_(int *, int *, int *);
void twh_logical_exsched_(int *, int *, int *);
void twh_logical_moveme_(int *, int *, int *);
void twh_logical_periodic_(int *, int *, int *);
void twh_moltyp_(int *, int *, int *);
void twh_nboxi_(int *, int *, int *);
void twh_parall_(int *, int *, int *, int *);
void twh_chainlist_(int *, int *, int *, int *, int *);
