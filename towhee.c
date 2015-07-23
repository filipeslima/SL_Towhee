/*     $Id: towhee.c,v 1.27 2010/09/04 14:16:50 lperi Exp $ */
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
/* rewritten from fortran 05-08-2003 by M.G. Martin 
   last modified 07-02-2009 by M.G. Martin
*/

/* globals  */
#include "globalc.h"
#include <time.h>
#include <unistd.h>
#ifdef INTEL_VISUAL_FORTRAN
#include "getopt.h"
#endif

#define USERUNTIME

/* fire up the code and pass some variables to the fortran portion */
int main(int argc, char* argv[])
{
  user_args args;
  int myproc = -1;

  initialize_args( &args );
#ifndef DISABLE_CLA
/* If command-line argument parsing disabled, only the defaults set in the
 * call above will be operative  */
  parse_args(argc, argv, &args);
#endif

  myproc = initialize_MPI(argc, argv);

/* initialize data structures which pass data between the C and Fortran parts
 * of towhee, and set the towhee_output_file in this structure. */
  initialize_wrapper_strings();
  set_towhee_input_file(args.towhee_input_fn);

  switch (args.parstyle) {
  case NONE:
    do_none(&args);
    break;
  case TRAMONTO:
    do_tramonto(&args, myproc);
    break;
  case JOBFARM:
    do_jobfarm(&args, myproc);
    break;
  case REX:
    do_rex(&args, myproc);
    break;
  case LCAO:
    do_lcao(&args, myproc);
    break;
  default:
    printf("Bad mode, uncaught by parse_args.\n");
    printf("args.parstyle %d \n",args.parstyle);
    exit(1);
  }
  exit(0);
}

/*
 * Return this process' rank.  This is safe to call if MPI is not
 * defined (will return -1 if MPI is not defined)
 */
int initialize_MPI(int argc, char* argv[]) {
#ifndef USEMPI
  return -1;
#else
  int myproc;
  /* initialize mpi */
  MPI_Init( &argc , &argv );
  MPI_Comm_rank(MPI_COMM_WORLD,&myproc);
  MPI_Barrier(MPI_COMM_WORLD);
  return myproc;
#endif
}


void initialize_clocks() {
    /* initialize all clocks */
    int clocknum = 0;
    while ( clocknum < MAXTIMERS ) {
      timekeeper(0,clocknum,0);
      clocknum++;
    }
    /* start timing for clock 0 */
    timekeeper(1,0,0);
}

void finalize_clocks() {
    /* timing */
    timekeeper(2,0,0);
    printf("Simulations complete \n");
    /* output timing */
    timekeeper(3,0,0);
#ifdef USETRAMONTO
    timekeeper(3,1,0);
#endif
}

void do_none(user_args *args)
{
  /* last modified 06-25-2009 by M.G. Martin */
    const int zero = 0;
    const int readwrapper = LOG_FALSE;
    const int pstyle = NONE;
    /* enter the fortran portion of the code
     * output location is zero, indicating no io redirection
     * we are not reading from the wrapper data */
    towheemainloop_(&zero,&(args->quiet_mode),&readwrapper,&pstyle);
}

void do_tramonto(user_args *args, int myproc) {
#ifdef USEMPI
    /* this is the functionality for using Tramonto in parallel
       but nothing else in the code is parallel
       only node 0 enters the towhee part of the code */
  if ( myproc == 0 ) initialize_clocks();
  tramonto_control(args->quiet_mode,myproc);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( myproc == 0 ) finalize_clocks();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  printf("Warning: MPI not supported, parstyle TRAMONTO disabled.\n");
  printf("Doing nothing.\n");
  return;
#endif
}

void do_jobfarm(user_args *args, int myproc) {
#ifdef USEMPI
    /* this uses a manager/worker paradigm to farm out a group of 
       serial jobs to a batch of processors.  Note that node 0 
       is the manager and does no work, except that the manager will 
       pitch in once there is only one job left to do */
  FILE *fp;
  if ( myproc == 0 ) {
    initialize_clocks();
    open_towhee_parallel("towhee_parallel", &fp);
  }

  jobfarm(fp, args->quiet_mode, myproc);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( myproc == 0 ) {
    finalize_clocks();
  /* close towhee_parallel file here since we opened it here */
    fclose(fp);
  }


  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  printf("Warning: MPI not supported, parstyle JOBFARM disabled.\n");
  printf("Doing nothing.\n");
  return;
#endif
}

void do_rex(user_args *args, int mpi_rank) {
#ifdef USEMPI
  int mpi_size,allocsize,flag;
  double dvalue;
  rex_params p;
 
  if (read_rex_params(args->rex_fn, &p) == 1) {
    printf("Error reading REX Parameters.\n");
    exit(1);
  }

  /* print_rex_params(&p); */
  /*  Test to make sure universe is large enough. */
  MPI_Comm_size(MPI_COMM_WORLD,&mpi_size);
  if ((mpi_size-1) < p.num_nodes) {
    printf("There are %d MPI nodes; num_nodes must be %d or less\n",
        mpi_size, mpi_size-1);
    printf("    (currently num_nodes = %d)\n",p.num_nodes);
    exit(1);
  }

  /* Likewise, make sure that our rank isn't more than num_nodes.  (This could
   * happen since mpirun specifies number of instances, which is independent
   * of num_modes) */
  if (mpi_rank > p.num_nodes) {
    printf("My rank (%d) exceeds number of nodes (%d).\n", mpi_rank,
        p.num_nodes);
    exit(1);
  }


  /* Initialize random number generator (used for creating seeds for towhee) */
  srandom (time(0) + mpi_rank);

  /* allocate memory for the rex routines */
  allocsize = p.num_nodes;
  if ( allocsize == 0 ) allocsize++;
  flag = GLB_ALLOC;
  dvalue = 0.0;
  twh_wrap_foreign_energy_(&flag,&allocsize,&dvalue);
  twh_wrap_foreign_lambda_lj_(&flag,&allocsize,&dvalue);
  twh_wrap_foreign_lambda_c_(&flag,&allocsize,&dvalue);

  /* server runs for rank=0, except for singleton mode, where num_nodes=0. */
  if (mpi_rank == 0 && p.num_nodes != 0) 
    run_rex_server(args, &p);
  else
    run_rex_client(args, &p, mpi_rank);
  /* free memory for the rex routines */
  flag = GLB_FREE;
  twh_wrap_foreign_energy_(&flag,&allocsize,&dvalue);
  twh_wrap_foreign_lambda_lj_(&flag,&allocsize,&dvalue);
  twh_wrap_foreign_lambda_c_(&flag,&allocsize,&dvalue);

  MPI_Finalize();

#else
  printf("Warning: MPI not supported, parstyle REX disabled.\n");
  printf("Doing nothing.\n");
  return;
#endif
}

void do_lcao(user_args *args, int myproc) {
#ifdef USEMPI
    /* this is the functionality for using LCAO in parallel
       but nothing else in the code is parallel
       only node 0 enters the towhee part of the code */
  if ( myproc == 0 ) initialize_clocks();
  lcao_control(args->quiet_mode,myproc);
  MPI_Barrier(MPI_COMM_WORLD);
  if ( myproc == 0 ) finalize_clocks();
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Finalize();
#else
  printf("Warning: MPI not supported, parstyle LCAO disabled.\n");
  printf("Doing nothing.\n");
  return;
#endif
}


/*
 * parses the command line arguments and sets various global variables.
 */
void initialize_args( user_args *args) {
  /* last modified 05-02-2006 by M.G. Martin */
  /* First, initialize the args to default values */
  args->quiet_mode = FALSE;
  /* the default parstyle is JOBFARM for MPI, NONE otherwise */
#ifdef USEMPI
  args->parstyle = JOBFARM;
#else
  args->parstyle = NONE;
#endif
  (void)strcpy(args->rex_fn, "towhee_rex");
  (void)strcpy(args->towhee_input_fn, "towhee_input");
  args->random_seed = -1;  /* this value implies that a random "random seed"
			                   * will be generated on each client invocation. */
}

void parse_args(int argc, char* argv[], user_args *args) {
  /* based on,
   * http://www.gnu.org/software/libc/manual/html_node/Example-of-Getopt.html */
  char errmsg[256];
  /*  char *cvalue = NULL; */
  /* last modified 02-15-2008 by M.G. Martin */
  int index;
  int c;
  int i = 0;

  while ((c = getopt (argc, argv, "hqs:f:p:")) != -1)
    switch (c) {
      case 'h': /* help flag */
        print_usage(NULL);
        break;
      case 'q': /* quiet flag */
        args->quiet_mode = 1;
        break;
      case 's': /* random seed flag */
        args->random_seed = (int) strtol(optarg, NULL, 10);
        /* check for illegal argument */
        if (args->random_seed == 0 && errno == EINVAL) {
            sprintf(errmsg, "Illegal random seed argument %s", optarg);
            print_usage(errmsg);
        }
        break;
      case 'f': /* filename flag */
        (void)strncpy(args->rex_fn, optarg, MAXDIRLENGTH);
        break;
      case 'p': /* parstyle flag */
        if (!strcmp(optarg, "none")) args->parstyle = NONE;
        else if (!strcmp(optarg, "tramonto")) args->parstyle = TRAMONTO;
        else if (!strcmp(optarg, "jobfarm")) args->parstyle = JOBFARM;
        else if (!strcmp(optarg, "rex")) args->parstyle = REX;
        else if (!strcmp(optarg, "lcao")) args->parstyle = LCAO;
        else {
            sprintf(errmsg, "Unknown parstyle %s", optarg);
            print_usage(errmsg);
        }
        break;
      case '?': /*  unknown flag */
        print_usage(" ");
        break;
      default: /* weirdness */
        sprintf(errmsg, "Abnormal error parsing arguments.\n");
        print_usage(errmsg);
    }

  c  = FALSE;
  for (index = optind; index < argc; index++) {
    if (i == 0) 
      strncpy(args->towhee_input_fn, argv[index], MAXDIRLENGTH);
/*  elseif (i == 1) do something with second arg, if you want, etc. */
    else {
       printf ("Non-option argument: %s\n", argv[index]);
       c = TRUE;
    }
    i++;
  }
  if (c) print_usage(" ");

  /* Debugging */
/*
  printf("Arguments: \n");
  printf("REX filename: %s\n", args->rex_fn);
  printf("quiet_mode: %d\n", args->quiet_mode);
  printf("parstyle: %d\n", args->parstyle);
  printf("random_seed: %d\n", args->random_seed);
*/
}

/*
 * Print usage.
 * If error != NULL, error message is printed, followed by usage (meant to
 * indicate incorrect invocation).
 * This function exits the program, with appropriate error code.
 */
void print_usage(const char* error) {
char version[] = PACKAGE_VERSION;
char msg[] = "\
******************************************************************\n\
* MCCCS - Towhee: A Monte Carlo molecular simulation program     *\n\
* Copyright (C) 2003-2005 Marcus G. Martin                       *\n\
* see the file license.gpl for the full license information      *\n\
******************************************************************\n\
Version: %s\n\
Usage:\n\
towhee [-h] [-q] [-s seed] [-p mode] [-f rexfn] [filename]\n\
Options:\n\
  -h -- print this message \n\
  -q -- quiet mode, supresses most output\n\
  -p (none | tramonto | jobfarm | rex | lcao) -- parallel style\n\
     none: single processor run, output and input in current directory.\n\
     The following modes require compilation with MPI support\n\
     tramonto: tramonto parallel style\n\
     jobfarm: parallel execution in MPI universe of multiple jobs.\n\
            the file 'towhee_parallel' must be in the current directory.\n\
     rex: Replica Exchange run.\n\
     lcao: LCAO parallel style.\n\
     default mode is none if no MPI, jobfarm if MPI\n\
  -f rexfn -- Replica Exchange parameter filename.  Defaults to towhee_rex\n\
  -s seed -- set random number seed to pass to main towhee code for REX \n\
       parstyle.  In general, the random seed is set in towhee_input, and\n\
       this parameter will NOT override this.  For REX, however, a random\n\
       seed is automatically generated for every client run.  Setting this\n\
       parameter overrides the random seed generation.  The default value\n\
       of -1 enables random seed generation.\n\
  filename -- specifies the filename of the towhee_input file.  If not\n\
       specified, the default 'towhee_input' is used.  filename is ignored\n\
       if the parallel style is jobfarm, since towhee_input is specified\n\
       in the 'towhee_parallel' file for that case.\n\
\n\
Full documentation of MCCCS Towhee may be found here:\n\
http://towhee.sourceforge.net/index.html\n\
";
    if (error != NULL) {
        printf("Error: %s\n", error);
    }
    printf(msg, version);

    /* exit value 0 is success, otherwise a failure.
     * See: http://www.tldp.org/LDP/abs/html/exit-status.html */
    if (error != NULL) 
        exit(1);
    else 
        exit(0);
}
