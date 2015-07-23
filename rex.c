/*     $Id: rex.c,v 1.14 2010/09/04 14:16:49 lperi Exp $ */
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
   Implementation of Replica Exchange.
   Note that this code is relatively immature.  Specifically, only
   the relatively specialized case of swapping between replicas
   of differing lambda_LJ and lambda_C values (Scaled Lennard
   Jones potential).  Implementing the more generally useful case
   of temperature swaps would take place here.

   Replica Exchange is currently maintained in Towhee by Matt
   Wyczalkowski (maw2@cec.wustl.edu)
 
   written 09-22-2005 by MAW
   last modified  12-03-2008 by M.G. Martin
*/

#include "globalc.h"
#include <dirent.h>

#include "wrapperC.h"
#include <fcntl.h>
#include <math.h>
#include <limits.h>
#include <time.h>

/* where data will be written. */
#define DATA_DIR "run_data"

/* homebrew keys for the two MPI messages we send */
#define LAMBDA_KEY 1
#define ENERGY_KEY 2

int copy_file(const char*, const char*);

void run_rex_server(user_args *args, rex_params *paramwrap) {
#ifdef USEMPI
  /* Notation: k and l always refer to node number (aka client rank)
   * m and n always refer to lambda key (index into paramwrap->lambda_X[m]) */
  int i, k, l, m, round, attempt;
  int verbose = !args->quiet_mode;
  MPI_Status status;
  double dtmp[MAX_FOREIGN_LAMBDA];
  float temp = 298;
  char msg[256];
  /* yuck.  Clean this up. */
  FILE *fp;
  char swapfn[] = "swap.out";

  /* Dynamically allocate U as a num_nodes x num_nodes 2D array.
   * U stores the foreign energies, as reported by the nodes (in kcal / mol)
   * U[k][m] holds the energy of configuration k (as found on node rank k)
   * evaluated at lambdas given by m (i.e. lambda_lj[m], lambda_c[m])
   * Note that k ranges 1..num_nodes and m ranges 0..num_nodes-1 */
  double **U = (double**) twh_allocate2dMatrix(sizeof(double), paramwrap->num_nodes+1,
      paramwrap->num_nodes);

  /* Dynamically allocate J as array size num_nodes.
   * J maps the node number to the lambda value.  That is, J[k] = m, giving
   * the lambda index (m) for node k.
   * Note that k  ranges 1..num_nodes, m 0..num_nodes-1. */
  int* J = (int*) twh_allocateVector(sizeof(int), paramwrap->num_nodes);

  /* Initialize J.  for the first round, each node k gets lambda index k-1. */
  for (k = 1; k <= paramwrap->num_nodes; k++) J[k] = k-1;

  /* We assume a temperature for now -- important for swaps. */
  printf("NOTE: Assuming that temperature = 298K\n");

  /* Now open the file that all swap info will be written to. */
  if (!(fp=fopen(swapfn,"w"))) {
    sprintf(msg, "Cannot open file %s\n", swapfn);
    perror(msg);
    exit(1);
  }
  

  /* Begin the rounds. */
  for (round = 1; round <= paramwrap->num_rounds; round++) {
    /* Make sure we're all in this together. */
    MPI_Barrier(MPI_COMM_WORLD);

    timestamp(msg);
    fprintf(fp, "%s: *******  Round %d starting *******\n", msg, round);
    if (verbose) {
      printf("%s: *******  Round %d starting (Server)  *******\n", msg, round);
    }

    /* First, write the file which tells us what client has what lambda in 
     * this round. */
    if (write_lambda_file(J, paramwrap, DATA_DIR, round, "lambda_key.txt") != 0) {
      printf("Error creating lambda_key file.\n");
      exit(1);
    }

    /* send lambda key to each node.  lambda key m is given by J[k] */
    for (k = 1; k <= paramwrap->num_nodes; k++) {
      m = J[k]; 
      MPI_Send(&m,1,MPI_INT,k,LAMBDA_KEY,MPI_COMM_WORLD);
    }

    /* Now we hang out and wait for responses from clients.
     * They will send back data for 'foreign energy', which we'll collect into
     * the data structure U.  We don't know the order in which they'll report
     * back. */
    for (k = 0; k < paramwrap->num_nodes; k++) {
      MPI_Recv(dtmp,paramwrap->num_nodes,MPI_DOUBLE,MPI_ANY_SOURCE,ENERGY_KEY,MPI_COMM_WORLD,
	       &status);
      timestamp(msg);
      printf("%s: Got response from %d\n", msg, status.MPI_SOURCE);
      /* Copy the energies to U */
      for (i = 0; i < paramwrap->num_nodes; i++) 
        U[status.MPI_SOURCE][i] = dtmp[i];
    }

    /* Once we have all the energies,  perform the swaps.  
     * don't bother swapping on the last round, though. */
    if (round == paramwrap->num_rounds) continue;

    /* Also, num_clients == 1, don't perform swaps, since cannot get distinct
     * nodes to swap against. */
    if (paramwrap->num_nodes == 1) continue;

    for (attempt = 0; attempt < paramwrap->num_swaps; attempt++) {
      /* Select two nodes at random.  If they are the same, redo the loop
       * (without incrementing loop counter).  Remember -- k and l start at 1. */
      k = (random() % paramwrap->num_nodes) + 1;
      l = (random() % paramwrap->num_nodes) + 1;
      if (k == l) {
        attempt--;
        continue;
      }

      fprintf(fp, "%d: Attempting swap between nodes %d (%6.3f, %6.3f)", attempt,
          k, paramwrap->lambda_lj[ J[k] ], paramwrap->lambda_c[ J[k] ]);
      fprintf(fp," and %d (%6.3f, %6.3f)\n",
          l, paramwrap->lambda_lj[ J[l] ], paramwrap->lambda_c[ J[l] ]);

      /* See if we are two swap the replicas, and do so if necessary.
       * To swap, we switch the lambda assignments between nodes k and l */
      int swap = attempt_exchange(temp, k, l, U, J, paramwrap, fp);
      if (swap) {
        m = J[k];
        J[k] = J[l];
        J[l] = m;
      }
    }
  }

  fclose(fp);
  /* free all variables mallocked here. */
  for (k = 0; k < paramwrap->num_nodes+1; k++) {
    free(U[k]);
  }
  free(U);
  free(J);
#endif
}


/*
 * attempt_exchange returns true or false depending on whether a replica
 * exchange should take place between two replicas identified by the task IDs K
 * and L.  This is the heart of the replica exchange algorithm.
 *
 * Given two task IDs, K and L, we first find the corresponding lambda indices
 * M, N, using M=J(K), N=J(L).  Then, we find du as,
 *  du = U(K,N) + U(L,M) - U(K,M) - U(L,N)
 *
 * A replica exchange will take place with the probability, min{1, exp(-du*B)}
 * That is, if du <= 0, exchange occurs; if du>0, exchange occurs with a
 * Boltzmann probability as given.  B = 1/R*T, where R is the gas constant and
 * T is the temperature in Kelvin.
 *
 * Please see documentation in multiswap-rex.pdf.
 *
 * For a reference, see Nymeyer et al., "Atomic Simulations of Protein Folding,
 * Using the Replica Exchange Algorithm" in Methods in Enzymology, Vol. 383
 * (2004)
 * 
 * k and l are indices of nodes we wish to try to swap
 * Neither U nor J are modified
 * returns 1 if swap should take place, 0 if it should not.
 * outputs to file fp if it is not null.
 */
int attempt_exchange(float temp, int k, int l, double **U, int *J, 
    rex_params *p, FILE *fp) {
  int swap;
  double prob=0., rnum=0.;
  /* get lambda indices */
  int m = J[k];
  int n = J[l];


  /* Get energies */
  double Ukn = U[k][n];
  double Ulm = U[l][m];
  double Ukm = U[k][m];
  double Uln = U[l][n];
  
  double du = Ukn + Ulm - Ukm - Uln;

  if (du <= 0.) swap = TRUE;
  else {
    /* Reference: http://scienceworld.wolfram.com/physics/UniversalGasConstant.html */
    const float R = 0.0019872; /* kcal / (mol K) */
    prob = exp(-du / (R * temp));
    /* rnum is 0..1, inclusive */
#ifdef INTEL_VISUAL_FORTRAN
    rnum = (double) RANDOM() / (double) LONG_MAX;
#else
    rnum = (double) random() / (double) LONG_MAX;
#endif
    swap = (rnum < prob);
  }

  if (fp != NULL) {
    fprintf(fp,"dU = %-12.3f %-12.3f - (%-12.3f %-12.3f) = %-12.3f\n", Ukn, Ulm, Ukm, Uln, du);
    if (du >= 0.) 
      fprintf(fp,"    prob = %f, rnum = %f => %s\n", prob, rnum,
          swap ? "SWAP" : "STAY");
    else
      fprintf(fp,"    SWAP\n");
  }

  return swap;
}


void run_rex_client(user_args *args, rex_params *paramwrap, int rank)
{
  /* last modified 07-02-2009 by M.G. Martin */
#ifdef USEMPI
  const int readwrapper = LOG_TRUE;
  const int zero = 0;
  const int pstyle = REX;
  int m, loop_num_nodes;

  char output_dir[MAXDIRLENGTH];
  int am_singleton;
  int round, lambda_key;
  int flag;
  double dvalue;
  double dtmp[MAX_FOREIGN_LAMBDA];

  MPI_Status status;

  am_singleton = (rank == 0) ? TRUE : FALSE;
  /* Adjust the number of nodes to loop through for singleton case. */
  loop_num_nodes = am_singleton ? 1 : paramwrap->num_nodes;

  /* Initialize the foreign lambda list which will be passed to the towhee
   * mainloop.  Since this does not change a each round, do this now. */
  flag = GLB_SET;
  for (m = 1; m <= loop_num_nodes; m++)  {
    dvalue = paramwrap->lambda_lj[m];
    twh_wrap_foreign_lambda_lj_(&flag,&m,&dvalue);
    dvalue = paramwrap->lambda_c[m];
    twh_wrap_foreign_lambda_c_(&flag,&m,&dvalue);
  }
  wrap_.lredirect_stdout = TRUE;
  wrap_.num_steps = paramwrap->num_steps;

  /* Begin the rounds. */
  for (round = 1; round <= paramwrap->num_rounds; round++) {
    /* Make sure we're all in this together. */
    if (!am_singleton) MPI_Barrier(MPI_COMM_WORLD);

    /*    if (!args->quiet_mode)
     *      printf("*******  Round %d starting (client %d) *******\n", round, rank); */

    /* Obtain lambda_key.  If singleton, we know what it is.  Otherwise, get
     * it from the server. */
    if (am_singleton) lambda_key = 0;
    else 
      MPI_Recv(&lambda_key,1,MPI_INT,0,LAMBDA_KEY,MPI_COMM_WORLD,&status);

    /* Set the data which will be passed to mainloop */
    wrap_.lambda_lj = paramwrap->lambda_lj[lambda_key];
    wrap_.lambda_c = paramwrap->lambda_c[lambda_key];
    wrap_.num_foreign_lambda = am_singleton ? 1 : paramwrap->num_nodes;

    /* Create the random seed to pass, or read it from args->random_seed.
     * if args->random_seed == -1, create a random random seed. */
    if (args->random_seed == -1)
      wrap_.random_seed = (int) random();
    else {
      printf("Using passed random seed %d\n", args->random_seed);
      wrap_.random_seed = args->random_seed;
    }

    if (round == 1) {
      /* On the first round, do things a bit differently.  linit for towhee will
       * be true, and we will be starting from a towhee_coords file. */
      if (populate_directory_initial(output_dir, DATA_DIR, rank) == 1) {
        printf("Error reported in populate_directory_initial\n");
        exit(1);  
      }
      wrap_.linit = TRUE;
    } else {
      /* subsequent runs restart. */
     if (populate_directory_restart(output_dir, DATA_DIR, round, rank) == 1) {
        printf("Error reported in populate_directory_restart\n");
        exit(1);  
      }
      wrap_.linit = FALSE;
    }

    set_output_dir(output_dir);

    /* Enter the fortran portion of the code, where all the work gets done. */
    towheemainloop_(&zero,&(args->quiet_mode),&readwrapper,&pstyle); 

    /* Now print the foreign lambdas. */
    if (!args->quiet_mode) {
      /* printf("Client rank %d finished round %d \n", rank, round);
       * printf("  lambda key: %d  (native lambdas: LJ = %5.3f C = %5.3f)\n", lambda_key, 
       * p->lambda_lj[lambda_key], paramwrap->lambda_c[lambda_key]);
       * printf("Foreign Energies:\n");
       * printf("  key   energy\n");
       * for (m = 0; m < loop_num_nodes; m++) 
       * flag = GLB_GET;
       * twh_wrap_foreign_energy_(&flag,&m,&dtmp[m]);
       * printf("  %3d = %15.3f\n", m, dtmp[m]);
       * printf("\n"); */
    }

    /* if singleton mode, we go on for only one round (we can't swap!); so
     * exit the while loop here. */
    if (am_singleton) break;

    /* retrieve the foreign energies */
    flag = GLB_GET;
    for (m = 1; m <= paramwrap->num_nodes; m++) {
      twh_wrap_foreign_energy_(&flag,&m,&dtmp[m-1]);
    }

    /* Now send the foreign energy back to the server */
    MPI_Send(dtmp,paramwrap->num_nodes,MPI_DOUBLE,0,ENERGY_KEY,MPI_COMM_WORLD);

  }
#endif
}


/*
 * Write the "lambda map" file in the 'round'-level directory.  This file
 * tells us what lambda key is associated with each client.
 * Returns 0 on success, 1 on failure.
 * last modified 03-06-2014 by M.G. Martin
 */
int write_lambda_file(const int* J, const rex_params *p, const char* base, int round, const char* fn) {
  FILE *fp;
  char fullname[MAXDIRLENGTH];
  char msg[256];
  int m, k;

  /* First, create the round-level directory and get its name.  Set rank=-1 to
   * do this.  Here, fullname is the name of round directory. */
  if (get_rex_directory(fullname, base, round, -1, TRUE) == 1) {
    printf("Error creating directory for round %d\n", round);
    return 1;
  }

  /* Append "/" and fn to fullname */
  (void)strncat(fullname, "/", MAXDIRLENGTH - strlen(fullname) - 1);
  (void)strncat(fullname, fn, MAXDIRLENGTH - strlen(fullname) - 1);

  /* Now open the file */
  if (!(fp=fopen(fullname,"w"))) {
    sprintf(msg, "Cannot open file %s\n", fullname);
    perror(msg);
    return 1;
  }
  
  fprintf(fp, "  Client   Lambda Key (lambda_lj   lambda_c)\n");
  /* We are looping over the nodes k, indexed 1..num_nodes
   * m is the corresponding lambda key */
  for (k = 1; k <= p->num_nodes; k++) {
    m = J[k];
    fprintf(fp, "   %3d       %3d     %8.3f    %8.3f\n", k , m, 
        p->lambda_lj[m], p->lambda_c[m]);
  }
  fclose(fp);
  return 0;
}





/*
 * Create and populate a directory for an initial run (that is, a run which
 * does not restart a previous run).
 * This consists of these steps:
 *  - create a run directory with a round of 1 and given rank and base
 *  - copy tohwee_coords from current directory to the run directory
 * 
 * Note that it is assumed that we are executing in the directory where 
 * 'towhee_coords' exists.
 *
 * Be sure to set wrap_.linit to TRUE.
 * 
 * Return value is the error code: 0 for success, 1 for failure.
 * last modified 03-07-2014 by M.G. Martin
 */
int populate_directory_initial(char* dirname, const char* base, int rank) {
    char output_file[MAXDIRLENGTH];
    char *coords_fn = "towhee_coords";

    if (get_rex_directory(dirname, base, 1, rank, TRUE) == 1) {
        printf("populate_directory_initial: Error reported in get_rex_directory\n");
        return 1;
    }

    /* create output filename */
    (void)strncpy(output_file,dirname, MAXDIRLENGTH);
    (void)strncat(output_file, "/", MAXDIRLENGTH - strlen(output_file) - 1);
    (void)strncat(output_file, coords_fn, MAXDIRLENGTH - strlen(output_file) - strlen(coords_fn));
    if (copy_file(coords_fn, output_file) == 1) {
        printf("populate_directory_initial: Error reported in copying \
file %s to %s\n", coords_fn, output_file);
        return 1;
    }
    return 0;
}

/*
 * Create and populate a directory for a restart run (that is, a run which
 * restarts a previous run).
 * This consists of these steps:
 *  - create a run directory with the given round, rank and base
 *  - copy towhee_final from the previous run directory (the directory with
 *    the same base and rank, but round-1) to towhee_initial in the run
 *    directory.
 * 
 * Be sure to set wrap_.linit to FALSE.
 * 
 * Return value is the error code: 0 for success, 1 for failure.
 * last modified 03-07-2014 by M.G. Martin
 */
int populate_directory_restart(char* dirname, const char* base, int round, int rank) {
    char old_dirname[MAXDIRLENGTH];
    char output_file[MAXDIRLENGTH], old_output_file[MAXDIRLENGTH];
    char *in_fn = "towhee_final";
    char *out_fn = "towhee_initial";

    /* not creating anything here, just getting the name */
    if (get_rex_directory(old_dirname, base, round-1, rank, FALSE) == 1) {
      /* This is really weird (impossible, in present implementation) */
        printf("populate_directory_restart: Error reported in get_rex_directory\n");
        printf("(just getting old directory name.  Really weird.)\n");
        return 1;
    }

    if (get_rex_directory(dirname, base, round, rank, TRUE) == 1) {
        printf("populate_directory_restart: Error reported in get_rex_directory\n");
        return 1;
    }

    /* create input filename */
    (void)strncpy(old_output_file,old_dirname, MAXDIRLENGTH);
    (void)strncat(old_output_file, "/", MAXDIRLENGTH - strlen(old_output_file) - 1);
    (void)strncat(old_output_file, in_fn, MAXDIRLENGTH - strlen(old_output_file) - strlen(in_fn));

    /* create output filename */
    (void)strncpy(output_file,dirname, MAXDIRLENGTH);
    (void)strncat(output_file, "/", MAXDIRLENGTH - strlen(output_file) - 1);
    (void)strncat(output_file, out_fn, MAXDIRLENGTH - strlen(output_file) - strlen(out_fn));
    if (copy_file(old_output_file, output_file) == 1) {
        printf("populate_directory_initial: Error reported in copying \
file %s to %s\n", old_output_file, output_file);
        return 1;
    }
    return 0;
}

/*
 * Obtain the directory name where output of this round for this
 * rank will go, and optionally create the directory, too (as specified by
 * create_dir).
 * It is assumed that the directory specified by base exists; error if it
 * does not.
 * The directory at the round directory will be created if necessary; it is
 * assumed that the directory at the rank level does not exist, and it is an
 * error if it does.
 * As a special condition, if rank=-1, the rank directory will *not* be
 * created, and the dirname of the round directory will be returned (and
 * created, as specified by create_dir).
 * Directory name will be returned in dirname, and will be no longer than
 * MAXDIRLENGTH (it is assumed that dirname is that size or larger).
 * Currently, the dirname format is base/round/rank
 * Specifically, there is no terminating forward-slash
 * Returns 0 on success, 1 on failure.
 */
int get_rex_directory(char* dirname, const char* base, int round, int rank,
        int create_dir) {
    DIR *dd=NULL;
    char round_dn[MAXDIRLENGTH];

    /* setting rank to -1 is a special condition: we care only about the
     * round-level directory. */
    if (rank == -1) 
      (void) snprintf(dirname, MAXDIRLENGTH, "%s/N%02d", base, round);        
    else
      (void) snprintf(dirname, MAXDIRLENGTH, "%s/N%02d/K%02d", base, round, rank);        
    if (! create_dir) return 0;

    /* Make sure the base dir exists. */
    if ((dd = opendir(base)) == NULL) {
        printf("Error opening base directory %s: %s\n", base, strerror(errno));
        return 1;
    }
    closedir(dd);

    /* get the dir name at the round level, and see if it exists. */
    (void) snprintf(round_dn, MAXDIRLENGTH, "%s/N%02d", base, round);        
    if ((dd = opendir(round_dn)) == NULL) {
      /* directory does not exist.  Create it. */
        if (mkdir(round_dn, 0755) == -1) {
            printf("Error creating directory %s: %s\n", round_dn,
                strerror(errno));
            return 1;
        }
    } else closedir(dd);

    /* If rank == -1, we are done. */
    if (rank == -1) return 0;

    /* Finally, create bottom-level (rank) dir.  Check if it exists -- should
     * not. */
    if ((dd = opendir(dirname)) != NULL) {
        closedir(dd);
        printf("Error: directory %s exists.\n", dirname);
        return 1;
    } else {
        if (mkdir(dirname, 0755) == -1) {
            printf("Error creating directory %s: %s\n", dirname,
                strerror(errno));
            return 1;
        }
    }
    return 0;
}

/*
 * Due to the differnet ways C and Fortran handle strings, its good hygiene to
 * initialize all of the strings which will be shared between the two with
 * spaces.
 */
void initialize_wrapper_strings() {
    int i;
    for (i = 0; i < MAXDIRLENGTH; i++) {
        swrap_.towhee_input_file[i] = ' ';
        swrap_.output_dir[i] = ' ';
    }
}

/*
 * The following are convenience functions to set the wrapper string
 * variables.
 */
void set_towhee_input_file(const char* text) {
  /* Strings passed to fortran via wrapper must be copied to preserve the
   * whitespace padding, otherwise they print out weird.
   */
   int l = strlen(text) < MAXDIRLENGTH ? strlen(text) : MAXDIRLENGTH;
   (void)strncpy(swrap_.towhee_input_file, text, l);
}

void set_output_dir(const char* text) {
   int l = strlen(text) < MAXDIRLENGTH ? strlen(text) : MAXDIRLENGTH;
   (void)strncpy(swrap_.output_dir, text, l);
}

/*
 * Copy a file from src to dst.
 * Code based on: http://www.metalshell.com/view/source/100/
 * Returns 0 on success, 1 on failure.
 */
int copy_file(const char* src, const char* dst) {
  FILE *srcF, *dstF;
  char msg[1024];
  int c;

  if((srcF = fopen(src, "r")) == NULL) {
    sprintf(msg, "copy_file: Error opening source file %s", src);
    perror(msg);
    return 1;
  }

  if((dstF = fopen(dst, "w+")) == NULL) {
    sprintf(msg, "copy_file: Error opening destination file %s", dst);
    perror(msg);
    return 1;
  }

  while ((c = fgetc(srcF)) != EOF) 
    fputc(c, dstF);

  fclose(srcF);
  fclose(dstF);
  return 0;
}

/*
 * Get next non-comment data line.  Please make line buffer of size 128.
 * Successfully deals with blank lines (ignored), lines beginning with #
 * (ignored)
 * returns 0 on success, 1 on failure.
 */
int read_data_line(FILE *fp, char *line) {
  char tmp[128];

  while (TRUE) {
    /* read in a line */
    if (fgets(line, 128, fp) == NULL) {
      printf("Unexpected End of File.\n");
      return 1;
    }
    /* see if it is a comment */
    if (sscanf(line, " #%s", tmp) != 0) continue;
    /* non-comment found */
    return 0;
  }
}

/*
 * Checks to see that line contains expected label, and writes error message
 * if it does not.  Returns 1 if match fails, 0 if it succeeds.
 * fn is the filename.
 */
int check_label(const char *line, const char *label, const char *fn) {
  if (!strcmp(line, label)) {
    printf("Error in file %s.\n", fn);
    printf("Expected label %s; found %s\n", label, line);
    return 1;
  }
  return 0;
}

/*
 * read rex parameters from given filename, fill in p.
 * returns 0 on success, 1 on failure.
 *
 * the parameter file is of a format similar to towhee_input and such;
 * consists of a label on one line, followed by one or more lines with data.
 * The labels are in a proscribed sequence.  Lines beginning with the # symbol
 * are comments and are ignored.
 */
int read_rex_params(const char* rex_fn, rex_params *p) {
    FILE *fp;
    char msg[256], line[128];
    int n, num_lam;

    if (!(fp=fopen(rex_fn,"r"))) {
      sprintf(msg, "Cannot open file %s\n", rex_fn);
      perror(msg);
      return 1;
    }

    if (read_data_line(fp, line)) return 1;
    if (check_label(line, "num_rounds", rex_fn)) return 1;
    if (read_data_line(fp, line)) return 1;
    if (!sscanf(line, " %d", &(p->num_rounds))) {
      printf("Error reading num_rounds\n");
      return 1;
    }

    if (read_data_line(fp, line)) return 1;
    if (check_label(line, "num_steps", rex_fn)) return 1;
    if (read_data_line(fp, line)) return 1;
    if (!sscanf(line, " %d", &(p->num_steps))) {
      printf("Error reading num_steps\n");
      return 1;
    }

    if (read_data_line(fp, line)) return 1;
    if (check_label(line, "num_swaps", rex_fn)) return 1;
    if (read_data_line(fp, line)) return 1;
    if (!sscanf(line, " %d", &(p->num_swaps))) {
      printf("Error reading num_swaps\n");
      return 1;
    }

    if (read_data_line(fp, line)) return 1;
    if (check_label(line, "num_nodes", rex_fn)) return 1;
    if (read_data_line(fp, line)) return 1;
    if (!sscanf(line, " %d", &(p->num_nodes))) {
      printf("Error reading num_rounds\n");
      return 1;
    }

    /* now read in num_nodes many lambda pairs.
     * Note that if num_nodes = 0 (singleton mode), read in one pair of lambda's. */
    num_lam = p->num_nodes ? p->num_nodes : 1;
    for (n = 0; n < num_lam; n++) {
      if (read_data_line(fp, line)) return 1;
      if (sscanf(line, " %f %f", 
            &(p->lambda_lj[n]), &(p->lambda_c[n])) != 2) {
        printf("Error reading entry number %d of lambda values\n", n);
        return 1;
      }
    }
    fclose(fp);
    return 0;
}

void print_rex_params(rex_params *p) {
  int n, num_lam;
  printf("Replica Exchange Parameters:\n");
  printf("num_rounds = %d\n", p->num_rounds);
  printf("num_steps = %d\n", p->num_steps);
  printf("num_swaps = %d\n", p->num_swaps);
  printf("num_nodes = %d\n", p->num_nodes);
  if (p->num_nodes == 0)
    printf("   Singleton Mode -- No Server\n");

  printf("Lambda LJ    Lambda C\n");
  /* now read in num_nodes many lambda pairs, requiring one pair for the
   * special singleton case (num_nodes = 0). */
  num_lam = p->num_nodes ? p->num_nodes : 1;
  for (n = 0; n < num_lam; n++) 
    printf("%8.3f    %8.3f\n", p->lambda_lj[n], p->lambda_c[n]);
}

/*
 * Utility function for to get timestamp, with annoying newline stripped.
 */
void timestamp(char* s) {
  time_t tm = time(NULL);
  char *t = ctime(&tm);
  t[24] = ' ';
  sprintf(s,"%s", t);
}

