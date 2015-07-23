/*  $Id: jobfarm.c,v 1.12 2009/06/26 17:23:35 marcus_martin Exp $ */
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
  jobfarm runs the master and slave sides of the jobfarm
  this parallel function is used to run a list of jobs on one 
  or more processors.
  originally split out of towhee.c 08-08-2003 by M.G. Martin
  last modified 07-26-2013 by M.G. Martin
*/
#include "preproc.h"
#include "globalc.h"

void jobfarm(FILE * tp_fp, int outputmode, int myproc)
{
#ifdef USEMPI
  /* file pointers
     tp_fp: towhee_parallel
     fd: directoryfile  
     fff: force field file 
  */
  FILE *fd=NULL, *fff=NULL;
  /* MPI data types */
  MPI_Status status;
  int isendbuffer[2],irecvbuffer[2];
  int logicflag,numproc,workerrank,jobcount,oldjobcount,numjob;
  int dirlength,ffnumber,iii;
  double dsendbuffer[2],drecvbuffer[2];
  double dvalue;
  char directory[MAXDIRLENGTH], currentline[MAXDIRLENGTH];
  const int zero = LOG_FALSE;
  const int pstyle = JOBFARM;

  /* determine the number of processors */
  MPI_Comm_size(MPI_COMM_WORLD,&numproc);

  /* read in the number of jobs */
  if ( myproc == 0 ) {
    clearline(tp_fp);
    fscanf(tp_fp,"%d",&numjob);
    clearline(tp_fp);
    printf("Number of jobs: %d \n",numjob);
    /* move to start of directory list */
    clearline(tp_fp);
    clearline(tp_fp);
    clearline(tp_fp);
  }
  if ( myproc == 0 ) {
    /* Managers Section */
    jobcount = 0;
    while ( jobcount < numjob ) {
      jobcount++;
      logicflag=0;
      printf("Job # %d ",jobcount);
      /* initialize directory and get next input directory */
      directory[0] = '\0';
      fscanf(tp_fp,"%s",directory);   
      printf("%s \n",directory);
      if (feof(tp_fp)) {
        printf ("Insufficient entries in towhee_parallel file\n");
        break;
      }
      
      /* make sure directory ends with a '/' */
      dirlength = strlen(directory);
      if (directory[dirlength - 1] != '/')
      {
        if (dirlength + 1 > MAXDIRLENGTH)
        {
          printf("Maxiumum directory length exceeded for job # %d\n", 
              jobcount);
          break;
        }
        strcat(directory,"/");
	dirlength++;
      }
      
      /* append input file name to directory */
      if ( dirlength + 12 > MAXDIRLENGTH )
        {
          printf("Maxiumum directory length exceeded for job # %d\n", 
		 jobcount);
	  printf("Need to increase maxdirlength to at least %d \n"
		 ,dirlength+12);
          break;
        }
      strcat(directory,"towhee_input");
      
      /* make sure candidate input file can be opened */
      if (!(fd = fopen(directory,"r")))
      {
        printf ("towhee_input cannot be opened for job # %d\n", jobcount);
        continue;
      }
      
      /* look for forcefield files in input, and test */
      while (fgets(currentline,MAXDIRLENGTH,fd))
      {
        if (strstr(currentline,"ffnumber"))
        {
          /* read number of forcefield files using fgets 
             and convert to integer.
             (fscanf appears incompatible with fgets) */
          fgets(currentline,MAXDIRLENGTH,fd);
          ffnumber=atoi(currentline);
          /* move to first forcefield line */
          fgets(currentline,MAXDIRLENGTH,fd);
          printf("Job # %d has %d forcefield files\n",jobcount,ffnumber);
          for (iii = 0; iii < ffnumber; iii++)
          {
            fgets(currentline,MAXDIRLENGTH,fd);
            /* strip trailing '\n' from currentline */
            if (currentline[strlen(currentline) - 1] == '\n')
            {
              currentline[strlen(currentline) - 1] = '\0';
            }
            /* if file does not open, set flag and continue */
            if (!(fff=fopen(currentline,"r")))
            {
              printf("ForceField file %s cannot be opened for job # %d\n",
                  currentline,jobcount);
              logicflag = 1;
              continue;
            }
            fclose(fff);
          }
          /* continue since force field files found */
          break;
        }   
      }
      fclose(fd);
      if (logicflag) continue;   
      if ( numproc != 1 && jobcount != numjob ) {
        /* find the next available node */
        MPI_Recv(irecvbuffer,2,MPI_INT,MPI_ANY_SOURCE,0,MPI_COMM_WORLD
            ,&status);
        workerrank = irecvbuffer[0];
        oldjobcount = irecvbuffer[1];
        /* send this node the next jobcount */
        printf("Job # %d submitted to node %d \n",jobcount,workerrank);
        MPI_Send(&jobcount,1,MPI_INT,workerrank,1,MPI_COMM_WORLD);
        /* output the node that was given this task */
        if ( oldjobcount != 0 ) {
          /* we need to get timing information from the previous run of this node */
          MPI_Recv(drecvbuffer,2,MPI_DOUBLE,workerrank,2,MPI_COMM_WORLD
              ,&status);
	  dvalue = drecvbuffer[0];
          twh_time_array(GLB_SET,2,&dvalue);
          /* output the timing information */
          timekeeper(3,2,oldjobcount);
        }
      }
      else {
        /* last job or I'm the only worker - time to do it myself */
        /* initialize this clock */
        timekeeper(0,2,jobcount);
	/* start the clock */
        timekeeper(1,2,jobcount);
        printf("Job # %d submitted to node %d \n",jobcount,myproc);
        towheemainloop_(&jobcount,&outputmode,&zero,&pstyle);
        timekeeper(2,2,jobcount);
        /* output timing information */
        timekeeper(3,2,jobcount);
      }
    }
    /* master lets all the nodes know that we are done */
    jobcount = 0;
    workerrank = 1;
    while ( workerrank < numproc ) {
      MPI_Recv(irecvbuffer,2,MPI_INT,workerrank,0,MPI_COMM_WORLD
          ,&status);
      MPI_Send(&jobcount,1,MPI_INT,workerrank,1,MPI_COMM_WORLD);
      if ( irecvbuffer[1] != 0 ) {
        /* we need to get timing information from the previous run of this node */
        MPI_Recv(drecvbuffer,2,MPI_DOUBLE,workerrank,2,MPI_COMM_WORLD
            ,&status);
	dvalue = drecvbuffer[0];
        twh_time_array(GLB_SET,2,&dvalue);
        /* output the timing information */
        oldjobcount = irecvbuffer[1];
        timekeeper(3,2,oldjobcount);
      }
      workerrank++;
    }
  }
  else {
    jobcount = -1;
    isendbuffer[0] = myproc;
    isendbuffer[1] = 0;
    while ( jobcount != 0 ) {
      /* contact the manager for instructions */
      MPI_Send(isendbuffer,2,MPI_INT,0,0,MPI_COMM_WORLD);
      MPI_Recv(&jobcount,1,MPI_INT,0,1,MPI_COMM_WORLD,&status);
      if ( isendbuffer[1] != 0 ) {
        /* we need to send in the timing from our last job */
        MPI_Send(dsendbuffer,2,MPI_DOUBLE,0,2,MPI_COMM_WORLD);
        /* our message is sent, reset isendbuffer */
        isendbuffer[1] = 0;
      }
      if ( jobcount != 0 ) {
        /* initialize this clock */
        timekeeper(0,2,jobcount);
        /* start the clock */
        timekeeper(1,2,jobcount);
        towheemainloop_(&jobcount,&outputmode,&zero,&pstyle);
        /* end the clock */
        timekeeper(2,2,jobcount);
        /* store the elapsed time in the dsendbuffer */
	twh_time_array(GLB_GET,2,&dvalue);
        dsendbuffer[0] = dvalue;
        /* set flag so we alert the master that we have a timing to send */
        isendbuffer[1] = jobcount;
      }
    }
  }
#endif
}

/*
 * Open the towhee_parallel file of given name for reading, and return the
 * file pointer.  This file remains open, should be closed (and
 * possibly read) later.
 * Exits in case of error.
 */
void open_towhee_parallel(const char* fn, FILE **fp) {
    char msg[256];
    if (!(*fp=fopen(fn,"r"))) {
      sprintf(msg, "Cannot open file %s\n", fn);
      perror(msg);
      exit(1);
    }
    return;
}

