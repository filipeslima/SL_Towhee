/*     $Id: control.c,v 1.14 2009/06/26 17:23:35 marcus_martin Exp $ */
/* 
* MCCCS - Towhee: A Monte Carlo molecular simulation program
* Copyright (C) 2003-2006 Marcus G. Martin
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
/* this file contains several routines that control entry into 
   towhee when calling external C routines
   originally written to work with the Tramonto code as a library
   file last modified 11-08-2006 by M.G. Martin
*/
#include "globalc.h"

void tramonto_control(int outputmode, int myproc)
/*
  tramonto_control is called by towhee.c and is in charge of assigning 
  tasks to processors for a tramonto calculation.  the master enters 
  the mainloop while the slaves wait in tramontoloop
  originally writtin 08-08-2003 by M.G. Martin
  last modified 06-25-2009 by M.G. Martin 
*/
{
#ifdef USETRAMONTO
  int jobcount,logicflag;
  const int zero = LOG_FALSE;
  const int pstyle = TRAMONTO;
  if ( myproc == 0 ) { 
    /* set jobcount to 0 to signify an output to the screen */
    jobcount = 0;
    /* master enters towhee, does not return until everything is finished */
    towheemainloop_(&jobcount,&outputmode,&zero,&pstyle); 
    /* let the other nodes know we are done */
    logicflag = 0;
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&logicflag,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
  }
  else { 
    /* worker node enters the tramontoloop and remains there until the master 
       tells him to stop working */
    tramontoloop(); 
  }
#endif
}
void lcao_control(int outputmode, int myproc)
/*
  lcao_control is called by towhee.c and is in charge of assigning 
  tasks to processors for an LCAO calculation.  the master enters 
  the mainloop while the slaves wait in lcaoloop
  originally writtin 02-15-2006 by M.G. Martin
  last modified 06-25-2009 by M.G. Martin 
*/
{
  int jobcount;
  const int zero = LOG_FALSE;
  const int pstyle = LCAO;
  if ( myproc == 0 ) { 
    /* set jobcount to 0 to signify an output to the screen */
    jobcount = 0;
    /* master enters towhee, does not return until everything is finished */
    towheemainloop_(&jobcount,&outputmode,&zero,&pstyle); 
    /* let the other nodes know we are done */
    twhendquantum_();
  }
  else { 
    /* worker node enters the quantumloop and remains there until the master 
       tells him to stop working */
    twhquantumloop_(); 
  }
}

double tramonto_()
/* 
   tramonto_ is called by the master in the main fortran portion of 
   the code.  It then alerts all of the other nodes (that are waiting
   patiently in the tramontoloop) that it is time to call dftmain and
   get to work 
   originally written 06-20-2002 by M.G. Martin
   last modified 11-04-2004 by M.G. Martin 
*/
{
  double energy;
#ifdef USETRAMONTO
  int logicflag;
  /* update start of timing for timer 1 */
  timekeeper(1,1,0);
  logicflag=1;
  MPI_Barrier(MPI_COMM_WORLD);
  MPI_Bcast(&logicflag,1,MPI_INT,0,MPI_COMM_WORLD);
  MPI_Barrier(MPI_COMM_WORLD);
  dftmain(&energy);
  /* update end of timing for timer 1 */
  timekeeper(2,1,0);
#else
  energy = 0.0;
#endif
  return energy;
}

void tramontoloop()
/*
  tramontoloop is where the slave nodes wait for the command to enter 
  the external library dftmain and perform a parallel calculation
  originally written 06-20-2003 by M.G. Martin
  last modified 11-04-2004 by M.G. Martin
*/
{
#ifdef USETRAMONTO
  int logicflag;
  double dumb;
  logicflag = 1;
  while ( logicflag==1) {
    MPI_Barrier(MPI_COMM_WORLD);
    MPI_Bcast(&logicflag,1,MPI_INT,0,MPI_COMM_WORLD);
    MPI_Barrier(MPI_COMM_WORLD);
    if ( logicflag==1) dftmain(&dumb);
  }
  return;
#endif
}
    
void clearline(FILE * fp)
/* 
   clearline is a simple utility to read to a newline or end of file and then return
   originally written 06-20-2003 by M.G. Martin
   last modified 06-20-2003 by M.G. Martin
*/
{   
  int c;
  while ( ((c=getc(fp)) != EOF) && (c !='\n' ) &&  (c !='\r') );
}

void timekeeper( int command, int clocknum, int jobflag)
/*
  timekeeper is a utility to keep track of all the timing information
  it is still not as useful as I would like and some of the other 
  routines access the time_array as well
  originally written some time after 06-2003 by M.G. Martin
  last modified 04-25-2008 by M.G. Martin
*/
{
#ifdef USEMPI
  int days,hours,minutes;
  double dvalue;
  if ( command == 0 ) {
    /* initialize the clock */
    dvalue = 0.0;
    twh_time_array(GLB_SET,clocknum,&dvalue);
  }
  else if ( command == 1 ) {
    /* update the timer for the start of a timing by subtracting 
       the current time */
    dvalue =  - MPI_Wtime();
    twh_time_array(GLB_INCR,clocknum,&dvalue);
  }
  else if ( command == 2 ) {
    /* update the timer for the end of a timing by adding the current
       time */
    dvalue = MPI_Wtime();
    twh_time_array(GLB_INCR,clocknum,&dvalue);
  }
  else if (command == 3 ) {
    /* output the timing information for this clocknum */
    if ( clocknum == 0 ) {
      printf("Total Time: ");
    }
    else if ( clocknum == 1 ) {
      printf("Total time spent in Tramonto was ");
    }
    else if ( clocknum == 2 ) {
      printf("Time for Job # %d was ",jobflag);
    }
    /* output the time information at the end of the line */
    twh_time_array(GLB_GET,clocknum,&dvalue);
    /* convert into a days, hours, minutes, seconds format */
    if ( dvalue > 86400.0 ) {
      /* nonzero days */
      days = dvalue / 86400.0;
      printf("%d days, ",days);
      dvalue = dvalue - days*86400.0;
    }
    if ( dvalue > 3600.0 ) {
      /* nonzero hours */
      hours = dvalue / 3600.0;
      printf("%d hours, ",hours);
      dvalue = dvalue - hours*3600.0;
    }
    if ( dvalue > 60.0 ) {
      /* nonzero minutes */
      minutes = dvalue / 60.0;
      printf("%d minutes, ",minutes);
      dvalue = dvalue - minutes*60.0;
    }
    printf("%f seconds.\n",dvalue);
  }
#endif
}
