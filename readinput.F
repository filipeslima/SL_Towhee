#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
c     * see the file license.gpl for the full license information      *
c     *                                                                *
c     * This program is free software; you can redistribute it and/or  *
c     * modify it under the terms of the GNU General Public License    *
c     * as published by the Free Software Foundation; either version 2 *
c     * of the License, or (at your option) any later version.         *
c     *                                                                *
c     * This program is distributed in the hope that it will be useful,*
c     * but WITHOUT ANY WARRANTY; without even the implied warranty of *
c     * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the  *
c     * GNU General Public License for more details.                   *
c     *                                                                *
c     * You should have received a copy of the GNU General Public      *
c     * License along with this program; if not, write to the Free     *
c     * Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,*
c     * MA  02111-1307, USA.                                           *
c     ******************************************************************
      subroutine twh_readinput(lfailure,lfinish)
c     ******************************************************************
c     * opens towhee_input and reads in the inputformat                *
c     * then it calls the appropriate input reading routine            *
c     *                                                                *
c     * rewritten 05-15-2002 by M.G. Martin                            *
c     * last modified 08-12-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#include "wrapperF.h"

c     --- variables passed to/from the subroutine
      logical lfailure,lfinish
c     --- local variables
      logical loutput
      character*15 inputformat,idname
      character*(MAXDIRLENGTH) io_dir

      idname = 'readinput'
      lfinish = .false.

c     --- retrieve constants
      call twh_loutput(GLB_GET,loutput)

c     --- open the input file for towhee
c     --- This is quite klunky right now: for parstyle = jobfarm,
c     --- io_directory is set; for other parstyles, 
c     --- warp_towhee_input_file may be set from the command line to
c     --- allow reading from different input files.  It
c     --- would be good to simplify this...
      call twh_openfile(lfinish,4,wrap_towhee_input_file,'DEFAULT'
     &     ,'OLD','DEFAULT')
      if ( lfinish ) return
      if ( loutput ) then
         write(6,'(1x,a,1x,a)') "Reading from towhee_input file:", 
     &        wrap_towhee_input_file
         call twh_io_directory(GLB_GET,io_dir)
         if (io_dir .eq. '') then
            write(6,'(1x,a)') "in directory: current directory"
         else 
            write(6,'(1x,a,1x,a)') "in directory:",io_dir
         endif
      endif

c     --- open towhee_altinp
      call twh_openfile(lfinish,15,'towhee_altinp','DEFAULT','UNKNOWN'
     &     ,'DEFAULT')
      if ( lfinish ) return

c     --- read in inputformat
      call twh_read_labeled_string(lfailure,4,'inputformat',loutput
     &     ,idname,inputformat)
      if ( lfailure ) return

c     --- head off to the various read routines
      if ( inputformat .eq. 'LAMMPS' ) then
         call twh_readlammps(lfinish)
         if ( lfinish ) return
         write(6,*) 'In order to run Towhee with the LAMMPS generated'
         write(6,*) 'input files use the following commands and then'
         write(6,*) 'run towhee again'
         write(6,*) 'cp towhee_altinp towhee_input'
         write(6,*) 'cat lammps_coords_* > towhee_coords'
         lfinish = .true.
      elseif ( inputformat .eq. 'Towhee' ) then
c        --- read in the rest of the file according to everything 
c        --- below
         call twh_readtowhee(lfinish)
      else
         write(6,*) 'READINPUT: invalid inputformat:',inputformat
         lfinish = .true.
      endif

c     --- close towhee_input and towhee_altinp
      close(4)
      close(15)

      return
      end
