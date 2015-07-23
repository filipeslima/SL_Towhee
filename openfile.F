#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2008 Marcus G. Martin                       *
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
      subroutine twh_openfile (lfailure,number,name,directory,io_status
     &     ,io_access)
c     ******************************************************************
c     * opens files with the specified directory.                      *
c     * if the directory name is 'default' then this routine fetches   *
c     * the global io_directory and uses that                          *
c     * number is the unit number assigned to this file                *
c     * name is the name of the file                                   *
c     * io_status must be one of:                                      *
c     *  'OLD' -- file must exist                                      *
c     *  'NEW' -- file must not exist, will be created                 *
c     *  'SCRATCH' -- temporary file will be created, then deleted at  *
c     *           end of execution                                     *
c     *  'UNKNOWN' -- if file exists, use it, otherwise, create new    *
c     * io_status is required so that the user gets a comprehensible   *
c     * error message                                                  *
c     * The default value for io_status, if none specified, is         *
c     * 'UNKNOWN'. Therefore, use 'UNKNOWN' to duplicate functionality *
c     * from previous versions.                                        *
c     * is_access must be one of:                                      *
c     *   'APPEND' -- file opened to append data instead of overwrite  *
c     *   'DEFAULT' -- file opened without specifying an access.  If   *
c     *        you are not sure what access you want then use this one *
c     * lfailure: set to true if there is a fatal error                *
c     *                                                                *
c     * originally written 10-04-2002 by M.G. Martin                   *
c     * last modified 06-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_STRING_LENGTH
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer number
      character*(*) name,io_status,directory,io_access
c     --- local variables
      character * (MAXDIRLENGTH+30) filename
      character*(MAXDIRLENGTH) localdir
      integer dirlength

      if (io_status .ne. 'OLD' .and.
     &    io_status .ne. 'NEW' .and.
     &    io_status .ne. 'SCRATCH' .and.
     &    io_status .ne. 'UNKNOWN') then
         write(6,*) 'Illegal io_status passed to twh_openfile'
         write(6,*) '  (name = ',name,')'
         lfailure = .true.
         return
      endif

      if ( directory .eq. 'DEFAULT' ) then
c        --- get the global io_directory
         call twh_io_directory(GLB_GET,localdir)
      else
c        --- use the directory passed into this subroutine
         localdir = directory
      endif

      dirlength = twh_get_string_length(localdir)

c     --- set the filename
      if ( dirlength .eq. 0 ) then
         filename = name
      else
         filename = localdir(1:dirlength)//name
      endif
      if ( io_access .eq. 'DEFAULT' ) then
#ifdef INTEL_RECL
         open(number,file=filename,form='formatted',status=io_status
     &        ,RECL=512,ERR=100)
#else
         open(number,file=filename,form='formatted',status=io_status
     &        ,ERR=100)
#endif
      elseif ( io_access .eq. 'APPEND' ) then
#ifdef INTEL_RECL
         open(number,file=filename,form='formatted',status=io_status
     &        ,access=io_access,RECL=512,ERR=100)
#else
         open(number,file=filename,form='formatted',status=io_status
     &        ,access=io_access,ERR=100)
#endif
      else
         write(6,*) 'Illegal io_access passed to twh_openfile'
         write(6,*) '  (name = ',name,')'
         lfailure = .true.
         return
      endif

      return

c     --- get here in case of error
100   if (io_status .eq. 'NEW') then
         write(6,*) 'File ',filename,' already exists.'
      else
         write(6,*) 'Error opening file'
         write(6,'(a)') filename
      endif
      lfailure = .true.
      return
      end
