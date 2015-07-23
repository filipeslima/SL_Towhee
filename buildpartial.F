#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2011 Marcus G. Martin                       *
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
      subroutine twh_buildpartial(lfailure,imolty,testcord,ifile)
c     ******************************************************************
c     * reads in atoms from towhee_partial, determines which match the *
c     * expected atoms, and then prepares the data structures to build *
c     * the rest                                                       *
c     *                                                                *
c     * originally written 05-18-2005 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- debug toggle
#define DEBUG_BUILDPARTIAL 0
c     ---- variables passed to/from the subroutine
      logical lfailure
      integer imolty,testcord,ifile
c     --- local variables
      logical lreadnew
      integer count,iunit,icord,ichain,nunit
      character*3 key_shortname,aminoshort
      character*4 key_pdbname,pdbname
      
      double precision positions
      dimension positions(3)

#if DEBUG_BUILDPARTIAL
      write(6,*) 'start BUILDPARTIAL'
#endif
      count = 0
      ichain = 0
      lreadnew = .true.
c     --- read in the first line
      call twh_nunit(GLB_GET,imolty,nunit)
      do iunit = 1,nunit
         call twh_aminoshort(GLB_GET,imolty,iunit,aminoshort)
         call twh_pdbname(GLB_GET,imolty,iunit,pdbname)
#if DEBUG_BUILDPARTIAL
         write(6,*) 'iunit,pdbname:',pdbname
         write(6,*) 'aminoshort:',aminoshort
#endif
         if ( lreadnew ) then
            read(ifile,*) (positions(icord),icord=1,3),key_shortname
     &           ,key_pdbname
            lreadnew = .false.
#if DEBUG_BUILDPARTIAL
            write(6,*) 'key_shortname:',key_shortname
            write(6,*) 'key_pdbname:',key_pdbname
#endif
         endif
c        --- test against the key_shortname and key_pdbname
         if ( key_pdbname .eq. pdbname .and.
     &        key_shortname .eq. aminoshort ) then
c           --- we have a match, increment count
            count = count + 1
            call twh_coordinates(GLB_SET,testcord,ichain,iunit
     &           ,positions(1),positions(2),positions(3))
            call twh_cbmc_atomname_list(GLB_SET,imolty,count,iunit)
c           --- will need to read in another next loop
            lreadnew = .true.
         endif
      enddo
      call twh_cbmc_atomname_tot(GLB_SET,imolty,count)
c     --- make sure we had at least 1 atom
      if ( count .eq. 0 ) then
         write(6,*) 'BUILDPARTIAL: no matches found'
         write(6,*) 'Molecule type:',imolty
         lfailure = .true.
         return
      endif

#if DEBUG_BUILDPARTIAL
      write(6,*) 'finished BUILDPARTIAL'
#endif
      return
      end
