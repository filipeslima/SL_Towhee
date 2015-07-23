#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * see the file license.gpl for the full license information      *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_setinclude()
c     ******************************************************************
c     * this subroutine sets up the inclusion table which is used to   *
c     * determine whether two atoms on the same molecule have          *
c     * non-bonded interactions with each other                        *
c     * note that 1-4 coulombic and van der waals interactions are     *
c     * included in the torsional potential and are therefore not      *
c     * considered to be non-bonded interactions for the purposes of   *
c     * setting up this table                                          *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modifed 08-09-2011 by M.G. Martin                         *
c     ******************************************************************
      implicit none
#define FUNCTION_GET_IJBEN3
#define FUNCTION_GET_IJTOR4
#define FUNCTION_GET_INBEN
#define FUNCTION_GET_INTOR
#define FUNCTION_GET_INVIB
#include "functions.h"

c     --- local variables
      logical lexclude
      dimension lexclude(NUMAX)
      integer m,n,nb,imolty,count,nmolty,nunit

c     --- loop over all types of molecules
      call twh_nmolty(GLB_GET,nmolty)
      do imolty = 1, nmolty
         call twh_nunit(GLB_GET,imolty,nunit)
c        --- loop over all units 
         do m = 1, nunit
c           --- initiallize the array
            do n = 1,nunit
               lexclude(n) = .false.
            enddo
c           --- exclude self
            lexclude(m) = .true.
c           --- exclude all directly bonded beads (vibrations) 
            do n = 1, twh_get_invib(imolty,m)
               call twh_ijvib(GLB_GET,imolty,m,n,nb)
               lexclude(nb) = .true.
            enddo
c           --- exclude all next-nearest neighbor bonded beads (bending)
            do n = 1, twh_get_inben(imolty,m)
               nb = twh_get_ijben3(imolty,m,n)
               lexclude(nb) = .true.
            enddo
c           --- exclude all 3rd-nearest neighbor bonded beads(torsions) 
            do n = 1, twh_get_intor(imolty,m)
               nb = twh_get_ijtor4(imolty,m,n)
               lexclude(nb) = .true.
            enddo

c           --- determine the number of excluded atoms and build arrays
            count = 0
            do n = 1,nunit
               if ( lexclude(n) ) then
                  count = count+1
                  call twh_excludelist(GLB_SET,imolty,m,count,n)
               endif
            enddo
c           --- set the number of excluded atoms
            call twh_nexclude(GLB_SET,imolty,m,count)
         enddo
      enddo

      return
      end

      function twh_linclude(imolty,m,n)
c     ******************************************************************
c     * this function returns as true if atoms m and n for moltyp      *
c     * imolty have full nonbonded interactions with each other.       *
c     * It returns as false if the two atoms are the same, connected   *
c     * by a bond vibration, bending angle, or proper torsional angle  *
c     *                                                                *
c     * originally written 05-16-2001 by M.G. Martin                   *
c     * last modified 12-03-2007 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_LINCLUDE
#include "functions.h"

c     --- variables passed to/from the function
      integer imolty,m,n
c     --- local variables
      integer nb,excludecount,test

      call twh_nexclude(GLB_GET,imolty,m,excludecount)
      do nb = 1,excludecount
         call twh_excludelist(GLB_GET,imolty,m,nb,test)
         if ( test .eq. n ) then
            twh_linclude = .false.
            return
         endif
      enddo

c     --- if we made it down here then include these interactions
      twh_linclude = .true.
      return
      end
