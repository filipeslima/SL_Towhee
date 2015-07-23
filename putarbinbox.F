#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin,                      *
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
      subroutine twh_putarbinbox (arba,arbb,arbc)
c     ******************************************************************
c     * puts the atom inside the central simulation box                *
c     *                                                                *
c     * originally written 01-11-2002 by M.G. Martin                   *
c     * last modified 08-01-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      double precision arba,arbb,arbc
c     --- local variables
      logical ldone
      integer nshift,shiftmax
      parameter (shiftmax = 10)

      ldone = .false.
      nshift = 0
      do while ( .not. ldone .and. nshift .lt. shiftmax) 
         if ( arba .gt. 1.0d0 ) then
            arba = arba - 1.0d0
            nshift = nshift + 1
         elseif (arba .lt. 0.0d0 ) then
            arba = arba + 1.0d0
            nshift = nshift + 1
         else
            ldone = .true.
         endif
      enddo
      if ( nshift .gt. 1 ) then
         write(6,*) 'PUTARBINBOX: A-dim required multiple shifts',nshift
      endif

      ldone = .false.
      nshift = 0
      do while ( .not. ldone .and. nshift .lt. shiftmax) 
         if ( arbb .gt. 1.0d0 ) then
            arbb = arbb - 1.0d0
            nshift = nshift + 1
         elseif (arbb .lt. 0.0d0 ) then
            arbb = arbb + 1.0d0
            nshift = nshift + 1
         else
            ldone = .true.
         endif
      enddo
      if ( nshift .gt. 1 ) then
         write(6,*) 'PUTARBINBOX: B-dim required multiple shifts',nshift
      endif

      ldone = .false.
      nshift = 0
      do while ( .not. ldone .and. nshift .lt. shiftmax) 
         if ( arbc .gt. 1.0d0 ) then
            arbc = arbc - 1.0d0
            nshift = nshift + 1
         elseif (arbc .lt. 0.0d0 ) then
            arbc = arbc + 1.0d0
            nshift = nshift + 1
         else
            ldone = .true.
         endif
      enddo
      if ( nshift .gt. 1 ) then
         write(6,*) 'PUTARBINBOX: C-dim required multiple shifts',nshift
      endif

      return
      end
