c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_sphere ( x, y, z )
c     ******************************************************************
c     * calculates a random vector on the unit sphere                  *
c     * uses a rejection method to get the vector on the sphere        *
c     *                                                                *
c     * originally written 02-03-1998 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the subroutine
      double precision x, y, z

c     --- local variables
      integer ii
      double precision xi1,xi2,xisq

c     --- calculate random vector on the unit sphere 
      do ii = 1,100
         xi1 = ( 2.0d0* twh_random() ) - 1.0d0
         xi2 = ( 2.0d0* twh_random() ) - 1.0d0
         xisq = xi1**2 + xi2**2
         if ( xisq .lt. 1.0d0 ) then
            x = 2.0d0* xi1* dsqrt( 1.0d0 - xisq )
            y = 2.0d0* xi2* dsqrt( 1.0d0 - xisq )
            z = ( 1.0d0 - 2.0d0* xisq )
            return
         endif
      enddo

      stop 'exceeded 100 tries to get a vector in sphere'

      end
