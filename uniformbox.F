c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2011 Marcus G. Martin                       *
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
      subroutine twh_uniformbox(ibox,rx,ry,rz)
c     ******************************************************************
c     * selects an atom uniformly from the volume of a simulation box  *
c     *                                                                *
c     * originally written 12-13-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_RANDOM
#include "functions.h"

c     --- variables passed to/from the subroutine
      integer ibox
      double precision rx,ry,rz
c     --- local variables
      double precision arba,arbb,arbc

      if ( ibox .eq. 0 ) then
c        --- ideal gas box - just return origin
         rx = 0.0d0
         ry = 0.0d0
         rz = 0.0d0
      else
c        --- select uniformly from the hmatrix cordinates
         arba = twh_random()
         arbb = twh_random()
         arbc = twh_random()
c        --- transform the inverse coordinates into cartesian space
         call twh_arbtocart(ibox,arba,arbb,arbc,rx,ry,rz)
      endif

      return
      end
