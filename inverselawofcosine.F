c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2010-2010 Marcus G. Martin                       *
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
      function twh_inverselawofcosine(dist_a,dist_b,dist_c )
c     ******************************************************************
c     * uses the law of cosines to turn three distances into one of the*
c     * angles.  In this case, it turns the a, b, and c distances into *
c     * the Cos(angle) on the triangle that is opposite of the c side  *
c     *                                                                *
c     * originally written 07-25-2010 by M.G. Martin                   *
c     * last modified 11-16-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_INVERSELAWOFCOSINE
#include "functions.h"

c     --- variables passed to the function
      double precision dist_a,dist_b,dist_c

      twh_inverselawofcosine = ( (dist_a**2) + (dist_b**2)
     &     -(dist_c**2))/(2.0d0*dist_a*dist_b)

      return
      end
