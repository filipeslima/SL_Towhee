c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2006 Marcus G. Martin                       *
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
      function twh_lawofcosine(dist_a,dist_b,angcos )
c     ******************************************************************
c     * uses the law of cosines to turn two distances and the cosine of*
c     * the angle between them into the third distance of the triangle *
c     *                                                                *
c     * originally written 10-05-2005 by M.G. Martin                   *
c     * last modified 03-09-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_LAWOFCOSINE
#include "functions.h"

c     --- variables passed to the function
      double precision dist_a,dist_b,angcos
c     --- local variables
      double precision lawcosa,lawcosb,distance

      lawcosa = dist_a*dist_a + dist_b*dist_b
      lawcosb = 2.0d0 * dist_a * dist_b
      distance = lawcosa - lawcosb*angcos
      if ( twh_cmp_gt(distance, 0.0d0) ) then
         twh_lawofcosine = dsqrt(distance)
      else
         twh_lawofcosine = 0.0d0
      endif

      return
      end
