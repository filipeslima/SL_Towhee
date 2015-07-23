#include "preproc.h"
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
      function twh_dotproduct( xone,yone,zone,xtwo,ytwo,ztwo,distone
     &     ,disttwo )
c     ******************************************************************
c     * computes the dot product of Vone * Vtwo and returns Cos(theta) *
c     * and the distances                                              *
c     * xone: x for vector one (passed in)                             *
c     * yone: y for vector one (passed in)                             *
c     * zone: z for vector one (passed in)                             *
c     * xtwo: x for vector two (passed in)                             *
c     * ytwo: y for vector two (passed in)                             *
c     * ztwo: z for vector two (passed in)                             *
c     * distone: length of vector one (returned)                       *
c     * disttwo: length of vector two (returned)                       *
c     * dotproduct: the Cos[theta] (returned as the function)          *
c     *                                                                *
c     * originally written 09-21-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_DOTPRODUCT
#include "functions.h"

c     --- variables passed to/from the function
      double precision xone,yone,zone,xtwo,ytwo,ztwo
     &     ,distone,disttwo,prod

      distone = twh_distance(xone,yone,zone)
      disttwo = twh_distance(xtwo,ytwo,ztwo)
      prod = distone*disttwo
      if ( prod .ne. 0.0d0 ) then
         twh_dotproduct = (xone*xtwo + yone*ytwo + zone*ztwo)/(prod)
      else
         twh_dotproduct = 0.0d0
      endif
      return
      end
