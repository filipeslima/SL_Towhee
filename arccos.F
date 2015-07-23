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
      function twh_arccos( value )
c     ******************************************************************
c     * computes the arc cosine of a value, with safety checks to make *
c     * sure that value is between -1.0 and 1.0.  If value is outside  *
c     * of that range, it is set to the nearest extreme of the range   *
c     *                                                                *
c     * originally written 02-12-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- need the arccos function
#define FUNCTION_ARCCOS
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#include "functions.h"

c     --- local variables
      double precision value

      if ( twh_cmp_gt(value, -1.0d0) ) then
         if ( twh_cmp_lt(value, 1.0d0) ) then
            twh_arccos = dacos(value)
         else
            twh_arccos = 0.0d0
         endif
      else
c        --- set the arccos to Pi
         call twh_constant_pi(GLB_GET,twh_arccos)
      endif

      return
      end
