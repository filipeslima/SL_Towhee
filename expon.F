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
      function twh_expon( value )
c     ******************************************************************
c     * computes the exponential of a value with safety checks to avoid*
c     * underflow.                                                     *
c     * Any value that would result in an exponential smaller than     *
c     * SMALLEST gets returns as 0.0 instead                           *
c     *                                                                *
c     * originally written 11-21-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#include "functions.h"

c     --- variables passed to the function
      double precision value,underflow
      
      call twh_constant_underflow(GLB_GET,underflow)
      if ( value .lt. underflow ) then
         twh_expon = 0.0d0
      else
         twh_expon = dexp(value)
      endif
      
      return
      end
