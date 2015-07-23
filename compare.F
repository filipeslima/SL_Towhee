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
 
c     ******************************************************************
c     * This file defines a class of comparsion functions which are    *
c     * intended to provide reproducible results across a variety of   *
c     * architectures, to be used primarily for debugging and examples.*
c     * For production runs, standard .eq., .gt. and .lt. operators    *
c     * replace the function calls.                                    *
c     * The operation of this file is controlled by the SAFE_COMPARE   *
c     * compiler directive.  If SAFE_COMPARE is not defined, simple    *
c     * macro substitutions for .eq., .gt. and .lt. are made (see      *
c     * compare.h).                                                    *
c     *                                                                *
c     * Generally speaking, the SAFE_COMPARE directive optimizes       *
c     * cross-platform consistency over performance.                   *
c     *                                                                *
c     * See                                                            *
c     * http://www.lahey.com/float.htm                                 *
c     * and                                                            *
c     * http://www.cygnus-software.com/papers/comparingfloats/         *
c     *        comparingfloats.htm                                     *
c     * for discussion of safe compares.                               *
c     ******************************************************************

c     ******************************************************************
c     * originally written 2-01-2005 by MAW                            *
c     * last modified 08-09-2011  by M.G. Martin                       *
c     ******************************************************************
#ifdef SAFE_COMPARE

      function twh_cmp_eq(A,B)
c     ******************************************************************
c     * returns true if A is approximately equal to B                  *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_EQ
#include "functions.h"
      double precision A,B,RE
c     adapted from cygnus-software site (URL above)

c     compare for absolute error
      if (abs(A-B) .lt. CMP_EPS_ABS) then
         twh_cmp_eq = .true.
      else
c     compare for relative error
         if (abs(B) > abs(A)) then
            RE = abs((A - B) / B)
         else
            RE = abs((B - A) / A)
         endif

         if (RE .le. CMP_EPS_REL) then
            twh_cmp_eq = .true.
         else
            twh_cmp_eq = .false.
         endif
      endif
      return
      end

      function twh_cmp_gt(A,B)
c     ******************************************************************
c     * returns true if A is definitely greater than B                 *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_EQ
#define FUNCTION_CMP_GT
#include "functions.h"
      double precision A,B

      if (twh_cmp_eq(A,B)) then
         twh_cmp_gt = .false.
      elseif (A .lt. B) then
         twh_cmp_gt = .false.
      else
         twh_cmp_gt = .true.
      endif
      return
      end

      function twh_cmp_lt(A,B)
c     ******************************************************************
c     * returns true if A is definitely less than B                    *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_EQ
#define FUNCTION_CMP_LT
#include "functions.h"
      double precision A,B
      if (twh_cmp_eq(A,B)) then
         twh_cmp_lt = .false.
      elseif (A .lt. B) then
         twh_cmp_lt = .true.
      else
         twh_cmp_lt = .false.
      endif
      return
      end

#else
c     If SAFE_COMPARE flag is not set, nothing is defined in this file.  Some
c     compilers don't like this, so we define a dummy program unit.
      block data dummy
      end
#endif

