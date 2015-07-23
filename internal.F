#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2008 Marcus G. Martin                       *
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
      function twh_derfc(x)
c     ******************************************************************
c     * computes the complementary error function                      *
c     * normally just calls the fortran intrinsic derfc, but if the    *
c     * --enable-internal option is used during the configure it uses  *
c     * a polynomial expansion instead                                 *
c     *                                                                *
c     * originally written by Bin Chen in 1998 as erfunc               *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#ifdef USEINTERNAL
#define FUNCTION_EXPON
#endif
#define FUNCTION_DERFC
#include "functions.h"

c     --- variables passed to/from the function
      double precision x
c     --- local variables
#ifdef USEINTERNAL
      logical lsignswap
      double precision p,a1,a2,a3,a4,a5,tt,eee,posx
      parameter (p=0.3275911d0,a1=0.254829592d0,
     &     a2=-0.284496736d0,a3=1.421413741d0,
     &     a4=-1.453152027d0,a5=1.061405429d0)
#endif

#ifdef USEINTERNAL
      if ( x .lt. 0.0d0 ) then
c        --- switch the sign
         lsignswap = .true.
         posx = -x
      else
c        --- sign is fine
         lsignswap = .false.
         posx = x
      endif
      eee = twh_expon(-posx*posx)
      tt = 1.0d0/(1.0d0 + p*posx)
      twh_derfc = ((((a5*tt+a4)*tt+a3)*tt+a2)*tt+a1)*tt*eee
      if ( lsignswap ) twh_derfc = 2.0d0 - twh_derfc
#else
      twh_derfc = derfc(x)
#endif
      return
      end

      function twh_derf(x)
c     ******************************************************************
c     * computes the error function.  Normally uses the intrinsic, but *
c     * if configured with --enable-internal it instead uses the       *
c     * polynomial version of the complementary errory function and the*
c     * relation between these two functions                           *
c     *                                                                *
c     * originally written 01-05-2006 by M.G. Martin                   *
c     * last modified 03-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#ifdef USEINTERNAL
#define FUNCTION_DERFC
#endif
#define FUNCTION_DERF
#include "functions.h"
c     --- variables passed to/from the function
      double precision x

#ifdef USEINTERNAL
      twh_derf = 1.0d0 - twh_derfc(x)
#else
      twh_derf = derf(x)
#endif
      return
      end

      function twh_len_trim(teststring)
c     ******************************************************************
c     * returns the index of the last character before a blank         *
c     *                                                                *
c     * originally written 03-14-2006 by M.G. Martin                   *
c     * last modified 03-14-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_LEN_TRIM
#include "functions.h"
c     --- variable passed to/from the function
      character*(*) teststring
c     --- local variables
#ifdef USEINTERNAL
      character*1 testchar
      integer count,maxcount
      parameter (maxcount=1000)
      testchar = ' '
      count = 0
      do while ( teststring(count+1:count+1) .ne. testchar )
         count = count + 1
      enddo
      twh_len_trim = count
#else
      twh_len_trim = len_trim(teststring)
#endif

      return
      end
