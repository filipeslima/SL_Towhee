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
      function twh_gaussprob(value,mean,stddev,lolimit,hilimit)
c     ******************************************************************
c     * returns the probability density for a value in a gaussian      *
c     * distribution with mean of mean and standard deviation of stddev*
c     * this probability density is now normallized properly using the *
c     * 1/sqrt(2 Pi) prefactor and also the limits                     *
c     *                                                                *
c     * originally written 09-17-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DERF
#define FUNCTION_EXPON
#define FUNCTION_GAUSSPROB
#include "functions.h"
c     --- variables passed to/from the function
      double precision value,mean,stddev,lolimit,hilimit
c     --- local variables
      double precision normal,normalize,hiscale,loscale,twopi

      hiscale = (hilimit-mean)/(dsqrt(2.0d0)*stddev)
      loscale = (lolimit-mean)/(dsqrt(2.0d0)*stddev)
      call twh_constant_twopi(GLB_GET,twopi)
      normalize = 0.5d0 * ( twh_derf(hiscale) - twh_derf(loscale) )
     &     *stddev*dsqrt(twopi)
      normal = (value-mean)/stddev
      twh_gaussprob = twh_expon(-0.5d0*normal*normal)/normalize

      return
      end
