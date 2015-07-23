#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2008 Marcus G. Martin                       *
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
      function twh_integratedgaussprob(lovalue,hivalue
     &     ,mean,stddev,lolimit,hilimit)
c     ******************************************************************
c     * returns the integrated probability that a value falls in the   *
c     * range of lovalue to hivalue inside of a Gaussian distribution  *
c     * with the mean, standard deviation, and limits passed into this *
c     * function                                                       *
c     *                                                                *
c     * originally written 10-06-2005 by M.G. Martin                   *
c     * last modified 10-13-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_DERF
#define FUNCTION_INTEGRATEDGAUSSPROB
#include "functions.h"
c     --- variables passed to/from the function
      double precision lovalue,hivalue,mean,stddev,lolimit,hilimit
c     --- local variables
      double precision hilimitscale,lolimitscale,scale
      double precision lovaluescale,hivaluescale

      scale = 1.0d0/(dsqrt(2.0d0)*stddev)
      hilimitscale = (hilimit-mean)*scale
      lolimitscale = (lolimit-mean)*scale
      hivaluescale = (hivalue-mean)*scale
      lovaluescale = (lovalue-mean)*scale
      twh_integratedgaussprob = 
     &     (twh_derf(hivaluescale) -twh_derf(lovaluescale))
     &     /(twh_derf(hilimitscale) - twh_derf(lolimitscale))

      return
      end
