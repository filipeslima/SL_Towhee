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
      function twh_gaussian(lfailure,mean,stddev,lolimit,hilimit)
c     ******************************************************************
c     * this function computes two random numbers that come from a     *
c     * standard normal distribution.  This algorithm was taken from a *
c     * web page (www.taygeta.com/random/gaussian.html)                *
c     *                                                                *
c     * originally written 09-17-2001 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GAUSSIAN
#define FUNCTION_RANDOM
#include "functions.h"
c     --- variables passed to/from the function
      logical lfailure
      double precision mean,stddev,lolimit,hilimit
c     --- local variables
      logical lfound
      integer attempts
      double precision x1,x2,w,first,extra

      lfound = .false.
      attempts = 0
      do while ( .not. lfound .and. attempts .lt. 50)
         attempts = attempts + 1
c        --- get two random variables from the standard normal
         w = 2.0
         do while ( w .ge. 1.0d0 )
            x1 = 2.0d0 * twh_random() -1.0d0
            x2 = 2.0d0 * twh_random() -1.0d0
            w = x1*x1 + x2*x2
         enddo
        
         w = dsqrt( (-2.0d0*dlog(w))/w)
         first = x1*w
         extra = x2*w
c        --- see if either of these variable fall within the desired
c        --- range
         twh_gaussian = mean + stddev*first
         if ( twh_gaussian .gt. lolimit .and. twh_gaussian .lt. hilimit 
     &        ) then
            lfound = .true.
         else
c           --- try the extra
            twh_gaussian = mean + stddev*extra
            if ( twh_gaussian .gt. lolimit 
     &           .and. twh_gaussian .lt. hilimit ) lfound = .true.
         endif
      enddo
      if ( attempts .eq. 50 ) then
         write(6,*) 'failed to get the gaussian number in 50 attempts'
         write(6,*) 'mean,stddev,lolimit,hilimit',mean,stddev,lolimit
     &        ,hilimit
         lfailure = .true.
         return
      endif

      return
      end
