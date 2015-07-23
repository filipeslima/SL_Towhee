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
      subroutine twh_getweight(lfailure,style,vin,wout,dpscale,lterm )
c     ******************************************************************
c     * turns an energy into a rosenbluth weight that has been scaled  *
c     * dpscale times by scalelog in order that wout is nonzero within *
c     * the machine precision                                          *
c     * style: 0 input is an engery, returns weight and scale          *
c     * style: 1 input is a weight, returns adjusted weight and scale  *
c     *                                                                *
c     * originally written 05-21-2001 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure,lterm
      integer style
      double precision vin,wout,dpscale
c     --- parameters
      double precision maxint
      parameter (maxint = 1.0d9)
c     --- local variables
      integer isafety,iscale
      double precision vratio,beta,vequiv,scalelolog,scalehilog

      if ( style .eq. 0 ) then
c        --- input is an energy, turn into a weight and scale
c        --- make sure that the ratio of energies does not exceed
c        --- the maximum integer
         call twh_vequiv(GLB_GET,vequiv)
         vratio = vin/vequiv
         call twh_beta(GLB_GET,beta)
         if ( vin .le. 0.0d0 ) then
            if ( -vratio .gt. maxint ) then
c              --- energy is overwhelmingly negative: this really should
c              --- not happen and it messes up the code
               write(6,*) 'GETWEIGHT: vin is overwhelmingly negative'
               write(6,*) 'not sure how this could happen, vin',vin
               lfailure = .true.
               return
            endif
c           --- scale to keep the weight under control
            iscale = -idint(-vratio)
            dpscale = dble(iscale)
            wout = twh_expon(-beta*(vin-(dpscale*vequiv)))
         else
            if ( vratio .gt. maxint ) then
c              --- extremely high energy, treat as an overlap
               lterm = .true.
               return
            endif
c           --- scale to make sure we don't get negative weights
            iscale = idint(vratio)
            dpscale = dble(iscale)
            wout = twh_expon(-beta*(vin-(dpscale*vequiv)))
         endif
      elseif ( style .eq. 1 ) then
c        --- input is a weight, rescale the weight and return
c        --- if it is still possible to bring this weight back up
c        --- then do so.  If it takes 10 interations and the weight
c        --- is still not back then give up
         dpscale = 0.0d0
         isafety = 0
         call twh_scalelolog(GLB_GET,scalelolog)
         call twh_scalehilog(GLB_GET,scalehilog)
         do while ( wout .lt. scalelolog .and. isafety .ne. 10)
c           --- scale up the system
            dpscale = dpscale + 1.0d0
            isafety = isafety + 1
c           --- multiply by scalehilog to adjust
            wout = wout*scalehilog
         enddo
         do while ( wout .gt. scalehilog .and. isafety .ne. 10)
c           --- scale down the system
            dpscale = dpscale - 1.0d0
            isafety = isafety + 1
c           --- multiply by scalelolog to adjust
            wout = wout*scalelolog
         enddo
         if ( isafety .eq. 10 ) lterm = .true.
      else
         write(6,*) 'GETWEIGHT: invalid style ',style
         lfailure = .true.
         return
      endif

      return
      end
