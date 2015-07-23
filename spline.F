#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_splineset( splinetype,itype,jtype)
c     ******************************************************************
c     * sets up the cubic spline information for the embedded atom     *
c     * potential                                                      *
c     *                                                                *
c     * originally taken from Peter Schultz (Quest) on 11-18-2002      *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer itype,jtype,splinetype
c     --- local variables 
      integer i,ntot,cross
      double precision dx,dy,dx1,d2x,sig,p,dy1,qn,un,u,ddy0,ddy1
c     --- double precision arrays
      dimension u(MAXTABLE)
      double precision x
      dimension x(MAXTABLE)
      double precision y
      dimension y(MAXTABLE)
      double precision y2
      dimension y2(MAXTABLE)

      if ( splinetype .eq. 1 ) then
c        --- eam_dens, force to zero at far end
         ddy0 = 1.0d30
         ddy1 = 0.0d0
         call twh_eam_ndens(GLB_GET,itype,jtype,ntot)
c        --- set the x and y arrays
         do i = 1,ntot
            call twh_eam_dens(GLB_GET,1,itype,jtype,i,x(i))
            call twh_eam_dens(GLB_GET,2,itype,jtype,i,y(i))
         enddo
      elseif ( splinetype .eq. 2 ) then
c        --- eam_embed, force to zero at far end
         ddy0 = 0.0d0
         ddy1 = 0.0d0
         call twh_eam_nembed(GLB_GET,itype,ntot)
c        --- set the x and y arrays
         do i=1,ntot
            call twh_eam_embed(GLB_GET,1,itype,i,x(i))
            call twh_eam_embed(GLB_GET,2,itype,i,y(i))
         enddo
      elseif ( splinetype .eq. 3 ) then
c        --- table_pair, force derivaties to zero at far end
         ddy0 = 1.0d30
         ddy1 = 0.0d0
         cross = twh_getnbtype(itype,jtype)
         call twh_table_npair(GLB_GET,cross,ntot)
c        --- set the x and y arrays
         do i=1,ntot
            call twh_table_pair(GLB_GET,1,cross,i,x(i))
            call twh_table_pair(GLB_GET,2,cross,i,y(i))
         enddo
      else
         write(6,*) 'SPLINESET: invalid splinetype',splinetype
         stop
      endif

c     --- compute the spline information for table_pair for the spline
      if( ddy0 .gt. .99d30 ) then
        y2(1) = 0.d0
        u(1) = 0.d0
      else
        y2(1) = -.5d0
        dx = x(2) - x(1)
        dy = y(2) - y(1)
        u(1) = 3.d0*(dy/dx - ddy0)/dx
      endif
c     
      do i=2,ntot-1
         dx = x(i) - x(i-1)
         dx1 = x(i+1) - x(i)
         d2x = x(i+1) - x(i-1)
         sig = dx/d2x
         p = sig*y2(i-1) + 2.d0
         y2(i) = (sig - 1.d0)/p
         dy = y(i) - y(i-1)
         dy1 = y(i+1) - y(i)
         u(i) = ( 6.d0*(dy1/dx1-dy/dx)/d2x - sig*u(i-1) )/p
      enddo
c     
      if( ddy1 .gt. .99d30 )then
        qn = 0.d0
        un = 0.d0
      else
        qn = .5d0
        dx = x(ntot) - x(ntot-1)
        dy = y(ntot) - y(ntot-1)
        un = (ddy1 - dy/dx)*3.d0/dx
      endif
c     
      y2(ntot) = (un-qn*u(ntot-1)) / (qn*y2(ntot-1)+1.d0)
c     
c     --- back substitute
      do i=ntot-1,1,-1
         y2(i) = y2(i)*y2(i+1) + u(i)
      enddo

c     --- set the eam arrays
      do i = 1,ntot
         if ( splinetype .eq. 1 ) then
            call twh_eam_dens(GLB_SET,3,itype,jtype,i,y2(i))
         elseif ( splinetype .eq. 2 ) then
            call twh_eam_embed(GLB_SET,3,itype,i,y2(i))
         elseif ( splinetype .eq. 3 ) then
            call twh_table_pair(GLB_SET,3,cross,i,y2(i))
         endif
      enddo

      return
      end

