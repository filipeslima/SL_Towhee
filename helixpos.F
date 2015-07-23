c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2011 Marcus G. Martin                       *
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
      subroutine twh_helixpos(lfailure,helixradius,helixangle,bondlength
     &     ,helixphase,index,x)
c     ******************************************************************
c     * computes the Cartesian coordinates of any indexed position     *
c     * of a helix around the z-axis starting at the origin            *
c     *                                                                *
c     *   needs to receive the parameters:                             *
c     *   helixradius     = radius of the helix                        *
c     *   helixangle      = angle the helix makes wrt the z-axis       *
c     *   bondlength      = contour length between nodes in the helix  *
c     *   helixphase      = starting angle in degrees for helix        *
c     *   index           = position number in helix                   *
c     *                                                                *
c     *   returns the Cartesian coordinates of the indexed position in *
c     *   helix using the array x(3)                                   *
c     *                                                                *
c     * originally written 07-01-2004 by Steve Lustig                  *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer index
      double precision helixradius,helixangle,bondlength,helixphase,x
      dimension x(3)
c     --- local variables
      double precision pangle,b,q,cosangle,h,angle,di,pi_180
      double precision phase

      di = dble(index-1)
      if((bondlength*helixradius*helixangle.eq.0.0)
     &   .or.(helixangle.eq.90.0)) then
         write(6,*) 'HELIXPOS: This helix assignment is invalid'
         lfailure = .true.
         return
      endif
      pi_180 = dasin(1.0d0)/90.0d0
      pangle = helixangle*pi_180
      phase  = helixphase*pi_180
      b = bondlength*dcos(pangle)
      q = bondlength*dsin(pangle)/helixradius
      cosangle = 1.0d0-0.5d0*q*q
      h = b/dacos(cosangle)
      angle= di*b/h
      x(1) = helixradius*dcos(angle+phase)
      x(2) = helixradius*dsin(angle+phase)
      x(3) = h*angle
      return
      end
