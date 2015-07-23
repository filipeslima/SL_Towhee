c     $Id: coneangle.F,v 1.4 2006/02/09 20:47:32 mattwycz Exp $
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2005 Marcus G. Martin                       *
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
      subroutine twh_coneangle(thetaone,phione,thetatwo,phitwo,cosangle)
c     ******************************************************************
c     * takes two unit vectors in spherical coordinates and computes   *
c     * the angle between them.                                        *
c     * thetaone is the angle between vector one and the z-axis        *
c     * phione is the angle between vector one and the x-axis          *
c     * thetatwo is the angle between vector two and the z-axis        *
c     * phitwo is the angle between vector two and the x-axis          *
c     * cosangle is the cos(angle) between the two vectors             *
c     * x = r sin (theta) cos (phi)                                    *
c     * y = r sin (theta) sin (phi)                                    *
c     * z = r cos (theta)                                              *
c     *                                                                *
c     * originally written 02-04-1998 by M.G. Martin                   *
c     * last modified 06-15-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none

c     --- variables passed to/from the subroutine
      double precision thetaone,thetatwo,phione,phitwo,cosangle

c     --- local variables
      double precision sintheone,costheone,sinthetwo,costhetwo
     &     ,sinphione,cosphione,sinphitwo,cosphitwo

      sintheone = dsin(thetaone)
      costheone = dcos(thetaone)
      sinthetwo = dsin(thetatwo)
      costhetwo = dcos(thetatwo)
      sinphione = dsin(phione)
      cosphione = dcos(phione)
      sinphitwo = dsin(phitwo)
      cosphitwo = dcos(phitwo)
      

      cosangle = sintheone*cosphione*sinthetwo*cosphitwo
     &     + sintheone*sinphione*sinthetwo*sinphitwo
     &     + costheone*costhetwo

      return
      end



