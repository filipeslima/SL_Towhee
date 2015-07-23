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
      subroutine twh_arbtocart (ibox,arba,arbb,arbc,rxuij,ryuij,rzuij)
c     ******************************************************************
c     * transforms from the arbitrary coordinate system (on 0,1) to the*
c     * cartesian coordinate system                                    *
c     *                                                                *
c     * originally written 12-13-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      integer ibox
      double precision arba,arbb,arbc,rxuij,ryuij,rzuij
      double precision hm11,hm12,hm13,hm21,hm22,hm23,hm31,hm32,hm33
c     --- fetch the hmatrix
      call twh_hmatrix(GLB_GET,ibox,1,1,hm11)
      call twh_hmatrix(GLB_GET,ibox,1,2,hm12)
      call twh_hmatrix(GLB_GET,ibox,1,3,hm13)
      call twh_hmatrix(GLB_GET,ibox,2,1,hm21)
      call twh_hmatrix(GLB_GET,ibox,2,2,hm22)
      call twh_hmatrix(GLB_GET,ibox,2,3,hm23)
      call twh_hmatrix(GLB_GET,ibox,3,1,hm31)
      call twh_hmatrix(GLB_GET,ibox,3,2,hm32)
      call twh_hmatrix(GLB_GET,ibox,3,3,hm33)
c     --- transform from arbitrary to cartesian coordinates
      rxuij = hm11*arba + hm21*arbb + hm31*arbc
      ryuij = hm12*arba + hm22*arbb + hm32*arbc
      rzuij = hm13*arba + hm23*arbb + hm33*arbc

      return
      end
