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
      subroutine twh_inverthmatrix (ibox)
c     ******************************************************************
c     * inverts the hmatrix for a particular box                       *
c     * also computes the boxvolume and stores it                      *
c     *                                                                *
c     * originally written 12-13-2001 by M.G. Martin                   *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_CROSSPRODUCT
#define FUNCTION_DISTANCE
#include "functions.h"
c     --- variables passed to/from the subroutine
      integer ibox

c     --- local variables
      logical lrectbox
      integer idim,jdim
      double precision tolerance
      parameter(tolerance=1.0d-7)
      double precision determ,test,hm11,hm12,hm13,hm21,hm22,hm23
      double precision hm31,hm32,hm33,dvalue,boxvolume
c     --- double precision arrays
      double precision hinvtemp
      dimension hinvtemp(3,3)
      double precision wx,wy,wz
      dimension wx(3),wy(3),wz(3)

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

c     --- invert the hmatrix
      hinvtemp(1,1) = hm22*hm33 - hm23*hm32
      hinvtemp(2,1) = hm31*hm23 - hm33*hm21
      hinvtemp(3,1) = hm21*hm32 - hm22*hm31
      hinvtemp(1,2) = hm32*hm13 - hm33*hm12
      hinvtemp(2,2) = hm11*hm33 - hm13*hm31
      hinvtemp(3,2) = hm31*hm12 - hm32*hm11
      hinvtemp(1,3) = hm12*hm23 - hm13*hm22
      hinvtemp(2,3) = hm21*hm13 - hm23*hm11
      hinvtemp(3,3) = hm11*hm22 - hm12*hm21
      
c     --- compute determinant (boxvolume) and make sure it is non-zero
      determ = hm11*hinvtemp(1,1) + hm12*hinvtemp(2,1)
     &     + hm13*hinvtemp(3,1)
      call twh_boxvolume(GLB_SET,ibox,determ)
      if ( dabs(determ) .lt. tolerance ) then
         write(6,*) 'INVERTHMATRIX: zero boxvolume ',determ
      else
         do idim = 1,3
            do jdim = 1,3
               hinvtemp(idim,jdim) = hinvtemp(idim,jdim) / determ
            enddo
         enddo
      endif
c     --- set hinverse
      do idim = 1,3
         do jdim = 1,3
            dvalue = hinvtemp(idim,jdim)
            call twh_hinverse(GLB_SET,ibox,idim,jdim,dvalue)
         enddo
      enddo

c     --- determine if this is a rectangular parallelapiped that is
c     --- aligned with global coordinate system (ie diagonal hmatrix)
      call twh_lrectbox(GLB_SET,ibox,.true.)
      test = 0.0d0
      do idim = 1,3
         do jdim = 1,3
            if ( idim .ne. jdim ) then
               call twh_hmatrix(GLB_GET,ibox,idim,jdim,dvalue)
               test = test + dvalue**2
            endif
         enddo
      enddo

c     --- set halfbox and check rectbox
      if ( test .gt. tolerance ) then
         call twh_lrectbox(GLB_SET,ibox,.false.)
         call twh_crossproduct(hm21,hm22,hm23,hm31,hm32,hm33
     &        ,wx(1),wy(1),wz(1))
         call twh_crossproduct(hm11,hm12,hm13,hm31,hm32,hm33
     &        ,wx(2),wy(2),wz(2))
         call twh_crossproduct(hm11,hm12,hm13,hm21,hm22,hm23
     &        ,wx(3),wy(3),wz(3))
c        --- retrieve the current boxvolume
         call twh_boxvolume(GLB_GET,ibox,boxvolume)
         do jdim = 1,3
c           --- e.g. width perpendicular to BC plane = (A.BxC)/|BxC|
            dvalue = 0.5d0 * boxvolume
     &           / twh_distance(wx(jdim),wy(jdim),wz(jdim))
            call twh_halfbox(GLB_SET,ibox,jdim,dvalue)
         enddo
      else
         do idim = 1,3
            jdim = idim
            call twh_hmatrix(GLB_GET,ibox,idim,jdim,dvalue)
            dvalue = 0.5d0*dvalue
            call twh_halfbox(GLB_SET,ibox,idim,dvalue)
         enddo
      endif

c     --- if it is not a rectangular parallelapiped don't use 
c     --- the additional COM cutoff as the logic can break down in
c     --- this case
      call twh_lrectbox(GLB_GET,ibox,lrectbox)
      if ( .not. lrectbox ) call twh_lusecom(GLB_SET,.false.)

      return
      end
