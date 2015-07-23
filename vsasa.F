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
      function twh_vsasa(ichain,ibox,cordtype)
c     ******************************************************************
c     * Computes the colvation energy using the solvent accessible     *
c     * surface area (SASA) implicit water model as defined in         *
c     * Ferrara, Apostolakis, and Caflisch, Proteins 46, 24-33 (2002)  *
c     *                                                                *
c     * originally written 10-31-2002 by M.G. Martin                   *
c     * last update 08-09-2011 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_DISTANCE
#define FUNCTION_GET_IJVIB
#define FUNCTION_GET_INVIB
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_VSASA
#include "functions.h"
c     --- variables passed to/from the function
      integer ichain,cordtype,ibox
c     --- local variables
      integer j,typej,inn,iunit,typei,imolty,nunit
      double precision rprobe,sasaprod,xi,yi,zi,risasa
      double precision pisasa,sisasa,sigmaisasa,xj,yj,zj
      double precision rjsasa,pijsasa,bijsasa,deltax,deltay
      double precision deltaz,dist,compare,prod,vsasai,onepi
      logical lneigh

      call twh_constant_pi(GLB_GET,onepi)
      call twh_moltyp(GLB_GET,ichain,imolty)
      call twh_nunit(GLB_GET,imolty,nunit)

      twh_vsasa = 0.0d0
      rprobe = 1.4d0
      sasaprod = 1.0d0
      do iunit = 1,nunit
         call twh_coordinates(GLB_GET,cordtype,ichain,iunit,xi,yi,zi)
         call twh_ntype(GLB_GET,imolty,iunit,typei)
         risasa = twh_get_nbcoeff(5,typei)
         sisasa = 4*onepi*(risasa+rprobe)*(risasa+rprobe)
         pisasa = twh_get_nbcoeff(6,typei)
         sigmaisasa = twh_get_nbcoeff(7,typei)
c        --- calculate the solvent assessible surface area of particle i
         do j = 1,nunit
            if (j .ne. iunit) then
               call twh_coordinates(GLB_GET,cordtype,ichain,j
     &              ,xj,yj,zj)
               call twh_ntype(GLB_GET,imolty,j,typej)
               rjsasa = twh_get_nbcoeff(5,typej)
               lneigh = .false.
               do inn = 1, twh_get_invib(imolty,iunit)
                  if (.not. lneigh) then
                     if (j .eq. twh_get_ijvib(imolty,iunit,inn)) then
                        lneigh = .true.
                     endif
                  endif
               enddo
               if (lneigh) then
                  pijsasa = 0.8875
               else 
                  pijsasa = 0.3516
               endif
               call twh_mimage(ibox,xj,yj,zj,xi,yi,zi
     &              ,deltax,deltay,deltaz)
               dist = twh_distance(deltax,deltay,deltaz)
               compare = risasa + rjsasa + (2*rprobe)
               if (dist .gt. compare) then
                  bijsasa = 0.0
               else
                  bijsasa = onepi * (risasa + rprobe) * 
     1               (risasa + rjsasa + (2 * rprobe) - dist) * 
     2               ( 1 + ((rjsasa - risasa)/dist))
               endif
               prod = pisasa * pijsasa * bijsasa/sisasa
               sasaprod = sasaprod * (1 - prod)
            endif
         enddo
         vsasai = sigmaisasa * sisasa * sasaprod
         twh_vsasa = twh_vsasa + vsasai
      enddo

      return
      end
