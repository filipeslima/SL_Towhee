#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2010 Marcus G. Martin                       *
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
      subroutine twh_initialize_loop_closure(lfailure,ibox,initcoords
     &     ,len0,b_ang0,t_ang0,aa13_min_sqr,aa13_max_sqr,delta,lcxi
     &     ,len_aa,eta)
c     ******************************************************************
c     * initialize_loop_closure takes the initial coordinates of       *
c     * the atoms involved in the concerted rotation move and          *
c     * initializes arrays which are required to solve the loop        *
c     * closure problem                                                *
c     *                                                                *
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 11-16-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_DOTPRODUCT
#define FUNCTION_VANGLE
#define FUNCTION_VTORSION
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ibox
      double precision aa13_min_sqr,aa13_max_sqr
      double precision initcoords
      dimension initcoords(9,3)
      double precision len0
      dimension len0(6)
      double precision b_ang0
      dimension b_ang0(7)
      double precision t_ang0
      dimension t_ang0(2)
      double precision delta
      dimension delta(0:3)
      double precision lcxi
      dimension lcxi(3)
      double precision len_aa
      dimension len_aa(3)
      double precision eta
      dimension eta(3)
c     --- local varibles
      integer i,j,k,ii
     
      double precision xcord,ycord,zcord,distone,disttwo
      double precision temp,rr_c1,conrotaxis,crang,torz,dumang
      double precision crbondtemp,lena1n1,lenc3a3,lenab,cosbang
      double precision rr_a1,rr_n2,rr_c1a1,rr_n2a2_ref,crp,us,rr_a2
      double precision rr_a1a2,dr,len1,len2,onepi
      double precision a_min,a_max,rr_a2_temp,temp33,torx,tory      

      dimension xcord(3),ycord(3),zcord(3)
      dimension temp(3),rr_a2_temp(3),temp33(3,3),torx(4)
      dimension torz(4),tory(4)
      dimension conrotaxis(3),rr_c1(3),crbondtemp(3,8),rr_a1(3),rr_n2(3)
      dimension rr_c1a1(3),rr_n2a2_ref(3),crp(4),us(3,3),rr_a2(3)
      dimension rr_a1a2(3),dr(3)

      call twh_constant_pi(GLB_GET,onepi)

      temp(1) = 0.0
      temp(2) = 0.0
      temp(3) = 0.0
      lena1n1 = 0.0d0

c     --- calculate the bond lengths and the angles which are 
c     --- required to initialize the various arrays used to solve 
c     --- the loop closure problem
      do j = 2,9
         i = j-1
         ii = j-2
         do k = 1,3
            crbondtemp(k,i) = initcoords(j,k) - initcoords(i,k)
         enddo
         crang = twh_dotproduct(crbondtemp(1,i),crbondtemp(2,i),
     &          crbondtemp(3,i),crbondtemp(1,i),crbondtemp(2,i),
     &          crbondtemp(3,i),distone,disttwo)
         if (j .eq. 2) then
            lena1n1 = distone
         elseif (j .eq. 9) then
            lenc3a3 = distone
         else 
            len0(ii) = distone         
         endif
      enddo
  
      do j = 1,8
         if (j .eq. 1) then
            lenab = lena1n1
         else if (j .eq. 8) then
            lenab = lenc3a3
         else 
            i = j-1
            lenab = len0(i)
         endif
         do k = 1,3
            crbondtemp(k,j) = crbondtemp(k,j)/lenab
         enddo
      enddo

      do j = 2,8
         i = j - 1
         cosbang = twh_dotproduct(-1.0d0*crbondtemp(1,i),
     &           -1.0d0*crbondtemp(2,i),
     &           -1.0d0*crbondtemp(3,i),crbondtemp(1,j),
     &           crbondtemp(2,j),
     &           crbondtemp(3,j),distone,disttwo)
         b_ang0(i) = twh_arccos(cosbang)
      enddo

c     --- calculate the dihedrals which are required to initialize the 
c     --- various arrays used to solve the loop closure problem      

      torx(1) = initcoords(2,1)
      torx(2) = initcoords(3,1)
      torx(3) = initcoords(4,1)
      torx(4) = initcoords(5,1)
      tory(1) = initcoords(2,2)
      tory(2) = initcoords(3,2)
      tory(3) = initcoords(4,2)
      tory(4) = initcoords(5,2)
      torz(1) = initcoords(2,3)
      torz(2) = initcoords(3,3)
      torz(3) = initcoords(4,3)
      torz(4) = initcoords(5,3)
 
      t_ang0(1) = -1.0d0 * twh_vtorsion(lfailure,ibox,0,0,0,0
     &     ,torx,tory,torz,dumang,.false.)
      if ( lfailure ) return

      torx(1) = initcoords(5,1)
      torx(2) = initcoords(6,1)
      torx(3) = initcoords(7,1)
      torx(4) = initcoords(8,1)
      tory(1) = initcoords(5,2)
      tory(2) = initcoords(6,2)
      tory(3) = initcoords(7,2)
      tory(4) = initcoords(8,2)
      torz(1) = initcoords(5,3)
      torz(2) = initcoords(6,3)
      torz(3) = initcoords(7,3)
      torz(4) = initcoords(8,3)
 
      t_ang0(2) = -1.0d0 * twh_vtorsion(lfailure,ibox,0,0,0,0
     &     ,torx,tory,torz,dumang,.false.)
      if ( lfailure ) return

      rr_c1(1) = 0.0d0
      conrotaxis(1) = 1.0d0
      rr_c1(2) = 0.0d0
      conrotaxis(2) = 0.0d0
      rr_c1(3) = 0.0d0
      conrotaxis(3) = 0.0d0

      do i = 0, 1
         rr_a1(1) = cos(b_ang0(3*i+2))*len0(3*i+1)
         rr_a1(2) = sin(b_ang0(3*i+2))*len0(3*i+1) 
         rr_a1(3) = 0.0d0
         rr_n2(1) = len0(3*i+2)
         rr_n2(2) = 0.0d0
         rr_n2(3) = 0.0d0

         rr_c1a1(1) = rr_a1(1) - rr_c1(1)
         rr_c1a1(2) = rr_a1(2) - rr_c1(2)
         rr_c1a1(3) = rr_a1(3) - rr_c1(3)
         rr_n2a2_ref(1) = -cos(b_ang0(3*i+3))*len0(3*i+3)
         rr_n2a2_ref(2) = sin(b_ang0(3*i+3))*len0(3*i+3)
         rr_n2a2_ref(3) = 0.0d0
         call twh_quaternion(conrotaxis, t_ang0(i+1)*0.25d0, crp)
         call twh_rotation_matrix(crp, us)

         do k = 1,3
            do j = 1,3
               temp33(k,j) = 0.0d0
            enddo
         enddo
         temp33(1,1) = rr_n2a2_ref(1)
         temp33(2,1) = rr_n2a2_ref(2)
         temp33(3,1) = rr_n2a2_ref(3)
         call twh_crmatmul(us,temp33,0,temp,rr_a2_temp)
         rr_a2(1) = rr_a2_temp(1) + rr_n2(1)
         rr_a2(2) = rr_a2_temp(2) + rr_n2(2)
         rr_a2(3) = rr_a2_temp(3) + rr_n2(3)
         rr_a1a2(1) = rr_a2(1) - rr_a1(1)
         rr_a1a2(2) = rr_a2(2) - rr_a1(2)
         rr_a1a2(3) = rr_a2(3) - rr_a1(3)
         dr(1) = rr_a1a2(1)
         dr(2) = rr_a1a2(2)
         dr(3) = rr_a1a2(3)
         len2 = (dr(1)*dr(1))+(dr(2)*dr(2))+(dr(3)*dr(3))
         len1 = sqrt(len2)
c        ---! len_aa
         len_aa(i+2) = len1


c        ---! xi
         xcord(1) = rr_a1(1)
         ycord(1) = rr_a1(2)
         zcord(1) = rr_a1(3)
         xcord(2) = rr_a2(1)
         ycord(2) = rr_a2(2)
         zcord(2) = rr_a2(3)
         xcord(3) = rr_n2(1)
         ycord(3) = rr_n2(2)
         zcord(3) = rr_n2(3)
         lcxi(i+2) = twh_vangle(lfailure,.false.,ibox,0,0,0,-1,dumang
     &        ,distone,disttwo,xcord,ycord,zcord)
         if ( lfailure ) return

c        ---! eta
         xcord(1) = rr_a2(1)
         ycord(1) = rr_a2(2)
         zcord(1) = rr_a2(3)
         xcord(2) = rr_a1(1)
         ycord(2) = rr_a1(2)
         zcord(2) = rr_a1(3)
         xcord(3) = rr_c1(1)
         ycord(3) = rr_c1(2)
         zcord(3) = rr_c1(3)
         eta(i+1) = twh_vangle(lfailure,.false.,ibox,0,0,0,-1,dumang
     &        ,distone,disttwo,xcord,ycord,zcord)
         if ( lfailure ) return

c        ! delta: pi -  dih of N(1)CA(1)CA(3)C(3)

         torx(1) = rr_c1(1)
         torx(2) = rr_a1(1)
         torx(3) = rr_a2(1)
         torx(4) = rr_n2(1)
         tory(1) = rr_c1(2)
         tory(2) = rr_a1(2)
         tory(3) = rr_a2(2)
         tory(4) = rr_n2(2)
         torz(1) = rr_c1(3)
         torz(2) = rr_a1(3)
         torz(3) = rr_a2(3)
         torz(4) = rr_n2(3)
   
         delta(i+1) = onepi + twh_vtorsion(lfailure,ibox,0,0,0,0
     &        ,torx,tory,torz,dumang,.false.)
         if ( lfailure ) return
      enddo

      a_min = b_ang0(4) - (lcxi(2) + eta(2))
      a_max = b_ang0(4) + (lcxi(2) + eta(2))

c     ---! min/max of base length
      aa13_min_sqr = len_aa(2)**2 + len_aa(3)**2 - 
     &                2.0d0*len_aa(2)*len_aa(3)*cos(a_min)
      aa13_max_sqr = len_aa(2)**2 + len_aa(3)**2 - 
     &                2.0d0*len_aa(2)*len_aa(3)*cos(a_max)

      return
      end

