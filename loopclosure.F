#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2008 Marcus G. Martin                       *
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
      subroutine twh_solv_3pep_poly(lfailure,max_soln,ibox,rootunit
     &     ,len0,b_ang0,solncoords,nsoln,aa13_min_sqr,aa13_max_sqr
     &     ,delta,lcxi,len_aa,eta)
c     ******************************************************************
c     * solv_3pep_poly takes the coordinates of the root atoms used    *
c     * in the moved (rootunit), the bond lengths (len0), and          *
c     * bond angles (b_ang0) used to define the chain in between the   *
c     * two rooted ends and then determines the number of solutions    *
c     * (nsoln) and the cartesian coordinates of the "nsoln" possible  *
c     * loop closure solutions (solncoords)                            *
c     *                                                                *
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- functions
#define FUNCTION_ARCCOS
#define FUNCTION_DOTPRODUCT
#define FUNCTION_VANGLE
#define FUNCTION_VTORSION
#include "functions.h"
      double precision twh_calc_t2,twh_calc_t1
c     --- paramaters passed to/from the subroutine
      integer max_soln
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision rootunit,solncoords,len0,b_ang0
      dimension rootunit(4,3),solncoords(9,3,max_soln),len0(6),b_ang0(7)
      integer nsoln,ibox
      double precision aa13_min_sqr,aa13_max_sqr
      double precision delta
      dimension delta(0:3)
      double precision lcxi
      dimension lcxi(3)
      double precision len_aa
      dimension len_aa(3)
      double precision eta
      dimension eta(3)
c     --- local parameters
      integer deg_pol
      parameter(deg_pol=16)
c     --- local variables
c     --- integer scalars
      integer i_soln,iii,jjj,kkk
c     --- double precision scalars
      double precision dr_sqr,dumang,dist1,dist2,onepi
      double precision lensq_n1a1,lensq_c3a3
      double precision lcangle,sig1_init
      double precision ht,tmp,sig1
      double precision lez
      double precision dcdl
      double precision distone,disttwo,cdlcos
c     --- double precision arrays
      double precision poly_coeff
      dimension poly_coeff(0:16)
      double precision roots
      dimension roots(max_soln)
      double precision r_a1a3
      dimension r_a1a3(3)
      double precision xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
      double precision torx,tory,torz
      dimension torx(4),tory(4),torz(4)
      double precision ex,ey,ez,b_a1a2,b_a3a2,r_tmp
      dimension ex(3),ey(3),ez(3),b_a1a2(3),b_a3a2(3),r_tmp(3)
      double precision p_s,s1,s2,p_t,t1,t2
      dimension p_s(3,3),s1(3,3),s2(3,3),p_t(3,3),t1(3,3),t2(3,3)
      double precision p_s_c,s1_s,s2_s,p_t_c,t1_s
      dimension p_s_c(3,3),s1_s(3,3),s2_s(3,3),p_t_c(3,3),t1_s(3,3)
      double precision half_tan,ex_tmp
      dimension half_tan(3),ex_tmp(3)
      double precision cos_tau,sin_tau,cos_sig
      dimension cos_tau(0:3),sin_tau(0:3),cos_sig(3)
      double precision sin_sig
      dimension sin_sig(3)
      double precision r_c,p,us
      dimension r_c(3,3),p(4),us(3,3)
      double precision r_s,r_t,r0,r_n,r_a
      dimension r_s(3),r_t(3),r0(3),r_n(3,3),r_a(3,3)
      double precision s11_tmp
      dimension s11_tmp(3)
      double precision rtmp
      dimension rtmp(3)
      double precision t2_s
      dimension t2_s(3,3)
      double precision lcalpha
      dimension lcalpha(3)
      double precision lctheta
      dimension lctheta(3)
      double precision cos_alpha
      dimension cos_alpha(3)
      double precision cos_theta
      dimension cos_theta(3)
      double precision sin_alpha
      dimension sin_alpha(3)
      double precision cos_delta
      dimension cos_delta(0:3)
      double precision sin_delta
      dimension sin_delta(0:3)
      double precision cos_xi
      dimension cos_xi(3)
      double precision cos_eta
      dimension cos_eta(3)
      double precision sin_xi
      dimension sin_xi(3)
      double precision sin_eta
      dimension sin_eta(3)
      double precision r_a1n1
      dimension r_a1n1(3)
      double precision r_a3c3
      dimension r_a3c3(3)
      double precision len_na
      dimension len_na(3)
      double precision len_ac
      dimension len_ac(3)
      double precision c0
      dimension c0(0:2,3)
      double precision c1
      dimension c1(0:2,3)
      double precision c2
      dimension c2(0:2,3)
      double precision crq
      dimension crq(0:16,0:4)
      double precision lcr
      dimension lcr(0:16,0:2)

c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)

c     --- start of get input angles section
c     --- determines the relavent bond lengths, bond angles, and
c     --- dihedral angles from the coordinates of the rooted units and
c     --- the bond lengths, bond angles, and torsion angles
      nsoln = max_soln
c     ---! virtual bond
      r_a1a3(1) = rootunit(3,1) - rootunit(2,1)
      r_a1a3(2) = rootunit(3,2) - rootunit(2,2)
      r_a1a3(3) = rootunit(3,3) - rootunit(2,3)
      dr_sqr = (r_a1a3(1)*r_a1a3(1))+(r_a1a3(2)*r_a1a3(2))
     &             +(r_a1a3(3)*r_a1a3(3))
      len_aa(1) = sqrt(dr_sqr)

      if (dr_sqr .lt. aa13_min_sqr .or. dr_sqr .gt. aa13_max_sqr) then
         nsoln = 0
         return
      end if
 
c     ---! bond lengths
      r_a1n1(1) = rootunit(1,1) - rootunit(2,1)
      r_a1n1(2) = rootunit(1,2) - rootunit(2,2)
      r_a1n1(3) = rootunit(1,3) - rootunit(2,3)
      lensq_n1a1 = (r_a1n1(1)*r_a1n1(1))+(r_a1n1(2)*r_a1n1(2))
     &             +(r_a1n1(3)*r_a1n1(3))
      len_na(1) = sqrt(lensq_n1a1)
      len_na(2) = len0(3)
      len_na(3) = len0(6)
      r_a3c3(1) = rootunit(4,1) - rootunit(3,1)
      r_a3c3(2) = rootunit(4,2) - rootunit(3,2)
      r_a3c3(3) = rootunit(4,3) - rootunit(3,3)
      lensq_c3a3 = (r_a3c3(1)*r_a3c3(1))+(r_a3c3(2)*r_a3c3(2))
     &             +(r_a3c3(3)*r_a3c3(3))
      len_ac(1) = len0(1)
      len_ac(2) = len0(4)
      len_ac(3) = sqrt(lensq_c3a3)

c     ---! delta(3): dih of N(1)CA(1)CA(3)C(3)
      torx(1) = rootunit(1,1)
      torx(2) = rootunit(2,1)
      torx(3) = rootunit(3,1)
      torx(4) = rootunit(4,1)
      tory(1) = rootunit(1,2)
      tory(2) = rootunit(2,2)
      tory(3) = rootunit(3,2)
      tory(4) = rootunit(4,2)
      torz(1) = rootunit(1,3)
      torz(2) = rootunit(2,3)
      torz(3) = rootunit(3,3)
      torz(4) = rootunit(4,3)

      delta(3) = -1.0d0 *twh_vtorsion(lfailure,ibox,0,0,0,0
     &     ,torx,tory,torz,dumang,.false.)
      if ( lfailure ) return

      delta(0) = delta(3)

c     ---! lcxi(1)
      xcord(1) = rootunit(3,1)
      xcord(2) = rootunit(2,1)
      xcord(3) = rootunit(1,1)

      ycord(1) = rootunit(3,2)
      ycord(2) = rootunit(2,2)
      ycord(3) = rootunit(1,2)

      zcord(1) = rootunit(3,3)
      zcord(2) = rootunit(2,3)
      zcord(3) = rootunit(1,3)
      lcxi(1) = onepi - twh_vangle(lfailure,.false.,ibox,0,0,0,-1
     &     ,dumang,dist1,dist2,xcord,ycord,zcord)
      if ( lfailure ) return

c     ---! eta(3)
      xcord(1) = rootunit(2,1)
      xcord(2) = rootunit(3,1)
      xcord(3) = rootunit(4,1)

      ycord(1) = rootunit(2,2)
      ycord(2) = rootunit(3,2)
      ycord(3) = rootunit(4,2)

      zcord(1) = rootunit(2,3)
      zcord(2) = rootunit(3,3)
      zcord(3) = rootunit(4,3)
      eta(3) = onepi - twh_vangle(lfailure,.false.,ibox,0,0,0,-1,dumang
     &     ,dist1,dist2,xcord,ycord,zcord)
      if ( lfailure ) return

      do iii = 1, 3
         cos_delta(iii) = cos(delta(iii))
         sin_delta(iii) = sin(delta(iii))
         cos_xi(iii) = cos(lcxi(iii))
         sin_xi(iii) = sin(lcxi(iii))
         cos_eta(iii) = cos(eta(iii))
         sin_eta(iii) = sin(eta(iii))
      enddo
      cos_delta(0) = cos_delta(3)
      sin_delta(0) = sin_delta(3)

c     ---! theta (N, CA, C) bond angle
      lctheta(1) = b_ang0(1)
      lctheta(2) = b_ang0(4)
      lctheta(3) = b_ang0(7)
      do iii = 1, 3
         cos_theta(iii) = cos(lctheta(iii))
      enddo

c     ! alpha
      cos_alpha(1) = -(len_aa(1)**2 + len_aa(2)**2 - 
     &                len_aa(3)**2)/(2.0d0*len_aa(1)*len_aa(2))
      lcalpha(1) = acos(cos_alpha(1))
      sin_alpha(1) = sin(lcalpha(1))
      cos_alpha(2) = (len_aa(2)**2 + len_aa(3)**2 - 
     &                len_aa(1)**2)/(2.0d0*len_aa(2)*len_aa(3))
      lcalpha(2) = acos(cos_alpha(2))
      sin_alpha(2) = sin(lcalpha(2))
      lcalpha(3) = onepi - lcalpha(1) + lcalpha(2)
      cos_alpha(3) = cos(lcalpha(3))
      sin_alpha(3) = sin(lcalpha(3))

c     --- check for existence of soln
      do iii = 1,3
         call twh_test_two_cone_exist_soln(max_soln,lctheta(iii)
     &        ,lcxi(iii),eta(iii),lcalpha(iii),nsoln)
         if (nsoln .eq. 0) return
      enddo
c     --- end of get input angles section

      call twh_get_poly_coeff(lfailure,poly_coeff,cos_alpha,cos_theta
     &     ,sin_alpha,cos_delta,sin_delta,cos_xi,cos_eta,sin_xi,sin_eta
     &     ,c0,c1,c2,crq,lcr)
      if ( lfailure ) return
 
      call solvesturm(deg_pol, nsoln, poly_coeff, roots)
      if (nsoln .eq. 0) return

c     --- start of the coord from poly roots section
c     --- Determines the solution coordinates from the roots of the
c     --- aforementioned 16th order polynomial
c     --- local variables 

      if (nsoln .eq. 0) return

c     --- Define body frame (ex, ey, ez)
      do kkk = 1,3
         ex(kkk) = r_a1a3(kkk)/len_aa(1)
      enddo
      call twh_crossproduct(r_a1n1(1),r_a1n1(2),r_a1n1(3),ex(1),ex(2),
     &               ex(3),ez(1),ez(2),ez(3))
c      call cross(r_a1n1, ex, ez)
      lez = dsqrt((ez(1)*ez(1))+(ez(2)*ez(2))+(ez(3)*ez(3)))

      do kkk = 1,3
         ez(kkk) = ez(kkk)/lez
      enddo
      call twh_crossproduct(ez(1),ez(2),ez(3),ex(1),ex(2),ex(3),ey(1),
     &                  ey(2),ey(3))
c      call cross(ez, ex, ey)
c     --- vertual bond vectors in the reference plane
      do kkk=1,3
         b_a1a2(kkk) = -1.0d0*cos_alpha(1)*ex(kkk)
     &        + sin_alpha(1)*ey(kkk)
         b_a3a2(kkk) = cos_alpha(3)*ex(kkk) + sin_alpha(3)*ey(kkk)
c     --- Define cone coordinates for each angle joint.
c     --- (p_s,s1,s2) and (p_t,t1,t2):  Right Orthonormal systems
c     --- residue 1
         p_s(kkk,1) = -1.0d0*ex(kkk)
c        --- (p_s)X(p_t)/||(p_s)X(p_t)||
         s1(kkk,1)  = ez(kkk)  
c        --- p_s X s1
         s2(kkk,1)  = ey(kkk)  
         p_t(kkk,1) = b_a1a2(kkk)
         t1(kkk,1)  = ez(kkk)  
c        --- p_t X t1
         t2(kkk,1)  = sin_alpha(1)*ex(kkk) + cos_alpha(1)*ey(kkk) 
c        --- residue 2
         p_s(kkk,2) = -1.0d0*b_a1a2(kkk)
c        --- sina1*ex(:) + cosa1*ey(:)
         s1(kkk,2)  = -1.0d0*ez(kkk)
         s2(kkk,2)  = t2(kkk,1)  
         p_t(kkk,2) = -1.0d0*b_a3a2(kkk)
         t1(kkk,2)  = -1.0d0*ez(kkk)
         t2(kkk,2)  = sin_alpha(3)*ex(kkk) - cos_alpha(3)*ey(kkk) 
c        --- residue 3
         p_s(kkk,3) = b_a3a2(kkk)
c        --- sina3*ex(:) + cosa3*ey(:) 
         s2(kkk,3)  = t2(kkk,2)   
         s1(kkk,3)  = ez(kkk)  
         p_t(kkk,3) = ex(kkk)
         t1(kkk,3) =  ez(kkk) 
         t2(kkk,3) = -1.0d0*ey(kkk)
c        --- scale vectors
         do iii = 1, 3
            p_s_c(kkk,iii) = p_s(kkk,iii)*cos_xi(iii)
            s1_s(kkk,iii)  =  s1(kkk,iii)*sin_xi(iii)
            s2_s(kkk,iii)  =  s2(kkk,iii)*sin_xi(iii)
            p_t_c(kkk,iii) = p_t(kkk,iii)*cos_eta(iii)
            t1_s(kkk,iii)  =  t1(kkk,iii)*sin_eta(iii)
            t2_s(kkk,iii)  =  t2(kkk,iii)*sin_eta(iii)
         end do
c         --- initial sig(1)
         r_tmp(kkk) = (r_a1n1(kkk)/len_na(1) - p_s_c(kkk,1))/sin_xi(1)
         s11_tmp(kkk) = s1(kkk,1)
      enddo

      cdlcos = twh_dotproduct(s11_tmp(1),s11_tmp(2),s11_tmp(3),
     &          r_tmp(1),r_tmp(2),r_tmp(3),distone,disttwo)

      lcangle = twh_arccos(cdlcos)

      dcdl = (r_tmp(1)*s2(1,1))+(r_tmp(2)*s2(2,1))
     &       +(r_tmp(3)*s2(3,1))
      sig1_init = sign(lcangle, dcdl)
      

      do kkk = 1,3      
c        --- CA
         r_a(kkk,1) = rootunit(2,kkk)
         r_a(kkk,2) = rootunit(2,kkk) + (len_aa(2)*b_a1a2(kkk))
         r_a(kkk,3) = rootunit(3,kkk)
         r0(kkk) = rootunit(2,kkk)
      enddo
      do i_soln = 1, nsoln 
         half_tan(3) = roots(i_soln)
         half_tan(2) = twh_calc_t2(half_tan(3),crq,lcr)
         half_tan(1) = twh_calc_t1(half_tan(3), half_tan(2),c0,c1,c2)
         do iii = 1, 3
            ht = half_tan(iii)
            tmp = 1.0d0 + ht*ht
            cos_tau(iii) = (1.0d0 - ht*ht)/tmp
            sin_tau(iii) = 2.0d0*ht/tmp
         enddo
         cos_tau(0) = cos_tau(3)
         sin_tau(0) = sin_tau(3)
         do iii = 1, 3
            jjj = iii - 1
            cos_sig(iii) = cos_delta(jjj)*cos_tau(jjj) + 
     &                   sin_delta(jjj)*sin_tau(jjj)
            sin_sig(iii) = sin_delta(jjj)*cos_tau(jjj) - 
     &                   cos_delta(jjj)*sin_tau(jjj)
         end do
         do kkk = 1,3
            do iii = 1, 3
               r_s(kkk) = p_s_c(kkk,iii) + cos_sig(iii)*s1_s(kkk,iii)+
     &                  sin_sig(iii)*s2_s(kkk,iii)
               r_t(kkk) = p_t_c(kkk,iii) + cos_tau(iii)*t1_s(kkk,iii)+
     &                  sin_tau(iii)*t2_s(kkk,iii) 
               r_n(kkk,iii) = r_s(kkk)*len_na(iii) + r_a(kkk,iii)
               r_c(kkk,iii) = r_t(kkk)*len_ac(iii) + r_a(kkk,iii)
            enddo
         enddo
      
c        --- rotate back atoms by -(sig(1) - sig1_init) around -ex
         sig1 = atan2(sin_sig(1), cos_sig(1))
         do kkk = 1,3
            ex_tmp(kkk) = -1.0d0*ex(kkk)
         enddo 
         call twh_quaternion(ex_tmp, -1.0d0*(sig1 - sig1_init)*0.25d0,p)
         call twh_rotation_matrix(p, us)
     
         do kkk = 1,3
            solncoords(1,kkk,i_soln) = rootunit(1,kkk)
            solncoords(2,kkk,i_soln) = rootunit(2,kkk)
            solncoords(8,kkk,i_soln) = rootunit(3,kkk)
            solncoords(9,kkk,i_soln) = rootunit(4,kkk)
         end do

         call twh_crmatmul(us,r_c,1,r0,rtmp)
         do kkk = 1,3
            solncoords(3,kkk,i_soln) = rtmp(kkk)
         enddo
         call twh_crmatmul(us,r_n,2,r0,rtmp)
         do kkk = 1,3
            solncoords(4,kkk,i_soln) = rtmp(kkk)
         enddo
         call twh_crmatmul(us,r_a,2,r0,rtmp)
         do kkk = 1,3
            solncoords(5,kkk,i_soln) = rtmp(kkk)
         enddo
         call twh_crmatmul(us,r_c,2,r0,rtmp)
         do kkk = 1,3
            solncoords(6,kkk,i_soln) = rtmp(kkk)
         enddo
         call twh_crmatmul(us,r_n,3,r0,rtmp)
         do kkk = 1,3
            solncoords(7,kkk,i_soln) = rtmp(kkk)
         enddo         
      end do
c     --- end of the coord from poly roots section

      return
      end 

      subroutine twh_test_two_cone_exist_soln(max_soln,tt,lckx,et,ap
     &     ,nsoln)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- parameters passed to/from the subroutine
      integer max_soln
c     --- variables passed to/from the subroutine
      double precision tt, lckx, et, ap
      integer nsoln
c     --- local variables
      character*2 case_type
      double precision at, ex, abs_at, ap1, kx1, et1,onepi
      double precision cos_tx1, cos_tx2, cos_te1, cos_te2, cos_ea1
      double precision cos_ea2, cos_xa1, cos_xa2, halfpi
      logical s1, s2, t1, t2

c     --- retrieve constants
      call twh_constant_halfpi(GLB_GET,halfpi)
      call twh_constant_pi(GLB_GET,onepi)
      nsoln = max_soln

      if (lckx .gt. halfpi) then
         case_type = '1'
      else if (et .gt. halfpi) then
         case_type = '2'
      else
         case_type = '0'
      end if

      ap1 = ap
      kx1 = lckx
      et1 = et
      if (case_type .ne. '0') then
         ap1 = onepi - ap1
         if (case_type .eq. '1') then
            kx1 = onepi - kx1
         else
            et1 = onepi - et1
         end if
      end if
  
      at = ap1 - tt
      ex = kx1 + et1
      abs_at = abs(at)

c     --- case of no soln
      if (abs_at .gt. ex) then
         nsoln = 0
         return
      end if

c     --- find type of intersection
      cos_tx1 = cos(tt+kx1)
      cos_tx2 = cos(tt-kx1)
      cos_te1 = cos(tt+et1)
      cos_te2 = cos(tt-et1)
      cos_ea1 = cos(et1+ap1)
      cos_ea2 = cos(et1-ap1)
      cos_xa1 = cos(kx1+ap1)
      cos_xa2 = cos(kx1-ap1)
      s1 = .false.
      s2 = .false.
      t1 = .false. 
      t2 = .false. 
      if ((cos_te1-cos_xa2)*(cos_te1-cos_xa1) .le. 0.0d0) s1 = .true.
      if ((cos_te2-cos_xa2)*(cos_te2-cos_xa1) .le. 0.0d0) s2 = .true.
      if ((cos_tx1-cos_ea2)*(cos_tx1-cos_ea1) .le. 0.0d0) t1 = .true.
      if ((cos_tx2-cos_ea2)*(cos_tx2-cos_ea1) .le. 0.0d0) t2 = .true.

      if (.not.s1 .and. .not.s2 .and. .not.t1 .and. .not.s2) then
         nsoln = 0
         return
      end if
      
      return

      end

      subroutine twh_get_poly_coeff(lfailure,poly_coeff,cos_alpha
     &     ,cos_theta,sin_alpha,cos_delta,sin_delta,cos_xi,cos_eta
     &     ,sin_xi,sin_eta,c0,c1,c2,crq,lcr)
c     ******************************************************************
c     * Solves for the coefficients of the 16th order polynomial that  *
c     * used to determine the possible loop closure solutions          *
c     *                                                                *
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      double precision poly_coeff
      dimension poly_coeff(0:16)
      double precision cos_alpha
      dimension cos_alpha(3)
      double precision cos_theta
      dimension cos_theta(3)
      double precision sin_alpha
      dimension sin_alpha(3)
      double precision cos_delta
      dimension cos_delta(0:3)
      double precision sin_delta
      dimension sin_delta(0:3)
      double precision cos_xi
      dimension cos_xi(3)
      double precision cos_eta
      dimension cos_eta(3)
      double precision sin_xi
      dimension sin_xi(3)
      double precision sin_eta
      dimension sin_eta(3)
      double precision c0
      dimension c0(0:2,3)
      double precision c1
      dimension c1(0:2,3)
      double precision c2
      dimension c2(0:2,3)
      double precision crq
      dimension crq(0:16,0:4)
      double precision lcr
      dimension lcr(0:16,0:2)
c     --- local variables
      integer i, j
      double precision a0,a1,a2,a3,a4,a21,a22,a31,a32,a41,a42
      double precision b0,b1,b2,b3,b4,b5,b6,b7,b8
      dimension b0(3),b1(3),b2(3),b3(3),b4(3),b5(3),b6(3),b7(3),b8(3)
      double precision u11,u12,u13,u31,u32,u33,um1,um2,um3,um4,um5,um6
      double precision q_tmp
      dimension u11(0:4,0:4),u12(0:4,0:4),u13(0:4,0:4),u31(0:4,0:4)
      dimension u32(0:4,0:4),u33(0:4,0:4),um1(0:4,0:4),um2(0:4,0:4)
      dimension um3(0:4,0:4),um4(0:4,0:4),um5(0:4,0:4),um6(0:4,0:4)
      dimension q_tmp(0:4,0:4)
      double precision f1,f2,f3,f4,f5,f6,f7,f8,f9,f10,f11,f12,f13,f14
      double precision f15,f16,f17,f18,f19,f20,f21,f22,f23,f24,f25,f26
      double precision f27,temp1,temp2,temp3,temp4,temp5,temp6,temp7
      double precision temp8
      dimension f1(0:16),f2(0:16),f3(0:16),f4(0:16),f5(0:16),f6(0:16)
      dimension f7(0:16),f8(0:16),f9(0:16),f10(0:16),f11(0:16),f12(0:16)
      dimension f13(0:16),f14(0:16),f15(0:16),f16(0:16),f17(0:16)
      dimension f18(0:16),f19(0:16),f20(0:16),f21(0:16),f22(0:16)
      dimension f23(0:16),f24(0:16),f25(0:16),f26(0:16),f27(0:16)
      dimension temp1(0:16),temp2(0:16),temp3(0:16),temp4(0:16)
      dimension temp5(0:16),temp6(0:16),temp7(0:16),temp8(0:16)
      integer p1,p3,p_um1,p_um2,p_um3,p_um4,p_um5,p_um6,p_q
      dimension p1(2),p3(2),p_um1(2),p_um2(2),p_um3(2),p_um4(2),p_um5(2)
      dimension p_um6(2),p_q(2)
      integer p2, p4, p_f1, p_f2, p_f3, p_f4, p_f5, p_f6, p_f7,p_f8,p_f9
      integer p_f10,p_f11,p_f12,p_f13,p_f14,p_f15,p_f16,p_f17,p_f18
      integer p_f20,p_f21,p_f22,p_f23,p_f24,p_f25,p_f26,p_final
      integer ii,jj,p_f19

c     --- Initialize arrays
      do i = 0,4
         do ii = 0,4
            u11(i,ii) = 0.0d0
            u12(i,ii) = 0.0d0
            u13(i,ii) = 0.0d0
            u31(i,ii) = 0.0d0
            u32(i,ii) = 0.0d0
            u33(i,ii) = 0.0d0
            um1(i,ii) = 0.0d0
            um2(i,ii) = 0.0d0
            um3(i,ii) = 0.0d0
            um4(i,ii) = 0.0d0
            um5(i,ii) = 0.0d0
            um6(i,ii) = 0.0d0
            q_tmp(i,ii) = 0.0d0
         enddo
      enddo

      do i = 0,16
         f1(i) = 0.0d0
         f2(i) = 0.0d0
         f3(i) = 0.0d0
         f4(i) = 0.0d0
         f5(i) = 0.0d0
         f6(i) = 0.0d0
         f7(i) = 0.0d0
         f8(i) = 0.0d0
         f9(i) = 0.0d0
         f19(i) = 0.0d0
         f11(i) = 0.0d0
         f12(i) = 0.0d0
         f13(i) = 0.0d0
         f14(i) = 0.0d0
         f15(i) = 0.0d0
         f16(i) = 0.0d0
         f17(i) = 0.0d0
         f18(i) = 0.0d0
         f19(i) = 0.0d0
         f20(i) = 0.0d0
         f21(i) = 0.0d0
         f22(i) = 0.0d0
         f23(i) = 0.0d0
         f24(i) = 0.0d0
         f25(i) = 0.0d0
         f26(i) = 0.0d0
         f27(i) = 0.0d0
      enddo

c     --- A0, B0
      do i = 1, 3
         a0 = cos_alpha(i)*cos_xi(i)*cos_eta(i) - cos_theta(i)
         a1 = -sin_alpha(i)*cos_xi(i)*sin_eta(i)
         a2 = sin_alpha(i)*sin_xi(i)*cos_eta(i)
         a3 = sin_xi(i)*sin_eta(i)
         a4 = a3*cos_alpha(i)
         j = i - 1
         a21 = a2*cos_delta(j)
         a22 = a2*sin_delta(j)
         a31 = a3*cos_delta(j)
         a32 = a3*sin_delta(j)
         a41 = a4*cos_delta(j)
         a42 = a4*sin_delta(j)
         b0(i) = a0 + a22 + a31
         b1(i) = 2.0d0*(a1 + a42)
         b2(i) = 2.0d0*(a32 - a21)
         b3(i) = -4.0d0*a41
         b4(i) = a0 + a22 - a31
         b5(i) = a0 - a22 - a31
         b6(i) = -2.0d0*(a21 + a32)
         b7(i) = 2.0d0*(a1 - a42)
         b8(i) = a0 - a22 + a31
      end do
c     --- c0i
      i = 1
      c0(0,i) = b0(i)
      c0(1,i) = b2(i)
      c0(2,i) = b5(i)
      c1(0,i) = b1(i)
      c1(1,i) = b3(i)
      c1(2,i) = b7(i)
      c2(0,i) = b4(i)
      c2(1,i) = b6(i)
      c2(2,i) = b8(i)
      do i = 2, 3
         c0(0,i) = b0(i)
         c0(1,i) = b1(i)
         c0(2,i) = b4(i)
         c1(0,i) = b2(i)
         c1(1,i) = b3(i)
         c1(2,i) = b6(i)
         c2(0,i) = b5(i)
         c2(1,i) = b7(i)
         c2(2,i) = b8(i)
      end do

c     --- first determinant
      do i = 0, 2
         u11(i,0) = c0(i,1)
         u12(i,0) = c1(i,1)
         u13(i,0) = c2(i,1)
         u31(0,i) = c0(i,2)
         u32(0,i) = c1(i,2)
         u33(0,i) = c2(i,2)
      end do

      p1(1) = 2
      p1(2) = 0
      p3(1) = 0
      p3(2) = 2


      call twh_poly_mul_sub2(u32,u32,u31,u33,p3,p3,p3,p3,um1,p_um1)

      call twh_poly_mul_sub2(u12,u32,u11,u33,p1,p3,p1,p3,um2,p_um2)

      call twh_poly_mul_sub2(u12,u33,u13,u32,p1,p3,p1,p3,um3,p_um3)

      call twh_poly_mul_sub2(u11,u33,u31,u13,p1,p3,p3,p1,um4,p_um4)

      call twh_poly_mul_sub2(u13,um1,u33,um2,p1,p_um1,p3,p_um2,um5,
     &                   p_um5)

      call twh_poly_mul_sub2(u13,um4,u12,um3,p1,p_um4,p1,p_um3,um6,
     &                   p_um6)

      call twh_poly_mul_sub2(u11,um5,u31,um6,p1,p_um5,p3,p_um6,
     &                   q_tmp,p_q)

      do ii = 0,4
         do jj = 0,4
            crq(ii,jj) = q_tmp(ii,jj)
         enddo
      enddo

c     --- second determinant
      do ii = 0,16
         do jj = 0,2
            lcr(ii,jj) = 0.0d0
         enddo
      enddo
      do ii = 0,2
         lcr(ii,0) = c0(ii,3)
         lcr(ii,1) = c1(ii,3)
         lcr(ii,2) = c2(ii,3)
      enddo
      p2 = 2
      p4 = 4

      do ii = 0,16
         temp1(ii) = lcr(ii,1) 
         temp2(ii) = lcr(ii,0) 
         temp3(ii) = lcr(ii,2)
         temp4(ii) = crq(ii,0)
         temp5(ii) = crq(ii,1)
         temp6(ii) = crq(ii,2)
         temp7(ii) = crq(ii,3)
         temp8(ii) = crq(ii,4)
      enddo
      call twh_poly_mul_sub1(temp1,temp1,temp2,temp3,p2,p2,p2,
     &                   p2,f1,p_f1)
      call twh_poly_mul1(temp1, temp3, p2, p2, f2, p_f2)

      call twh_poly_mul_sub1(temp1,f1,temp2,f2,p2,p_f1,p2,p_f2,f3,
     &                   p_f3)
      call twh_poly_mul1(temp3, f1, p2, p_f1, f4, p_f4)
      call twh_poly_mul_sub1(temp1,f3,temp2,f4,p2,p_f3,p2,p_f4,f5,
     &                   p_f5)

      call twh_poly_mul_sub1(temp5,temp1,temp4,temp3,p4,p2,
     &                   p4,p2,f6,p_f6)
      call twh_poly_mul_sub1(temp6,f1,temp3,f6,p4,p_f1,p2,p_f6,f7,
     &                   p_f7)
      call twh_poly_mul_sub1(temp7,f3,temp3,f7,p4,p_f3,p2,p_f7,f8,
     &                   p_f8)
      call twh_poly_mul_sub1(temp8,f5,temp3,f8,p4,p_f5,p2,p_f8,f9,
     &                   p_f9)

      call twh_poly_mul_sub1(temp7,temp1,temp8,temp2,p4,p2,
     &                   p4,p2,f10,p_f10)
      call twh_poly_mul_sub1(temp6,f1,temp2,f10,p4,p_f1,p2,
     &                   p_f10,f11,p_f11)
      call twh_poly_mul_sub1(temp5,f3,temp2,f11,p4,p_f3,p2,p_f11,
     &                   f12,p_f12)

      call twh_poly_mul_sub1(temp6,temp1,temp5,temp3,p4,p2,
     &                   p4,p2,f13,p_f13) 
      call twh_poly_mul_sub1(temp7,f1,temp3,f13,p4,p_f1,p2,p_f13,
     &                   f14,p_f14)
      call twh_poly_mul_sub1(temp7,temp1,temp6,temp3,p4,p2,
     &                   p4,p2,f15,p_f15)
      call twh_poly_mul_sub1(temp8,f1,temp3,f15,p4,p_f1,p2,p_f15,
     &                   f16,p_f16)
      call twh_poly_mul_sub1(temp5,f14,temp4,f16,p4,p_f14,p4,
     &                   p_f16,f17,p_f17)

      call twh_poly_mul_sub1(temp6,temp3,temp7,temp1,p4,p2,
     &                   p4,p2,f18,p_f18)
      call twh_poly_mul_sub1(temp5,temp3,temp7,temp2,p4,p2,
     &                   p4,p2,f19,p_f19)
      call twh_poly_mul_sub1(temp7,f19,temp6,f18,p4,p_f19,p4,
     &                   p_f18,f20,p_f20)
      call twh_poly_mul_sub1(temp5,temp1,temp6,temp2,p4,
     &                   p2,p4,p2,f21,p_f21)
      call twh_poly_mul1(temp8,f21,p4,p_f21,f22,p_f22)
      call twh_poly_sub1(f20,f22,p_f20,p_f22,f23,p_f23)
      call twh_poly_mul1(temp2,f23,p2,p_f23,f24,p_f24)
      call twh_poly_sub1(f17, f24, p_f17, p_f24, f25, p_f25)
      call twh_poly_mul_sub1(temp8,f12,temp3,f25,p4,p_f12,p2,
     &                   p_f25,f26,p_f26)
      call twh_poly_mul_sub1(temp4,f9,temp2,f26,p4,p_f9,p2,p_f26,
     &                   poly_coeff, p_final)

      do ii = 0,16
         lcr(ii,1) = temp1(ii) 
         lcr(ii,0) = temp2(ii)
         lcr(ii,2) = temp3(ii)
         crq(ii,0) = temp4(ii)
         crq(ii,1) = temp5(ii)
         crq(ii,2) = temp6(ii)
         crq(ii,3) = temp7(ii)
         crq(ii,4) = temp8(ii)
      enddo


      if (p_final .ne. 16) then
         write(6,*) 'GET_POLY_CLOSURE: Degree of polynomial is not 16'
         lfailure = .true.
         return
      endif

      if (poly_coeff(16) .lt. 0.0d0) then
         do ii = 0,16
            poly_coeff(ii) = -poly_coeff(ii)
         enddo 
      end if

c      print*, 'poly_coeff'
c      do i = 0, 16
c         write(*,"(i5,e15.6)") i, poly_coeff(16)
c      enddo

      return
      end 


      subroutine twh_poly_mul_sub2(u1,u2,u3,u4,p1,p2,p3,p4,u5,p5)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision u1, u2, u3, u4,u5
      dimension u1(0:4,0:4), u2(0:4,0:4), u3(0:4,0:4), u4(0:4,0:4)
      dimension u5(0:4,0:4)
      integer p1, p2, p3, p4,p5
      dimension p1(2),p2(2),p3(2),p4(2),p5(2)
c     --- local variables 
      double precision d1, d2
      dimension d1(0:4,0:4),d2(0:4,0:4)
      integer pd1, pd2
      dimension pd1(2),pd2(2)


      call twh_poly_mul2(u1, u2, p1, p2, d1, pd1)

      call twh_poly_mul2(u3, u4, p3, p4, d2, pd2)

      call twh_poly_sub2(d1, d2, pd1, pd2, u5, p5)
      
      return
      end


      subroutine twh_poly_mul2(u1, u2, p1, p2, u3, p3)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 03-02-2006 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision u1,u2,u3
      dimension u1(0:4,0:4),u2(0:4,0:4),u3(0:4,0:4)
      integer p1,p2,p3 
      dimension p1(2),p2(2),p3(2)

c     --- local variables
      integer i1,j1,i2,j2,i3,j3,p11,p12,p21,p22
      integer j,k
      double precision u1ij
      

      do j = 1,2
         p3(j) = p1(j) + p2(j)
      enddo
      do j = 0,4
         do k = 0,4
            u3(j,k) = 0.0d0
         enddo
      enddo

      p11 = p1(1)
      p12 = p1(2)
      p21 = p2(1)
      p22 = p2(2)

      do i1 = 0, p12
         do j1 = 0, p11
            u1ij = u1(j1,i1)
            do i2 = 0, p22
               i3 = i1 + i2
               do j2 = 0, p21
                  j3 = j1 + j2
                  u3(j3,i3) = u3(j3,i3) + u1ij*u2(j2,i2)
               end do
            end do
         end do
      end do
 

      
      return
      end

      subroutine twh_poly_sub2(u1, u2, p1, p2, u3, p3)
c     ******************************************************************
c     *                                                                *
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-08-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision u1,u2,u3
      dimension u1(0:4,0:4),u2(0:4,0:4),u3(0:4,0:4)
      integer p1,p2,p3
      dimension p1(2),p2(2),p3(2)

c     --- local variables
      integer i, j, p11, p12, p21, p22
      logical i1_ok, i2_ok

      p11 = p1(1)
      p12 = p1(2)
      p21 = p2(1)
      p22 = p2(2)
      p3(1) = max(p11,p21)
      p3(2) = max(p12,p22)

      do i = 0, p3(2)
         i1_ok = (i .gt. p12)
         i2_ok = (i .gt. p22)
         do j = 0, p3(1)
            if (i2_ok .or. (j .gt. p21)) then
               u3(j,i) = u1(j,i)
            else if (i1_ok .or. (j .gt. p11)) then
               u3(j,i) = -u2(j,i)
            else
               u3(j,i) = u1(j,i) - u2(j,i)
            end if
         end do
      end do
 
      return
      end


      subroutine twh_poly_mul_sub1(u1,u2,u3,u4,p1,p2,p3,p4,u5,p5)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision u1,u2,u3,u4,u5
      dimension u1(0:16),u2(0:16),u3(0:16),u4(0:16),u5(0:16)
      integer p1,p2,p3,p4,p5

c     --- local variables
      double precision d1,d2
      dimension d1(0:16),d2(0:16)
      integer pd1, pd2

      call twh_poly_mul1(u1, u2, p1, p2, d1, pd1)
      call twh_poly_mul1(u3, u4, p3, p4, d2, pd2)
      call twh_poly_sub1(d1, d2, pd1, pd2, u5, p5)
 
      return
      end

      subroutine twh_poly_mul1(u1, u2, p1, p2, u3, p3)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision u1,u2,u3
      dimension u1(0:16),u2(0:16),u3(0:16)
      integer p1,p2,p3
c     --- local variables
      integer i1, i2, i3,j
      double precision u1i

      p3 = p1 + p2
      do j = 0,16
         u3(j) = 0.0d0
      enddo

      do i1 = 0, p1
         u1i = u1(i1)
         do i2 = 0, p2
            i3 = i1 + i2
            u3(i3) = u3(i3) + u1i*u2(i2)
         end do
      end do

      return
      end


      subroutine twh_poly_sub1(u1, u2, p1, p2, u3, p3)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision u1,u2,u3
      dimension u1(0:16),u2(0:16),u3(0:16)
      integer p1,p2,p3

c     --- local variables
      integer i

      p3 = max(p1, p2)

      do i = 0, p3
         if (i .gt. p2) then
            u3(i) = u1(i)
         else if (i .gt. p1) then
            u3(i) = -u2(i)
         else
            u3(i) = u1(i) - u2(i)
         end if
      end do
 
      return
      end

      function twh_calc_t2(t0,crq,lcr)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision t0
      double precision crq
      dimension crq(0:16,0:4)
      double precision lcr
      dimension lcr(0:16,0:2)
c     --- local variables 
      double precision twh_calc_t2
      double precision b0,b1,b2,a0,a1,a2,a3,a4,b2_2,b2_3
      double precision k0,k1,k2,k3,t0_2,t0_3,t0_4

      t0_2 = t0*t0
      t0_3 = t0_2*t0
      t0_4 = t0_3*t0

      a0 = crq(0,0)+crq(1,0)*t0+crq(2,0)*t0_2+crq(3,0)*t0_3
     &             +crq(4,0)*t0_4
      a1 = crq(0,1)+crq(1,1)*t0+crq(2,1)*t0_2+crq(3,1)*t0_3
     &             +crq(4,1)*t0_4
      a2 = crq(0,2)+crq(1,2)*t0+crq(2,2)*t0_2+crq(3,2)*t0_3
     &             +crq(4,2)*t0_4
      a3 = crq(0,3)+crq(1,3)*t0+crq(2,3)*t0_2+crq(3,3)*t0_3
     &             +crq(4,3)*t0_4
      a4 = crq(0,4)+crq(1,4)*t0+crq(2,4)*t0_2+crq(3,4)*t0_3
     &             +crq(4,4)*t0_4

      b0 = lcr(0,0) + lcr(1,0)*t0 + lcr(2,0)*t0_2
      b1 = lcr(0,1) + lcr(1,1)*t0 + lcr(2,1)*t0_2
      b2 = lcr(0,2) + lcr(1,2)*t0 + lcr(2,2)*t0_2

      b2_2 = b2*b2
      b2_3 = b2_2*b2

      k0 = a2*b2 - a4*b0
      k1 = a3*b2 - a4*b1
      k2 = a1*b2_2 - k1*b0
      k3 = k0*b2 - k1*b1
  
      twh_calc_t2 = (k3*b0 - a0*b2_3)/(k2*b2 - k3*b1)

      return
      end


      function twh_calc_t1(t0,t2,c0,c1,c2)
c     ******************************************************************
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision  t0, t2
      double precision c0
      dimension c0(0:2,3)
      double precision c1
      dimension c1(0:2,3)
      double precision c2
      dimension c2(0:2,3)
c     --- local variables
      double precision twh_calc_t1
      double precision u11,u12,u13,u31,u32,u33
      double precision t0_2,t2_2

      t0_2 = t0*t0
      t2_2 = t2*t2

      u11 = c0(0,1) + c0(1,1)*t0 + c0(2,1)*t0_2
      u12 = c1(0,1) + c1(1,1)*t0 + c1(2,1)*t0_2
      u13 = c2(0,1) + c2(1,1)*t0 + c2(2,1)*t0_2
      u31 = c0(0,2) + c0(1,2)*t2 + c0(2,2)*t2_2
      u32 = c1(0,2) + c1(1,2)*t2 + c1(2,2)*t2_2
      u33 = c2(0,2) + c2(1,2)*t2 + c2(2,2)*t2_2

      twh_calc_t1 = (u31*u13-u11*u33)/(u12*u33-u13*u32)
  
      return
      end

      subroutine twh_quaternion(axis, quarter_ang, p)
c     ******************************************************************
c     * calculate quaternion, given rotation axis and angle.           *
c     *                                                                *
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision axis,quarter_ang,p
      dimension p(4),axis(3)

c     --- local variables
      double precision tan_w, tan_sqr, tan1, cosine, sine

      tan_w = tan(quarter_ang)
      tan_sqr = tan_w * tan_w
      tan1 = 1.0d0 + tan_sqr
      cosine = (1.0d0 - tan_sqr)/tan1
      sine = 2.0d0*tan_w/tan1
      p(1) = cosine
      p(2) = axis(1) * sine
      p(3) = axis(2) * sine
      p(4) = axis(3) * sine

      return
      end

      subroutine twh_rotation_matrix(qr, u)
c     ******************************************************************
c     * constructs rotation matrix u from quaternion qr.               *
c     *                                                                *
c     * originally rewritten from E. Coutsias 05-01-2003 by C.D. Lorenz*
c     * last modified 10-09-2005 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision qr,u
      dimension qr(4),u(3,3)
c     --- local variables
      double precision q0,q1,q2,q3,b0,b1,b2,b3,q00,q01,q02,q03,q11,q12
      double precision q13,q22,q23,q33
  
      q0 = qr(1)
      q1 = qr(2) 
      q2 = qr(3)
      q3 = qr(4)
      b0 = 2.0d0*q0
      b1 = 2.0d0*q1
      q00 = b0*q0-1.0d0
      q02 = b0*q2
      q03 = b0*q3
      q11 = b1*q1
      q12 = b1*q2
      q13 = b1*q3  
      b2 = 2.0d0*q2
      b3 = 2.0d0*q3
      q01 = b0*q1
      q22 = b2*q2
      q23 = b2*q3
      q33 = b3*q3 
      u(1,1) = q00+q11
      u(1,2) = q12-q03
      u(1,3) = q13+q02
      u(2,1) = q12+q03
      u(2,2) = q00+q22
      u(2,3) = q23-q01
      u(3,1) = q13-q02
      u(3,2) = q23+q01
      u(3,3) = q00+q33

      return
      end

      subroutine twh_crmatmul(r1, r2, r2flag, r3, r4)
c     ******************************************************************
c     * matrix multiplication for a (3x3) matrix * a vector (3x1)      *
c     * r1 = (3x3) matrix, r2 = (3x3) matrix, r3 = input vector        *
c     * r4 = resulting vector                                          *
c     * first takes difference of the (r2flag)th column of the 3x3     *
c     * matrix and r3 and then does multiplication and then adds r3    *
c     * back in and returns r4 as the solution                         *
c     *                                                                *
c     * originally written 05-01-2003 by C.D. Lorenz                   *
c     * last modified 10-09-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      double precision r1,r2,r3,r4
      dimension r1(3,3), r2(3,3),r3(3),r4(3)

      integer r2flag

c     --- local variable
      integer i
      double precision r2t
      dimension r2t(3)

      if (r2flag .eq. 0) then
         do i = 1,3
            r2t(i) = r2(i,1)
         enddo
         r4(1) = (r1(1,1)*r2t(1))+(r1(1,2)*r2t(2))+(r1(1,3)*r2t(3))
         r4(2) = (r1(2,1)*r2t(1))+(r1(2,2)*r2t(2))+(r1(2,3)*r2t(3))
         r4(3) = (r1(3,1)*r2t(1))+(r1(3,2)*r2t(2))+(r1(3,3)*r2t(3))

      else 

         do i = 1,3
            r2t(i) = r2(i,r2flag) - r3(i) 
         enddo

         r4(1) = (r1(1,1)*r2t(1))+(r1(1,2)*r2t(2))+(r1(1,3)*r2t(3))
         r4(2) = (r1(2,1)*r2t(1))+(r1(2,2)*r2t(2))+(r1(2,3)*r2t(3))
         r4(3) = (r1(3,1)*r2t(1))+(r1(3,2)*r2t(2))+(r1(3,3)*r2t(3))

         do i = 1,3
            r4(i) = r4(i) + r3(i)
         enddo
      endif
      
      return
      end

