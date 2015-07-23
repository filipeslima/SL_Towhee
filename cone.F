#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 1999 Bin Chen, Marcus G. Martin,                 *
c     * J. Ilja Siepmann, John Stubbs, and Collin D. Wick              *
c     * Copyright (C) 2000-2011 Marcus G. Martin                       *
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
      subroutine twh_cone(lfailure,iinit,xxx,yyy,zzz,alpha,gamma
     &,ux,uy,uz)
c     ******************************************************************
c     * if iinit = 1 then it sets up the rotation matrix for the cone  *
c     * using xxx,yyy,zzz as a unit vector pointing in the +z direction*
c     ******************************************************************
c     * if iinit = 2 then it creates a unit vector that has an angle   *
c     * of alpha from the +z direction (previous vector) and an angle  *
c     * of gamma (0,2Pi) around the cone circle and returns this as    *
c     * ux,uy,uz                                                       *
c     ******************************************************************
c     * if iinit = 3 then this computes gamma (0,2Pi) for a unit       *
c     * vector ux,uy,uz with angle alpha from the -z direction         *
c     * (Seems gamma is on [-Pi, Pi] -- confirm and change docs)       *
c     ******************************************************************
c     * note that the rotation matrix is saved after each call so it   *
c     * needs to be reset when you wish to use another cone            *
c     *                                                                *
c     * originally written prior to 1995                               *
c     * last modified 08-09-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_CMP_EQ
#include "functions.h"

c     --- variables passed to/from the subroutine
      logical lfailure
      integer iinit
      double precision xxx,yyy,zzz,alpha,gamma,ux,uy,uz

c     --- local variables
      double precision a11,a12,a13,a21,a22,a31,a32,a33
      double precision sinthe,costhe,sinpsi,cospsi,singamma,cosgamma
      double precision uxtemp,uytemp,uztemp,sinalph,cosalph
      double precision determ
      save a11,a12,a13,a21,a22,a31,a32,a33,cosalph,sinalph

c      write(6,*) 'start CONE'

      if ( iinit .eq. 1 ) then      
c       --- setup the unit cone
c       --- assume x,y,z is a rotation of vector (0,0,1)
c       --- conversions from xyz to theta,psi
c           xxx = -dsin(theta)*dcos(psi)
c           yyy = dsin(theta)*dsin(psi)
c           zzz = dcos(theta)
c       --- rotation matrix (note that sin(phi) = 0, cos(phi) = 1)
c       a11 = cos(psi) cos(theta) cos(phi) - sin(psi) sin(phi)
c           = cos(psi) cos(theta)
c       a12 = -sin(psi) cos(theta) cos(phi) - cos(psi) sin(phi)
c           = -sin(psi) cos(theta)
c       a13 = sin(theta) cos(phi) 
c           = sin(theta)
c       a21 = cos(psi) cos(theta) sin(phi) + sin (psi) cos(phi)
c           = sin(psi)
c       a22 = -sin(psi) cos(theta) sin(phi) + cos (psi) cos(phi)
c           = cos(psi)
c       a23 = sin(theta) sin(phi) 
c           = 0.0
c       a31 = -cos(psi) sin(theta)
c           = xxx
c       a32 = sin(psi)  sin (theta)
c           = yyy
c       a33 = cos(theta)
c           = zzz
         costhe = zzz
         sinthe = dsqrt(1.0d0 - costhe*costhe)
         if (sinthe .ne. 0.0d0) then
            sinpsi=(yyy)/sinthe
            cospsi=(-xxx)/sinthe
         else
            sinpsi = 0.0d0
            cospsi = 1.0d0
         endif
c        --- rotation matrix
         a11 = cospsi*costhe
         a12 = -(sinpsi*costhe)
         a13 = sinthe
         a21 = sinpsi
         a22 = cospsi
         a31 = xxx
         a32 = yyy
         a33 = zzz
      elseif ( iinit .eq. 2 ) then
c        --- find the vector on a cone:
c        --- alpha is the angle with the cone axis (theta)
         sinalph = dsin(alpha)
         cosalph = dcos(alpha)
c        --- gamma is now passed to cone - (phi)
         uztemp = -cosalph
         uytemp = sinalph*dcos(gamma)
         uxtemp = sinalph*dsin(gamma)
         ux = a11*uxtemp + a21*uytemp + a31*uztemp
         uy = a12*uxtemp + a22*uytemp + a32*uztemp
         uz = a13*uxtemp +              a33*uztemp
      elseif ( iinit .eq. 3 ) then
c        --- find gamma for a given unit vector on the cone
c        --- use cramer's rule where a23 has already been replaced by 0
c        --- we don't need the uztemp so it is not computed
         determ = 
     &        ( a11*a22*a33
     &        + a21*(a32*a13 - a12*a33)
     &        - a31*a22*a13 )

         uxtemp =  
     &        ( ux*a22*a33
     &        + a21*(a32*uz - uy*a33)
     &        - a31*a22*uz ) / determ

         uytemp = 
     &        ( a11*(uy*a33 - a32*uz)
     &        + ux*(a32*a13 - a12*a33)
     &        + a31*(a12*uz - uy*a13) ) / determ

#ifdef SAFE_COMPARE
c       if u*temp are not numerically different from zero, set them to 
c       zero (small errors change singamma, hence sign of gamma) 
        if (twh_cmp_eq(uxtemp, 0.0d0)) uxtemp = 0.0d0
        if (twh_cmp_eq(uytemp, 0.0d0)) uytemp = 0.0d0
#endif
         sinalph = dsin(alpha)
         singamma = uxtemp/sinalph
         cosgamma = uytemp/sinalph
c        --- now need to find the gamma on [-Pi,Pi] that satisfies cos 
c        --- and sine
         gamma = twh_arccos(cosgamma)
         if ( singamma .lt. 0.0d0 ) gamma = -gamma

      else
         write(6,*) 'iinit ',iinit
         write(6,*) 'CONE: non valid iinit'
         lfailure = .true.
         return
      endif
      
c      write(6,*) 'finish CONE'

      return
      end


