#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2010 Marcus G. Martin                       *
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
      function twh_vimproper(ibox,imolty,unit_1,impind,xcord
     &     ,ycord,zcord,phic,lbadstereo)
c     ******************************************************************
c     * Computes the improper torsional energy for a particular        *
c     * improper type.                                                 *
c     * ibox: the box the molecule is in - needed for minimum image    *
c     * imolty: molecule type                                          *
c     * unit_1: unit numbers of the starting atom of the torsion       *
c     * impind: the index of the torsion                               *
c     * x,y,z cord : coordinates of the four atoms in the torsion      *
c     * phic: cosine of the torsional angle (oop impropers)            *
c     * phic: is instead the angle from -pi to +pi for stereo impropers*
c     *                                                                *
c     * rewritten from vtorsion 03-21-2002 by M.G. Martin              *
c     * last modified 11-16-2010 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_DOTPRODUCT
#define FUNCTION_GET_IMPCOEFF
#define FUNCTION_GET_IMPFORM
#define FUNCTION_GET_IMPSTYLE
#define FUNCTION_GET_ITIMPROP
#define FUNCTION_VIMPROPER
#include "functions.h"
c     --- variables passed to and from function
      logical lbadstereo
      integer imolty,impind,unit_1,ibox
      double precision phic
      double precision xcord,ycord,zcord
      dimension xcord(4),ycord(4),zcord(4)
c     --- local variables
c     --- integer scalars
      integer style,ivec,improperform,impffnum
c     --- double precision scalars
      double precision phi,onepi,temperature
      double precision xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2
      double precision xaa3,yaa3,zaa3,daa3
      double precision xcc,ycc,zcc,tcc
      double precision dot1,dot2,dot3,halfpi
c     --- double precision arrays
      double precision xvec,yvec,zvec
      dimension xvec(4),yvec(4),zvec(4)
      double precision costheta,invstheta
      dimension costheta(3),invstheta(3)
      double precision phiarray
      dimension phiarray(3)

c      write(6,*) 'Start VIMPROPER'

c     --- retrieve constants
      call twh_constant_halfpi(GLB_GET,halfpi)

c     --- initialize lbadstereo and vimproper
      lbadstereo = .false.
      twh_vimproper = 0.0d0
      phi = 0.0d0
      do ivec = 1,3
         phiarray(ivec) = 0.0d0
      enddo

c     --- determine the improper forcefield number and style
      impffnum = twh_get_itimprop(imolty,unit_1,impind)
c     --- determine the improper form
      improperform = twh_get_impform(impffnum)
c     --- determine the improper style
      style = twh_get_impstyle(impffnum)
      if ( improperform .eq. 1 ) then
c        --- amber/stereo improper angle
         do ivec = 1,3
            call twh_mimage(ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(1),ycord(1),zcord(1)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- we need a fourth vector for the scalar triple product
         call twh_mimage(ibox,xcord(1),ycord(1),zcord(1),xcord(2)
     &        ,ycord(2),zcord(2),xvec(4),yvec(4),zvec(4))
c        --- cross product i-k x k-l (3-1 x 1-2)
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
c        --- cross product j-k x k-l (4-1 x 1-2)
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(3),yvec(3)
     &        ,zvec(3),xaa2,yaa2,zaa2)
c        --- calculate dot product of cross products
         phic = twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2)
         phi = twh_arccos(phic)
c        --- calculate cross product of cross products
         call twh_crossproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,xcc,ycc
     &        ,zcc)
c        --- calculate scalar triple product***
         tcc = xcc*(xvec(4)) + ycc*(yvec(4)) + zcc*(zvec(4))
c        --- compute the angle over a range -pi to +pi
         if ( tcc .lt. 0.0d0 ) phi = -phi
      elseif ( improperform .eq. 2 ) then
c        --- compute the standard oop improper
c        --- defined so that a cis bond is at 0.0 phi
         do ivec = 1,3
            call twh_mimage(ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(ivec),ycord(ivec),zcord(ivec)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- calculate cross products 1-2 x 2-3 
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
c        --- calculate cross products 3-4 x 2-3
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(2),yvec(2)
     &        ,zvec(2),xaa2,yaa2,zaa2)
c        --- calculate dot product of cross products
         phic = -twh_dotproduct(xaa1,yaa1,zaa1,xaa2,yaa2,zaa2,daa1,daa2)
      elseif ( improperform .eq. 3 ) then
c        --- wilson oop, average of the three possible oop terms
c        --- defined so that a cis bond is at 0.0 phi
c        --- compute the vectors to the center
         do ivec = 1,3
            call twh_mimage(ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(1),ycord(1),zcord(1)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- compute angles
         costheta(1) = twh_dotproduct(xvec(1),yvec(1),zvec(1),xvec(2)
     &        ,yvec(2),zvec(2),daa1,daa2)
         costheta(2) = twh_dotproduct(xvec(2),yvec(2),zvec(2),xvec(3)
     &        ,yvec(3),zvec(3),daa2,daa3)
         costheta(3) = twh_dotproduct(xvec(1),yvec(1),zvec(1),xvec(3)
     &        ,yvec(3),zvec(3),daa1,daa3)
c        --- compute inverse sin of those angles
         do ivec = 1,3
            invstheta(ivec) = 1.0d0/dsin(twh_arccos(costheta(ivec)))
c           --- divide invstheta by the product of the three distances
            invstheta(ivec) = invstheta(ivec)/(daa1*daa2*daa3)
         enddo
c        --- cross products
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(1),yvec(1)
     &        ,zvec(1),xaa2,yaa2,zaa2)
         call twh_crossproduct(xvec(2),yvec(2),zvec(2),xvec(3),yvec(3)
     &        ,zvec(3),xaa3,yaa3,zaa3)
c        --- calculate dot product of cross products with original vecs
c        --- want the unscaled dotproducts
         dot1 = twh_dotproduct(xaa3,yaa3,zaa3,xvec(1),yvec(1),zvec(1)
     &        ,daa1,daa2)
         dot1 = dot1*daa1*daa2
         dot2 = twh_dotproduct(xaa2,yaa2,zaa2,xvec(2),yvec(2),zvec(2)
     &        ,daa1,daa2)
         dot2 = dot2*daa1*daa2
         dot3 = twh_dotproduct(xaa1,yaa1,zaa1,xvec(3),yvec(3),zvec(3)
     &        ,daa1,daa2)
         dot3 = dot3*daa1*daa2
c        --- average phi over the three values
         phi = (asin(dot1*invstheta(1))
     &        + asin(dot2*invstheta(2))
     &        + asin(dot3*invstheta(3)))/3.0d0
      elseif ( improperform .eq. 4 ) then
c        --- UFF out of plane definition
c        --- All atoms (j,k,l) are bonded to i
c        --- want the angle that IL axis makes with the IJK plane
c        --- get this from w = y - Pi
c        --- where y is the angle between IL and the normal of the plane

c        --- defined so that all atoms in the same plane is zero angle
         do ivec = 2,4
            call twh_mimage(ibox,xcord(ivec),ycord(ivec),zcord(ivec)
     &           ,xcord(1),ycord(1),zcord(1)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- take the cross product of IJ (2) with IK (3) to get the 
c        --- normal (1)
         call twh_crossproduct(xvec(2),yvec(2),zvec(2),xvec(3),yvec(3)
     &        ,zvec(3),xvec(1),yvec(1),zvec(1))
c        --- normal(1) dot IL(4) gives cosine of the angle with the 
c        --- normal
         phic = twh_dotproduct(xvec(1),yvec(1),zvec(1),xvec(4),yvec(4)
     &        ,zvec(4),daa1,daa2)
         phi = twh_arccos(phic)
c        --- subtract pi/2 to get angle of IL with the IJK plane
         phi = phi - halfpi
      elseif ( improperform .eq. 5 ) then
c        --- wilson oop, vector of angles
c        --- defined so that a cis bond is at 0.0 phi
c        --- compute the vectors to the center
         do ivec = 1,3
            call twh_mimage(ibox,xcord(ivec+1),ycord(ivec+1)
     &           ,zcord(ivec+1),xcord(1),ycord(1),zcord(1)
     &           ,xvec(ivec),yvec(ivec),zvec(ivec))
         enddo
c        --- compute angles
         costheta(1) = twh_dotproduct(xvec(1),yvec(1),zvec(1),xvec(2)
     &        ,yvec(2),zvec(2),daa1,daa2)
         costheta(2) = twh_dotproduct(xvec(2),yvec(2),zvec(2),xvec(3)
     &        ,yvec(3),zvec(3),daa2,daa3)
         costheta(3) = twh_dotproduct(xvec(1),yvec(1),zvec(1),xvec(3)
     &        ,yvec(3),zvec(3),daa1,daa3)
c        --- compute inverse sin of those angles
         do ivec = 1,3
            invstheta(ivec) = 1.0d0/dsin(twh_arccos(costheta(ivec)))
c           --- divide invstheta by the product of the three distances
            invstheta(ivec) = invstheta(ivec)/(daa1*daa2*daa3)
         enddo
c        --- cross products
         call twh_crossproduct(xvec(1),yvec(1),zvec(1),xvec(2),yvec(2)
     &        ,zvec(2),xaa1,yaa1,zaa1)
         call twh_crossproduct(xvec(3),yvec(3),zvec(3),xvec(1),yvec(1)
     &        ,zvec(1),xaa2,yaa2,zaa2)
         call twh_crossproduct(xvec(2),yvec(2),zvec(2),xvec(3),yvec(3)
     &        ,zvec(3),xaa3,yaa3,zaa3)
c        --- calculate dot product of cross products with original vecs
c        --- want the unscaled dotproducts
         dot1 = twh_dotproduct(xaa3,yaa3,zaa3,xvec(1),yvec(1),zvec(1)
     &        ,daa1,daa2)
         dot1 = dot1*daa1*daa2
         dot2 = twh_dotproduct(xaa2,yaa2,zaa2,xvec(2),yvec(2),zvec(2)
     &        ,daa1,daa2)
         dot2 = dot2*daa1*daa2
         dot3 = twh_dotproduct(xaa1,yaa1,zaa1,xvec(3),yvec(3),zvec(3)
     &        ,daa1,daa2)
         dot3 = dot3*daa1*daa2
c        --- compute the three values
         phiarray(1) = asin(dot1*invstheta(1))
         phiarray(2) = asin(dot2*invstheta(2))
         phiarray(3) = asin(dot3*invstheta(3))
      else
         write(6,*) 'VIMPROPER: invalid improperform',improperform
         stop
      endif

      if ( style .eq. 1 ) then
c        --- simple harmonic potential 
c        --- used for TraPPE double bonds (cis vs trans) 
c        --- used for Gromos improper torsions
c        --- used for Charmm improper torsions
         phi = twh_arccos(phic)
         twh_vimproper = twh_get_impcoeff(impffnum,0)
     &        *(phi-twh_get_impcoeff(impffnum,1))**2.0d0
      elseif ( style .eq. 2  ) then
c        --- OPLS style cosine series 
c        --- used in the OPLS-AA force field
c        --- stero version so phi is phic
         twh_vimproper = twh_get_impcoeff(impffnum,1)*(1.d0 + dcos(phi))
     &        + twh_get_impcoeff(impffnum,2)*(1.d0 - dcos(2.0d0*phi))
     &        + twh_get_impcoeff(impffnum,3)*(1.d0 + dcos(3.0d0*phi))
      elseif ( style .eq. 3 ) then
c        --- potentials for enforcing stereochemistry
         if ( phi .lt. 0.0d0 ) then
            lbadstereo = .true.
            if ( phi .lt. -halfpi ) then
               call twh_constant_pi(GLB_GET,onepi)
               tcc = onepi + phi
            else
               tcc = -phi
            endif
            call twh_temperature(GLB_GET,temperature)
            twh_vimproper = (twh_get_impcoeff(impffnum,0) 
     &           + twh_get_impcoeff(impffnum,1)*tcc)*temperature
         endif
      elseif ( style .eq. 4 ) then
c        --- Amber improper torsion
         twh_vimproper = twh_get_impcoeff(impffnum,1)* 
     &        ( 1.0d0 + dcos(twh_get_impcoeff(impffnum,2)*phi
     &        - twh_get_impcoeff(impffnum,3)))
      elseif ( style .eq. 5 ) then
c        --- wilson harmonic
         twh_vimproper = twh_get_impcoeff(impffnum,0)
     &        *(phi-twh_get_impcoeff(impffnum,1))**2.0d0
      elseif ( style .eq. 6 ) then
c        --- UFF cosine Fourier expansion
         twh_vimproper = twh_get_impcoeff(impffnum,3)
     &        *(twh_get_impcoeff(impffnum,0) 
     &        + twh_get_impcoeff(impffnum,1)*dcos(phi)
     &        + twh_get_impcoeff(impffnum,2)*dcos(2.0d0*phi) )
      elseif ( style .eq. 7 ) then
c        --- k[1 - Cos(phi)]
         twh_vimproper = twh_get_impcoeff(impffnum,1)
     &        *(1.0d0 - dcos(phi))
      elseif ( style .eq. 8 ) then
c        --- harmonic over the three phi angles
         twh_vimproper = twh_get_impcoeff(impffnum,1)*(phiarray(1)**2
     &        + phiarray(2)**2 + phiarray(3)**2)
      else
         stop 'you picked a non-defined improper type'
      endif

c      write(6,*) 'finished VIMPROPER'

      return
      end
