#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2001-2013 Marcus G. Martin                       *
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
      subroutine twh_tail()
c     ******************************************************************
c     * this subroutine sets up all of the constants needed to compute *
c     * the van der waals tail corrections for any two molecule types  *
c     * the constants are stored in the two arrays.                    *
c     * See Allen and Tildesley section 2.8 for more information about *
c     * long-range corrections                                         *
c     * coru: constants for vdw energy tail corrections                *
c     *       these constants do not include the N*rho prefactor where *
c     *       N is the number of molecules, but they do include all    *
c     *       cross terms for the different atom types in the two      *
c     *       molecule types                                           *
c     * coru = 2 Pi Integral(rcut,infinity) r*r*v(r) dr                *
c     * corp: constants for vdw pressure tail corrections              *
c     * corp = -(2/3) Pi Integral(rcut,infinity) r*r*w(r)dr            *
c     *      where w(r) = r dv(r)/dr                                   *
c     *                                                                *
c     * originally written 02-13-2001 by M.G. Martin                   *
c     * last modified 07-19-2013 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_HBONDCOEFF
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- local variables
c     --- character string scalars
      character*30 localcp
c     --- integer scalars
      integer imolty,jmolty,ii,jj,ntii,ntjj,ntij,hbtype
      integer nmolty,numiunit,numjunit
c     --- double precision scalars
      double precision rci3,overra,onepi,aaa,bbb,twopi,dvalue
      double precision rcut,sigma,powrn,aofn,fecofn

c     --- get the classical potential
      localcp = twh_get_classical_potential()
c     --- retrieve constants
      call twh_constant_pi(GLB_GET,onepi)
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_nmolty(GLB_GET,nmolty)
      call twh_rcut(GLB_GET,rcut)

      do imolty = 1,nmolty
         call twh_nunit(GLB_GET,imolty,numiunit)
         do jmolty = 1,nmolty
            call twh_nunit(GLB_GET,jmolty,numjunit)
c           --- initialize coru and corp
            call twh_coru(GLB_SET,imolty,jmolty,0.0d0)
            call twh_corp(GLB_SET,imolty,jmolty,0.0d0)

c           --- compute coru and corp
            do ii = 1, numiunit
               call twh_ntype(GLB_GET,imolty,ii,ntii)
               do jj = 1, numjunit 
                  call twh_ntype(GLB_GET,jmolty,jj,ntjj)
                  ntij = twh_getnbtype(ntii,ntjj)
                  if ( localcp .eq. 'Lennard-Jones'
     &                 .or. localcp .eq. 'UFF 12-6'
     &                 .or. localcp .eq. '12-6 plus solvation' 
     &                 .or. localcp .eq. 'Drukker H-bond'
     &                 ) then
c                    --- Lennard-Jones 12-6 tail correction
c                    --- Note: twh_get_nbcoeff(13,ij) is sigma^6
c                    --- Note: twh_get_nbcoeff(14,ij) is 4*epsilon
                     rci3 = dsqrt(twh_get_nbcoeff(13,ntij)) / rcut**3
                     dvalue = 2.0d0* onepi*twh_get_nbcoeff(14,ntij)* 
     &                    dsqrt(twh_get_nbcoeff(13,ntij))* 
     &                    (rci3* rci3* rci3 / 9.0d0 - rci3 / 3.0d0)
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                     dvalue = 8.0d0* onepi*twh_get_nbcoeff(14,ntij) 
     &                    * dsqrt(twh_get_nbcoeff(13,ntij))*
     &                    ( rci3*rci3*rci3/9.0d0 - rci3/6.0d0)
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                  elseif ( localcp 
     &                    .eq. 'Scaled Lennard-Jones' ) then
c                   --- analytical tail correction not implemented for
c                   --- Scaled Lennard Jones (soft core, in particular)
                        
                     write(6,*) 'Analytical tail correction (ltailc = 
     & .true) not implemented for Scaled Lennard-Jones potential.'
                     stop
                  elseif ( localcp 
     &                    .eq. '12-6 plus 12-10 H-bond' ) then
c                    --- 12-6 as above plus a 12-10
c                    --- Note: twh_get_nbcoeff(13,ij) is sigma^6
c                    --- Note: twh_get_nbcoeff(14,ij) is 4*epsilon
                     rci3 = dsqrt(twh_get_nbcoeff(13,ntij)) / rcut**3
c                    --- 12-6 LJ portion
                     dvalue = 2.0d0* onepi*twh_get_nbcoeff(14,ntij)* 
     &                    dsqrt(twh_get_nbcoeff(13,ntij))* 
     &                    (rci3* rci3* rci3 / 9.0d0 - rci3 / 3.0d0)
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                     dvalue = 8.0d0* onepi*twh_get_nbcoeff(14,ntij) 
     &                    * dsqrt(twh_get_nbcoeff(13,ntij))*
     &                    ( rci3*rci3*rci3/9.0d0 - rci3/6.0d0)
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
c                    --- 12-10 H-bond portion (if needed)
                     call twh_nb_to_hbtype(GLB_GET,ntij,hbtype)
                     if ( hbtype .ne. 0 ) then
                        dvalue = twopi * (twh_get_hbondcoeff(hbtype,1)
     &                       / (9.0d0 * rcut**9)
     &                       - twh_get_hbondcoeff(hbtype,2)
     &                       /(7.0d0 * rcut**7) )
                        call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                        dvalue = - (twopi/3.0d0) * ( -(12.0d0/9.0d0)
     &                       *twh_get_hbondcoeff(hbtype,1)
     &                       /(rcut**9)
     &                       +(10.0d0/7.0d0)
     &                       *twh_get_hbondcoeff(hbtype,2)/(rcut**7) )
                        call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                     endif
                  elseif ( localcp .eq. '9-6' ) then
c                    --- Lennard-Jones 9-6 tail correction
c                    --- Note: twh_get_nbcoeff(1,ij) is sigma^3
                     rci3 = twh_get_nbcoeff(1,ntij) / rcut**3
                     dvalue = 2.0d0* onepi*twh_get_nbcoeff(2,ntij)* 
     &                    twh_get_nbcoeff(1,ntij)* ((rci3* rci3)/3.0d0
     &                    - rci3)
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                     dvalue = onepi*twh_get_nbcoeff(2,ntij)
     &                    *twh_get_nbcoeff(1,ntij)
     &                    *( 2.0d0*rci3*rci3 - 4.0d0*rci3 )
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                  elseif ( localcp .eq. 'Exponential-6' 
     &                    ) then
c                    --- exponential-6 tail correction
                     dvalue = 2.0d0*onepi*(twh_get_nbcoeff(1,ntij)
     &                    /(3.0d0*(rcut**3))
     &                    +(-rcut*rcut+2.0d0*rcut
     &                    /twh_get_nbcoeff(3,ntij)
     &                    -2.0d0/(twh_get_nbcoeff(3,ntij)
     &                    *twh_get_nbcoeff(3,ntij)))
     &                    *twh_get_nbcoeff(2,ntij)
     &                    *twh_expon(twh_get_nbcoeff(3,ntij)
     &                    *rcut)/twh_get_nbcoeff(3,ntij))
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                     dvalue = (2.0d0/3.0d0)*onepi*(2.0d0
     &                    *twh_get_nbcoeff(1,ntij)
     &                    /(rcut*rcut*rcut)+twh_get_nbcoeff(2,ntij)
     &                    *twh_expon(twh_get_nbcoeff(3,ntij)*rcut)
     &                    *(-6.0d0/(twh_get_nbcoeff(3,ntij)
     &                    *twh_get_nbcoeff(3,ntij)
     &                    *twh_get_nbcoeff(3,ntij))+6.0d0*rcut
     &                    /(twh_get_nbcoeff(3,ntij)
     &                    *twh_get_nbcoeff(3,ntij))
     &                    -3.0d0*rcut*rcut/twh_get_nbcoeff(3,ntij)
     &                    +rcut*rcut*rcut))
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                  elseif ( localcp .eq. 'Exponential-12-6' 
     &                    ) then
c                    --- exp-6 plus LJ 12
                     dvalue = 2.0d0*onepi*(twh_get_nbcoeff(1,ntij)
     &                    /(3.0d0*(rcut**3))
     &                    +(-rcut*rcut+2.0d0*rcut
     &                    /twh_get_nbcoeff(4,ntij)
     &                    -2.0d0/(twh_get_nbcoeff(4,ntij)
     &                    *twh_get_nbcoeff(4,ntij)))
     &                    *twh_get_nbcoeff(3,ntij)
     &                    *twh_expon(twh_get_nbcoeff(4,ntij)
     &                    *rcut)/twh_get_nbcoeff(4,ntij))
     &                    + twopi*twh_get_nbcoeff(2,ntij)
     &                    /(9.0d0*(rcut**9))
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                     dvalue = (2.0d0/3.0d0)*onepi*(2.0d0
     &                    *twh_get_nbcoeff(1,ntij)
     &                    /(rcut*rcut*rcut)+twh_get_nbcoeff(3,ntij)
     &                    *twh_expon(twh_get_nbcoeff(4,ntij)*rcut)
     &                    *(-6.0d0/(twh_get_nbcoeff(4,ntij)
     &                    *twh_get_nbcoeff(4,ntij)
     &                    *twh_get_nbcoeff(4,ntij))+6.0d0*rcut
     &                    /(twh_get_nbcoeff(4,ntij)
     &                    *twh_get_nbcoeff(4,ntij))
     &                    -3.0d0*rcut*rcut/twh_get_nbcoeff(4,ntij)
     &                    +rcut*rcut*rcut))
     &                    + 8.0d0*onepi/(9.0d0*(rcut**9))
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                  elseif ( localcp .eq. '12-9-6' ) then
                     dvalue = twopi*(
     &                    twh_get_nbcoeff(1,ntij)/(9.0d0*(rcut**9))
     &                    + twh_get_nbcoeff(2,ntij)/(6.0d0*(rcut**6))
     &                    + twh_get_nbcoeff(3,ntij)/(3.0d0*(rcut**3)))
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                     dvalue = (2.0d0/3.0d0)*onepi*
     &                    ( (4.0d0*twh_get_nbcoeff(1,ntij))
     &                    /(3.0d0*rcut**9)
     &                    + (3.0d0*twh_get_nbcoeff(2,ntij))
     &                    /(2.0d0*rcut**6)
     &                    + (2.0d0*twh_get_nbcoeff(3,ntij))/(rcut**3) )
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                  elseif ( localcp .eq. 'Gordon n-6' ) then
                     call twh_nbcoeff(GLB_GET,1,ntij,sigma)
                     call twh_nbcoeff(GLB_GET,3,ntij,powrn)
                     call twh_nbcoeff(GLB_GET,4,ntij,aofn)
                     call twh_nbcoeff(GLB_GET,5,ntij,fecofn)
                     overra = 1.0d0/(rcut - aofn)
c                    --- from Peter Gordon 11-29-2004
                     dvalue = twopi*fecofn*((sigma**powrn)*(
     &                    (1.0d0/(powrn-3.0d0))*(overra**(powrn-3.0d0))
     &                    + ( (2.0d0*aofn)/(powrn-2.0d0) )
     &                    * (overra**(powrn-2.0d0))
     &                    + ((aofn**2)/(powrn-1.0d0)
     &                    )*(overra**(powrn-1.0d0)) )
     &                    - (sigma**6)*((1.0d0/3.0d0)*(overra**3)
     &                    + (aofn/2.0d0)*(overra**4)
     &                    +((aofn**2)/5.0d0)*(overra**5) ))
                     call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
c                    --- from Peter Gordon 
                     dvalue = (twopi * fecofn/3.0d0)*(
     &                    powrn*(sigma**powrn)*((1.0d0/(powrn-3.0d0))
     &                    *(overra**(powrn-3))
     &                    +(3.0d0*aofn/(powrn-2.0d0)
     &                    )*(overra**(powrn-2.0d0))
     &                    + (3.0d0*(aofn**2)/(powrn-1.0d0))
     &                    *(overra**(powrn-1.0d0))
     &                    + ((aofn**3)/powrn)*(overra**powrn) )
     &                    - 6.0d0*(sigma**6) *(
     &                    (1.0d0/3.0d0)*overra**3
     &                    + (3.0d0*aofn/4.0d0)*overra**4
     &                    + (3.0d0*(aofn**2)/5.0d0)*overra**5
     &                    + ((aofn**3)/6.0d0)*overra**6 ) )
                     call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                  elseif ( localcp .eq. 'Double Exponential' ) then
c                    --- repulsive part
                     aaa = twh_get_nbcoeff(1,ntij)
                     bbb = twh_get_nbcoeff(2,ntij)
                     if ( bbb .ne. 0.0d0 ) then
                        dvalue = 2.0d0*onepi*(aaa/bbb)
     &                       *twh_expon(-bbb*rcut)
     &                       * (rcut*rcut + 2.0d0*rcut/bbb + 2.0d0/bbb)
                        call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                        dvalue = (2.0d0/3.0d0)*onepi*aaa
     &                       *twh_expon(-bbb*rcut)
     &                       *(rcut**3 + ((3.0d0*rcut**2)/bbb)
     &                       + ((6.0d0*rcut)/(bbb**2))
     &                       + ((6.0d0)/(bbb**3)))
                        call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                     endif
c                    --- attractive part
                     aaa = -twh_get_nbcoeff(1,ntij)
                     bbb = twh_get_nbcoeff(2,ntij)
                     if ( bbb .ne. 0.0d0 ) then
                        dvalue = 2.0d0*onepi*(aaa/bbb)
     &                       *twh_expon(-bbb*rcut)
     &                       * (rcut*rcut + 2.0d0*rcut/bbb + 2.0d0/bbb)
                        call twh_coru(GLB_INCR,imolty,jmolty,dvalue)
                        dvalue = (2.0d0/3.0d0)*onepi*aaa
     &                       *twh_expon(-bbb*rcut)
     &                       *(rcut**3 + ((3.0d0*rcut**2)/bbb)
     &                       + ((6.0d0*rcut)/(bbb**2))
     &                       + ((6.0d0)/(bbb**3)))
                        call twh_corp(GLB_INCR,imolty,jmolty,dvalue)
                     endif
                  else
                     stop 'undefined classical potential in TAIL'
                  endif
               enddo
            enddo
         enddo
      enddo

      return
      end
