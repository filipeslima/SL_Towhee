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
      subroutine twh_recipsum(lfailure,ibox,vrecip)
c     ******************************************************************
c     * calculates the total reciprocal space ewald-sum term for volume*
c     * moves                                                          *
c     *                                                                *
c     * originally written in 1998 by Bin Chen.                        *
c     * Last modified 08-04-2011 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_EXPON
#define FUNCTION_GET_NBOXI
#include "functions.h"
c     --- variables passed to/from subroutine
      logical lfailure
      integer ibox
      double precision vrecip
c     --- local variables
c     --- integer scalars
      integer nkx,nky,nkz,nky_min,nky_max,nkz_min,nkz_max
      integer ichain,iunit,imolty,ncount,nchain,nunit,kmax
c     --- double precision scalars
      double precision ksqr,sumr,sumi,arg,invvol,invalpsq4
      double precision n_one,n_two,n_thr,rxui,ryui,rzui,twopi,qqfact
      double precision charge,kx,ky,kz,prefact,calp
      double precision hinv1,hinv2,hinv3,boxvolume
c     --- retrieve constants
      call twh_constant_twopi(GLB_GET,twopi)
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_nchain(GLB_GET,nchain)

c     --- Set up the reciprocal space vectors
      ncount = 0
      vrecip = 0.0d0
      call twh_ewald_calp(GLB_GET,ibox,calp)
      call twh_ewald_kmax(GLB_GET,ibox,kmax)
      if ( kmax .eq. 0 ) then
c        --- no vectors
         call twh_ewald_numvect(GLB_SET,ibox,ncount)
         return
      endif

c     --- set up some ewald constants
      invalpsq4 = 1.0d0/(4.0d0*calp*calp)
c     --- set up the inverse volume along with the charge factor
      call twh_boxvolume(GLB_GET,ibox,boxvolume)
      invvol = qqfact*2.0d0*twopi/boxvolume

      do nkx = 0, kmax
c        --- set the double precision version of nkx
         n_one = dble(nkx)
         nky_max = int(dsqrt(dble(kmax)**2 - n_one**2))
         if (nkx .eq. 0) then
            nky_min = 0
         else
            nky_min = -nky_max
         endif
         do nky = nky_min,nky_max
c           --- set the double precision version of nky
            n_two = dble(nky)
            nkz_max = int( dsqrt(dble(kmax)**2 -
     &           n_one**2 - n_two**2) )
            if (nkx .eq. 0 .and. nky .eq. 0) then
               nkz_min = 1 
            else
               nkz_min = -nkz_max
            endif
            do nkz = nkz_min,nkz_max
c              --- set up the double versions of nkz
               n_thr = dble(nkz)
c              --- increment total number of vectors
               ncount = ncount + 1
c              --- safety check total number of vectors
               if ( ncount .gt. VECTORMAX ) then
                  write(6,*) 'RECIPSUM: ncount ',ncount
                  write(6,*) 'choose a larger vectormax in preproc.h'
                  lfailure = .true.
                  return
               endif
c              --- set up the k vectors, these are written as variable
c              --- h in the Dean Wheeler primer
c              --- note that we are taking the dot products of n with
c              --- hinverse-transpose
               call twh_hinverse(GLB_GET,ibox,1,1,hinv1)
               call twh_hinverse(GLB_GET,ibox,1,2,hinv2)
               call twh_hinverse(GLB_GET,ibox,1,3,hinv3)
               kx = twopi*(n_one*hinv1 + n_two*hinv2 + n_thr*hinv3)
               call twh_hinverse(GLB_GET,ibox,2,1,hinv1)
               call twh_hinverse(GLB_GET,ibox,2,2,hinv2)
               call twh_hinverse(GLB_GET,ibox,2,3,hinv3)
               ky = twopi*(n_one*hinv1 + n_two*hinv2 + n_thr*hinv3)
               call twh_hinverse(GLB_GET,ibox,3,1,hinv1)
               call twh_hinverse(GLB_GET,ibox,3,2,hinv2)
               call twh_hinverse(GLB_GET,ibox,3,3,hinv3)
               kz = twopi*(n_one*hinv1 + n_two*hinv2 + n_thr*hinv3)
               call twh_ewald_kvector(GLB_SET,CRD_REAL,ncount,ibox
     &              ,kx,ky,kz)
c              --- take the dot product of k with itself
               ksqr = kx*kx + ky*ky + kz*kz
c              --- compute the B(h) prefactor 
               prefact = invvol*twh_expon(-ksqr*invalpsq4)/ksqr
               call twh_ewald_prefact(GLB_SET,CRD_REAL,ncount,ibox
     &              ,prefact)
c              --- sum up q*cos and q*sin 
               sumr = 0.0d0
               sumi = 0.0d0
               do ichain = 1,nchain
                  call twh_moltyp(GLB_GET,ichain,imolty)
                  call twh_nunit(GLB_GET,imolty,nunit)
                  if ( twh_get_nboxi(ichain) .eq. ibox) then
                     do iunit = 1,nunit
                        call twh_coordinates(GLB_GET,CRD_REAL,ichain
     &                       ,iunit,rxui,ryui,rzui)
                        call twh_qqatom(GLB_GET,imolty,iunit,charge)
                        arg = kx*rxui + ky*ryui + kz*rzui
                        sumr = sumr + dcos(arg)*charge
                        sumi = sumi + dsin(arg)*charge
                     enddo
                  endif
               enddo
               call twh_ewald_ssumr(GLB_SET,CRD_REAL,ncount,ibox,sumr)
               call twh_ewald_ssumi(GLB_SET,CRD_REAL,ncount,ibox,sumi)
c              --- Potential energy
               vrecip = vrecip + (sumr*sumr + sumi*sumi) * prefact
            enddo
         enddo
      enddo

      call twh_ewald_numvect(GLB_SET,ibox,ncount)
      return 
      end








