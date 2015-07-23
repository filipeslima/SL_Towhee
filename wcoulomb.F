#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2010 Marcus G. Martin                       *
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
      function twh_wcoulomb( ibox,rij,imolty,iunit,jmolty,junit)
c     ******************************************************************
c     * Computes the coulombic pair virial function divided by r^2     *
c     * as defined by equation 2.60 in Allen and Tildesley             *
c     * w(r) = r dv(r)/dr                                              *
c     * w(r)/r^2 = (dv(r)/dr) /r                                       *
c     *                                                                *
c     * Units returned are [K/molecule]/[angstrom^2]                   *
c     *                                                                *
c     * originally written 05-16-2005 by M.G. Martin                   *
c     * last update 11-16-2010 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_DERFC
#define FUNCTION_EXPON
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_SCALE_ATOMS
#define FUNCTION_WCOULOMB
#include "functions.h"

c     --- variables passed to/from the function
      integer ibox, imolty, jmolty, iunit, junit
      double precision rij
c     --- local variables
      logical scp_scale_coulomb,lewald
      integer ntii, ntjj,isolvtype
      double precision rijsq,qqa,qqb,onepi,qqfact,calp,rcelect

c     --- retrieve constants
      call twh_constant_qqfact(GLB_GET,qqfact)
      call twh_qqatom(GLB_GET,imolty,iunit,qqa)
      call twh_qqatom(GLB_GET,jmolty,junit,qqb)
      if ( ibox .ne. 0 ) call twh_rcelect(GLB_GET,ibox,rcelect)
      call twh_lewald(GLB_GET,lewald)

c     --- initialize
      twh_wcoulomb = 0.0d0
      if ( ibox .eq. 0 .or. ( .not. lewald) ) then
c        --- ideal gas box of infinite extent, or just no ewald
c        --- either way compute explicit coulomb interaction
         twh_wcoulomb = (-1.0d0)*qqa*qqb*qqfact/(rij**3)
      else
         call twh_constant_pi(GLB_GET,onepi)
         if ( twh_cmp_gt(rcelect,rij) ) then
            rijsq = rij*rij
            call twh_ewald_calp(GLB_GET,ibox,calp)
            call twh_isolvtype(GLB_GET,isolvtype)
            if ( isolvtype .eq. SOLV_NONE ) then
               twh_wcoulomb = - (2.0d0*calp*
     &              twh_expon(-calp*calp*rijsq)/dsqrt(onepi)
     &              +twh_derfc(calp*rij)/rij)*qqa*qqb*qqfact/rijsq

            elseif ( isolvtype .eq. SOLV_SASA ) then
c              --- SASA implicit solvent model applies a distance 
c              --- dependent dielectric constant (specifically Ferrara,
c              --- Apostolakis, Caflisch use e(r) = 2r
               twh_wcoulomb = - (2.0d0*calp*
     &              twh_expon(-calp*calp*rijsq)/dsqrt(onepi)
     &              +(2.0d0)*twh_derfc(calp*rij)/rij)
     &              *(0.5d0)*qqa*qqb*qqfact/(rijsq*rij)

            elseif ( isolvtype .eq. SOLV_EEF1 ) then
c              --- EEF1 implicit solvent model applies a distance 
c              --- dependent dielectric constant (specifically Lazaridis
c              --- and Karplus use e(r) = r)
               twh_wcoulomb = - (2.0d0*calp*
     &              twh_expon(-calp*calp*rijsq)/dsqrt(onepi)
     &              +(2.0d0)*twh_derfc(calp*rij)/rij)
     &              *qqa*qqb*qqfact/(rijsq*rij)
            endif
         endif
      endif

c        --- do scaling by lambda here.
c        --- Scaling is linear in lambda; any solvation effects are 
c        --- scaled linearly as well.
      call twh_scp_scale_coulomb(GLB_GET,scp_scale_coulomb)
      if (scp_scale_coulomb) then
         call twh_ntype(GLB_GET,imolty,iunit,ntii)
         call twh_ntype(GLB_GET,jmolty,junit,ntjj)
c        --- scaling is done if one, but not both, atoms are scalable
         if (scp_scale_atoms(ntii,ntjj)) then 
            twh_wcoulomb = scp_get_native_lambda(LAMBDA_C)*twh_wcoulomb
         endif
      endif

      return
      end
