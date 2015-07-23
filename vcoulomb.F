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
      function twh_vcoulomb( ibox,rij,imolty,iunit,jmolty,junit,
     &   ltorsion, vforeign, dvcoulomb,lovrlap)
c     ******************************************************************
c     * Computes the nonbond coulombic energy for two charges in a     *
c     * simulation box (scaled if scp_scale_coulomb is true).          *
c     * dvcoulomb returns dU/d[lambda_c], if scp_eval_derivative true, *
c     * 0 otherwise.                                                   *
c     * vforeign returns the potential at 'foreign' lambda values,     *
c     * corresponding to a different hamiltonian at this configuration *
c     * if we are computing foreign lambdas.  It is set to 0.0         *
c     * otherwise.  (see scaled_potential.F)                           *
c     *                                                                *
c     * Value returned is in internal [EE] units, with distances in    *
c     * angstroms and charges in multiples of electron charge e.       *
c     * Multiply by eetokcal to obtain energies in [kcal/mol] and by   *
c     * qqfact for energies in [K/molecule].                           *
c     *                                                                *
c     * originally written 01-08-2002 by M.G. Martin                   *
c     * last update 08-09-2011 by M.G. Martin                          *
c     ******************************************************************
      implicit none
#define FUNCTION_CMP_GT
#define FUNCTION_CMP_LT
#define FUNCTION_DERFC
#define FUNCTION_VCOULOMB
#define FUNCTION_GET_FOREIGN_LAMBDA
#define FUNCTION_GET_NATIVE_LAMBDA
#define FUNCTION_SCALE_ATOMS
#include "functions.h"

c     --- variables passed to/from the function
      logical ltorsion,lovrlap
      integer ibox, imolty, jmolty, iunit, junit
      double precision rij, dvcoulomb
      double precision vforeign
      dimension vforeign(MAX_FOREIGN_LAMBDA)
c     --- local variables
      logical scp_scale_coulomb,lewald
      double precision qqa,qqb,calp,rcelect,rmin
      integer ntii, ntjj, l,isolvtype,scp_num_foreign_lambda

c     --- retrieve constants
      call twh_scp_num_foreign_lambda(GLB_GET,scp_num_foreign_lambda)
      call twh_qqatom(GLB_GET,imolty,iunit,qqa)
      call twh_qqatom(GLB_GET,jmolty,junit,qqb)
      call twh_lewald(GLB_GET,lewald)

      if ( ibox .eq. 0 ) then
c        --- we compute all interactions using a straight Coulomb
c        --- potential in the 0 box.  Set rcelect to be larger than rij
         rcelect = rij + 1.0d0
      else
c        --- get rcelect for this simulation box
         call twh_rcelect(GLB_GET,ibox,rcelect)
      endif

      dvcoulomb = 0.0d0

      if ( twh_cmp_gt(rcelect,rij) ) then
c        --- check hard overlap
         call twh_rmin(GLB_GET,rmin)
         if ( twh_cmp_lt(rij, rmin) ) lovrlap = .true.
c        --- compute coulombic interactions
         if ( ibox .eq. 0 .or. ( .not. lewald) .or. ltorsion ) then
c           --- straight coulombic interactions
            twh_vcoulomb = qqa*qqb/rij
         else
c           --- use ewald scaling
            call twh_ewald_calp(GLB_GET,ibox,calp)
            twh_vcoulomb = twh_derfc(calp*rij)*qqa*qqb/rij
         endif
         
         call twh_isolvtype(GLB_GET,isolvtype)
         if ( isolvtype .eq. SOLV_SASA ) then
c           --- SASA implicit solvent model applies a distance 
c           --- dependent dielectric constant (specifically Ferrara,
c           --- Apostolakis, Caflisch use e(r) = 2r
            twh_vcoulomb = 0.5d0*twh_vcoulomb/rij
         elseif ( isolvtype .eq. SOLV_EEF1 ) then
c           --- EEF1 implicit solvent model applies a distance 
c           --- dependent dielectric constant (specifically Lazaridis
c           --- and Karplus use e(r) = r)
            twh_vcoulomb = twh_vcoulomb/rij
         endif


c        If evaluating foreign energies, push vcoulomb into them now.
c        These will be scaled later if necessary.  Note that foreign
c        energies will contain all interaction energies, even those 
c        which are not scaled.
         do l = 1, scp_num_foreign_lambda
            vforeign(l) = twh_vcoulomb
         enddo
c        --- do scaling by lambda here.
c        --- Scaling is linear in lambda; any solvation effects are 
c        --- scaled linearly as well.
         call twh_scp_scale_coulomb(GLB_GET,scp_scale_coulomb)
         if (scp_scale_coulomb) then
            call twh_ntype(GLB_GET,imolty,iunit,ntii)
            call twh_ntype(GLB_GET,jmolty,junit,ntjj)
            if (scp_scale_atoms(ntii, ntjj) ) then 
c           --- lambda derivative is easy for linear scaling
               dvcoulomb = twh_vcoulomb
c              --- evaluate foreign energies first, if necessary
               do l = 1, scp_num_foreign_lambda
                  vforeign(l) = vforeign(l) * 
     &                 scp_get_foreign_lambda(LAMBDA_C, l)
               enddo
c              --- scale native energy
               twh_vcoulomb = scp_get_native_lambda(LAMBDA_C)
     &               *twh_vcoulomb
            endif
         endif
      else
         twh_vcoulomb = 0.0d0
         do l = 1, scp_num_foreign_lambda
            vforeign(l) = 0.0d0
         enddo
      endif

      return
      end
