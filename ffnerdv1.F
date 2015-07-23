#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2014 Marcus G. Martin                       *
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
c     *                                                                *
c     * See the file towhee.F for more information about the code      *
c     ******************************************************************
      subroutine twh_ffnerdv1(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the original     *
c     * version of the NERD forcefield                                 *
c     *                                                                *
c     * originally written 06-29-2004 by N.B. du Preez                 *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal,ntslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'NERDv1'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      type = 0

c     ******************************************************************
c     * NERD Version 1                                                 *
c     ******************************************************************
c     The main difference from NERD Version 1 and Version 2 is a change
c     in the torsional parameters for alkanes
c     Here we list just the NERD Version 1 references
c     -------------------------------------------------
c     jcp1998: S.K. Nath; F.A. Escobedo; J.J. de Pablo; "On the 
c     simulation of vapor-liquid equilibria for alkanes", J. Chem.
c     Phys. 108, 23, 9905-9911 (1998)
c     -------------------------------------------------
c     mp2000: S.K. Nath; J.J. de Pablo; "Simulation of vapour-liquid
c     equilibria for branched alkanes", Mol. Phys. 98, 4, 231-238 (2000)
c     -------------------------------------------------

c     --- we are not using the 1-4 nonbonded terms so set them to zero
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0

c     Ethane methyl group (CH3) united-atom [jcp1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.825d0
      nbclocal(2) = 100.6d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3eth'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     Propane methyl group (CH3) united-atom [jcp1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.857d0
      nbclocal(2) = 102.6d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3pro'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     The "generic" methyl group (CH3) united-atom that is not in 
c     one of the many other categories of "specific" types [jcp1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.91d0
      nbclocal(2) = 104.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3gen'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in isobutane (CH3) united-atom [mp2000]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.88d0
      nbclocal(2) = 78.23d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3isob'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in 2-methylbutane (CH3) united-atom [mp2000]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.90d0
      nbclocal(2) = 79.5d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp32mb'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in "general" methyl side-chains (CH3) united-atom 
c     [mp2000]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.85d0
      nbclocal(2) = 70.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3mesc'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in "general" ethyl side-chains (CH3) united-atom 
c     [mp2000]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.82d0
      nbclocal(2) = 83.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3etsc'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methylene group (CH2) united atom [jpc1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.93d0
      nbclocal(2) = 45.8d0
      masslocal = 14.0268d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH2sp3'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methine group (CH) united-atom [mp2000]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.85d0
      nbclocal(2) = 39.7d0
      masslocal = 13.0189d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHsp3'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ethene (CH2) united-atom [mp2000]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.79d0
      nbclocal(2) = 84.7d0
      masslocal = 14.0268d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH2sp2eth'
      cbnlocal = 'Csp2'
      canlocal = 'CH2sp2'
      ctnlocal = 'CH2sp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD bond lengths - simple harmonic                            *
c     ******************************************************************

c     C-C single bond between carbons
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (0.5d0)*(96500d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- both sp3 hybrids carbon single bond [jcp1998]
      vnlocal(1,1) = 'Csp3'
      vnlocal(1,2) = 'Csp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     C=C double bond between carbons
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (0.5d0)*(96500d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- both sp3 carbon double bond [mp2000]
      vnlocal(1,1) = 'Csp2'
      vnlocal(1,2) = 'Csp2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD angles - simple harmonic                                  *
c     ******************************************************************

c     C* - CH2sp3 - C* [jcp1998]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = (0.5d0)*(62500d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 6
c     CH3 - CH2 - CH3 (propane)
      bnlocal(1,1) = 'CH3sp3'
      bnlocal(1,2) = 'CH2sp3'
      bnlocal(1,3) = 'CH3sp3'
c     CH3 - CH2 - CH2 (linear alkane chains)
      bnlocal(2,1) = 'CH3sp3'
      bnlocal(2,2) = 'CH2sp3'
      bnlocal(2,3) = 'CH2sp3'
c     CH3 - CH2 - CH (linear bits of branched alkane side chains)
      bnlocal(3,1) = 'CH3sp3'
      bnlocal(3,2) = 'CH2sp3'
      bnlocal(3,3) = 'CHsp3'
c     CH2 - CH2 - CH2 (linear alkane chains/sides interiors)
      bnlocal(4,1) = 'CH2sp3'
      bnlocal(4,2) = 'CH2sp3'
      bnlocal(4,3) = 'CH2sp3'
c     CH2 - CH2 - CH 
      bnlocal(5,1) = 'CH2sp3'
      bnlocal(5,2) = 'CH2sp3'
      bnlocal(5,3) = 'CHsp3'
c     CH - CH2 - CH 
      bnlocal(6,1) = 'CHsp3'
      bnlocal(6,2) = 'CH2sp3'
      bnlocal(6,3) = 'CHsp3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     C* - CHsp3 - C* centered at a CH unit [mp2000]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.4d0
      bclocal(1) = (0.5d0)*(62500d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 6
c     CH3 - CH - CH3
      bnlocal(1,1) = 'CH3sp3'
      bnlocal(1,2) = 'CHsp3'
      bnlocal(1,3) = 'CH3sp3'
c     CH3 - CH - CH2
      bnlocal(2,1) = 'CH3sp3'
      bnlocal(2,2) = 'CHsp3'
      bnlocal(2,3) = 'CH2sp3'
c     CH3 - CH - CH
      bnlocal(3,1) = 'CH3sp3'
      bnlocal(3,2) = 'CHsp3'
      bnlocal(3,3) = 'CHsp3'
c     CH2 - CH - CH2
      bnlocal(4,1) = 'CH2sp3'
      bnlocal(4,2) = 'CHsp3'
      bnlocal(4,3) = 'CH2sp3'
c     CH2 - CH - CH
      bnlocal(5,1) = 'CH2sp3'
      bnlocal(5,2) = 'CHsp3'
      bnlocal(5,3) = 'CHsp3'
c     CH - CH - CH
      bnlocal(6,1) = 'CHsp3'
      bnlocal(6,2) = 'CHsp3'
      bnlocal(6,3) = 'CHsp3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD torsions using standard OPLS cosine series                *
c     * k0 + k1*(1 + Cos(phi)) + k2*(1 - Cos(2*phi))                   *
c     * + k3*(1 + Cos(3*phi))                                          *
c     ******************************************************************

c     C* - CH2sp3  - CH2sp3 - C* linear alkanes [jcp1998]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 0.0d0
      tclocal(1) = 355.04d0
      tclocal(2) = -68.19d0
      tclocal(3) = 701.32d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'CH2sp3'
      tnlocal(1,4) = 'Cwild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C* - CH2sp3  - CHsp3 - C* branched alkanes [mp2000]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 1416.3d0
      tclocal(1) = 398.3d0
      tclocal(2) = 139.12d0
      tclocal(3) = -901.2d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Cwild - CH2sp3 CHsp3 - Cwild
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'CHsp3'
      tnlocal(1,4) = 'Cwild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
