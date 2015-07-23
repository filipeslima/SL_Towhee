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
      subroutine twh_ffnerdv3(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for version 3 of the *
c     * NERD forcefield                                                *
c     *                                                                *
c     * rewritten from ffnerdv2 08-20-2004 by M.G. Martin              *
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
      integer nvslocal,nbslocal
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

      ffname = 'NERDv3'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      type = 0

c     ******************************************************************
c     * NERD Version 3                                                 *
c     ******************************************************************
c     The main difference from NERD Version 2 and Version 3 is a change
c     in the bonded parameters for many of the molecules.  Nonbonded 
c     parameters stayed essentially the same, although a few new types
c     were added
c     Here we list the NERD Version 3 references including those
c     that were already listed for previous versions
c     -NERDv1------------------------------------------------
c     jcp1998: S.K. Nath; F.A. Escobedo; J.J. de Pablo; "On the 
c     simulation of vapor-liquid equilibria for alkanes", J. Chem.
c     Phys. 108, 23, 9905-9911 (1998)
c     -------------------------------------------------
c     mp2000: S.K. Nath; J.J. de Pablo; "Simulation of vapour-liquid
c     equilibria for branched alkanes", Mol. Phys. 98, 4, 231-238 (2000)
c     -NERDv2------------------------------------------------
c     jcp2001a: S.K. Nath; B.J. Banaszak; J.J. de Pablo; "A new 
c     united atom force field for alpha-olefins", J. Chem. Phys. 114, 
c     8, 3612-3616 (2001)
c     -------------------------------------------------
c     jcp2001b: S.K. Nath; R. Khare; "New forcefield parameters for
c     branched hydrocarbons", J. Chem. Phys. 115, 23, 10837-10844 (2001)
c     -------------------------------------------------
c     jpcb2003a: S.K. Nath; "Molecular Simulation of Vapour-Liquid
c     Phase Equilibria of Hydrogen Sulfide and Its Mixtures with 
c     Alkanes", J. Phys. Chem. B 107, 9498-9504 (2003) 
c     -NERDv3------------------------------------------------
c     jpcb2004: R. Khare; A.K. Sum; S.K. Nath; J.J. de Pablo
c     "Simulation of Vapor-Liquid Phase Equilibria of Primary Alcohols
c     and Alcohol-Alkane Mixtures", 108, 10071-10076 (2004)

c     --- we are not using the 1-4 nonbonded terms so set them to zero
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0

c     Methane (CH4) united-atom [jcp2001b]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.73d0
      nbclocal(2) = 148.0d0
      masslocal = 16.0426d0
      elelocal = ' C'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

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
      nbnlocal = 'CH3sp3eth'
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
      nbnlocal = 'CH3sp3pro'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     The "generic" methyl group (CH3) united-atom that is not in 
c     one of the many other categories of "specific" types [jcp1998]
c     --- these units first published in [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.91d0
      nbclocal(2) = (0.8647d0)*CON_KJTOK
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3sp3gen'
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
      nbnlocal = 'CH3sp3isob'
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

c     methyl group in neopentane (CH3) united-atom [jcp2001b]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.91d0
      nbclocal(2) = 74.5d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3sp3neop'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in methyl side-chains of 2,2-dimethylbutane (CH3) 
c     united-atom [jcp2001b]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.87d0
      nbclocal(2) = 74.5d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3sp322db'
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
      nbnlocal = 'CH3sp3mesc'
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
      nbnlocal = 'CH3sp3etsc'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methyl group in propene (CH3) united-atom [jcp2001a]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.85d0
      nbclocal(2) = 100.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3sp3pene'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     methylene group (CH2) united atom [jpc1998]
c     --- current units first published in [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.93d0
      nbclocal(2) = (0.3808d0)*CON_KJTOK
      masslocal = 14.0268d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2sp3'
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
      nbnlocal = 'CH2sp2eth'
      cbnlocal = 'Csp2'
      canlocal = 'CH2sp2'
      ctnlocal = 'CH2sp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     primary alkene (CH2) united-atom [jcp2001a]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.72d0
      nbclocal(2) = 92.5d0
      masslocal = 14.0268d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2sp2prim'
      cbnlocal = 'Csp2'
      canlocal = 'CH2sp2'
      ctnlocal = 'CH2sp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     alkene (CH) united-atom [jcp2001a]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.77d0
      nbclocal(2) = 46.0d0
      masslocal = 13.0189d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHsp2'
      cbnlocal = 'Csp2'
      canlocal = 'CHsp2'
      ctnlocal = 'CHsp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     quaternary alkane carbon (C) atom [jcp2001b]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.91d0
      nbclocal(2) = 17.0d0
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Csp3'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = 'Csp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     hydrogen in SH2 (H) atom [jpcb2003]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.98d0
      nbclocal(2) = 3.9d0
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H_s'
      cbnlocal = 'H'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     hydrogen in alcohols (H) atom [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.98d0
      nbclocal(2) = (0.0324d0)*CON_KJTOK
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H_o'
      cbnlocal = 'H'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Oxygen **

c     oxygen in alcohol [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.98d0
      nbclocal(2) = (0.8979d0)*CON_KJTOK
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Osp3'
      cbnlocal = 'Osp3'
      canlocal = 'Osp3'
      ctnlocal = 'Osp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Sulfur **

c     sulfur in SH2 (S) atom [jpcb2003]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.72d0
      nbclocal(2) = 250.0d0
      masslocal = 32.066d0
      elelocal =  ' S'
      bplocal = 'p'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S_sh2'
      cbnlocal = 'S'
      canlocal = 'S'
      ctnlocal = 'S'
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
c     --- force constant changed in [jpcb2004]
      vclocal(1) = (0.5d0)*(2800d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- both sp3 hybrids carbon single bond [jcp1998]
      ii = 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Csp3'
c     --- sp3 carbon with sp3 carbon single bond [jcp2201a]
      ii = 2
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Csp2'
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

c     S-H in SH2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (0.5d0)*(96500d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- S-H [jpcb2003]
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     CH2-O in alcohols
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.428d0
      vclocal(1) = (0.5d0)*(3300d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Csp3 - Osp3 [jpcb2004]
      vnlocal(1,1) = 'Csp3'
      vnlocal(1,2) = 'Osp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     O-H in alcohols
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.961d0
      vclocal(1) = (0.5d0)*(5200d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Osp3 - H
      vnlocal(1,1) = 'Osp3'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * NERD angles - simple harmonic                                  *
c     ******************************************************************

c     C* - CH2sp3 - C* new parameters in [jpcb2004]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = (0.5d0)*(760d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CH3 - CH2 - CH3 (propane)
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CH3sp3'
c     CH3 - CH2 - CH2 (linear alkane chains)
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CH2sp3'
c     CH3 - CH2 - CH (linear bits of branched alkane side chains)
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH3 - CH2 - C [jpc2001b]
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Csp3'
c     CH2 - CH2 - CH2 (linear alkane chains/sides interiors)
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CH2sp3'
c     CH2 - CH2 - CH 
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH2 - CH2 - C 
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Csp3'
c     CH - CH2 - CH 
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH - CH2 - C
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Csp3'
c     C - CH2 - C 
      ii = ii + 1
      bnlocal(ii,1) = 'Csp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Csp3'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     C* - CHsp3 - C* centered at a CH [mp2000] or a C [jcp2001b]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.4d0
      bclocal(1) = (0.5d0)*(62500d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CH3 - CH - CH3
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'CH3sp3'
c     CH3 - CH - CH2
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'CH2sp3'
c     CH3 - CH - CH
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH3 - CH - C
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'Csp3'
c     CH2 - CH - CH2
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'CH2sp3'
c     CH2 - CH - CH
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH2 - CH - C
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'Csp3'
c     CH - CH - CH
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH - CH - C
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'Csp3'
c     C - CH - C
      ii = ii + 1
      bnlocal(ii,1) = 'Csp3'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'Csp3'
c     CH3 - C - CH3
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'CH3sp3'
c     CH3 - C - CH2
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'CH2sp3'
c     CH3 - C - CH
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH3 - C - C
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'Csp3'
c     CH2 - C - CH2
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'CH2sp3'
c     CH2 - C - CH
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH2 - C - C
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'Csp3'
c     CH - C - CH
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'CHsp3'
c     CH - C - C
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'Csp3'
c     C - C - C
      ii = ii + 1
      bnlocal(ii,1) = 'Csp3'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'Csp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     C-C=C for linear primary alkenes [jcp2001a]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.0d0
      bclocal(1) = (0.5d0)*(62500d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CH2sp2 - CHsp2 - CH3sp3
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp2'
      bnlocal(ii,2) = 'CHsp2'
      bnlocal(ii,3) = 'CH3sp3'
c     CH2sp2 - CHsp2 - CH2sp3
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp2'
      bnlocal(ii,2) = 'CHsp2'
      bnlocal(ii,3) = 'CH2sp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     H2S [jpcb2003]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 91.5d0
      bclocal(1) = (0.5d0)*(62500d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     H - S - H
      ii = ii + 1
      bnlocal(ii,1) = 'H'
      bnlocal(ii,2) = 'S'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     CH2 - CH2 - O in alcohols [jpcb2004]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.0d0
      bclocal(1) = (0.5d0)*(1000d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CH2sp3 - CH2-sp3 - Osp3
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Osp3'
c     CH3sp3 - CH2-sp3 - Osp3
c     --- almost certainly the angle term they used for ethanol as well
c     --- as the text indicates the intramolecular terms were fit to
c     --- pentanol and used for all other alcohols
      ii = ii + 1
      bnlocal(ii,1) = 'CH3sp3'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Osp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     CH2 - O - H in alcohols [jpcb2004]
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.5d0
      bclocal(1) = (0.5d0)*(460d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CH2-sp3 - Osp3 - H
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
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
      tclocal(3) = 791.32d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     Cwild - CH2sp3 - CH2sp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

c     CH3-CH2-CH2-CH2 [jpcb2004]
c     not sure if this is a special case, or supposed to replace 
c     all of the linear alkane torsions.  Implementing as special case 
c     unless I hear differently [MGM, 08-20-2004]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 0.0d0
      tclocal(1) = (2.952d0)*CON_KJTOK
      tclocal(2) = (-0.567d0)*CON_KJTOK
      tclocal(3) = (6.579d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     CH3sp3 - CH2sp3 - CH2sp3 - CH2sp3
      ii = ii + 1
      tnlocal(ii,1) = 'CH3sp3'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'CH2sp3'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

c     C* - CH*sp3  - CHsp3 - C* branched alkanes [mp2000]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 1416.3d0
      tclocal(1) = 398.3d0
      tclocal(2) = 139.12d0
      tclocal(3) = -901.2d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     Cwild - CH2sp3 CHsp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'CHsp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

c     C* - CH*sp3  - CHsp3 - C* branched alkanes [jcp2001b]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 0.0d0
      tclocal(1) = 0.0d0
      tclocal(2) = 0.0d0
      tclocal(3) = 1635.7d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     Cwild - CHsp3 CHsp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CHsp3'
      tnlocal(ii,3) = 'CHsp3'
      tnlocal(ii,4) = 'Cwild'
c     Cwild - CH2sp3 - Csp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'Csp3'
      tnlocal(ii,4) = 'Cwild'
c     Cwild - CHsp3 - Csp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CHsp3'
      tnlocal(ii,3) = 'Csp3'
      tnlocal(ii,4) = 'Cwild'
c     Cwild - Csp3 - Csp3 - Cwild
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'Csp3'
      tnlocal(ii,3) = 'Csp3'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

c     C-C-C=C linear alkenes [jcp2001a]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 47.97d0
      tclocal(1) = 86.31d0
      tclocal(2) = -109.71d0
      tclocal(3) = 282.08d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     Cwild - CH2sp3 - CHsp2 - CH2sp2
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'CHsp2'
      tnlocal(ii,4) = 'CH2sp2'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

c     CH2-CH2-CH2-O in alcohols [jpcb2004]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (1.929d0)*CON_KJTOK
      tclocal(1) = (2.963d0)*CON_KJTOK
      tclocal(2) = (-3.926d0)*CON_KJTOK
      tclocal(3) = (9.311d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     CH2sp3 - CH2sp3 - CH2sp3 - Osp3
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp3'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'CH2sp3'
      tnlocal(ii,4) = 'Osp3'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

c     CH2-CH2-O-H in alcohols [jpcb2004]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.000d0)*CON_KJTOK
      tclocal(1) = (2.987d0)*CON_KJTOK
      tclocal(2) = (0.491d0)*CON_KJTOK
      tclocal(3) = (1.836d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     CH2sp3 - CH2sp3 - Osp3 - H
      ii = ii + 1
      tnlocal(ii,1) = 'CH2sp3'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'Osp3'
      tnlocal(ii,4) = 'H'
c     --- it appears that this term was meant for ethanol as well
c     CH3sp3 - CH2sp3 - Osp3 - H
      ii = ii + 1
      tnlocal(ii,1) = 'CH3sp3'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'Osp3'
      tnlocal(ii,4) = 'H'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
