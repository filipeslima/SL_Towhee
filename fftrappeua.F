#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2014 Marcus G. Martin                       *
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
      subroutine twh_fftrappeua(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nbitypes,noftypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the trappe-ua    *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
     &        ,nbitypes,noftypes
c     --- local variables
      integer ii,type,isame
      character*15 ordertor
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
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)
      character*(FFNAMELEN) ofnlocal
      dimension ofnlocal(MAXOFNAME)
      double precision onepi

      ffname = 'TraPPE-UA'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     **************************************************************
c     * TraPPE-UA                                                  *
c     **************************************************************
c     References:
c     [TraPPE-1] 
c     M.G. Martin, J.I. Siepmann, 
c     `Transferable potentials for phase equilibria. 1. 
c     United-atom description of n-alkanes'
c     J. Phys. Chem. B 102, 2569-2577 (1998)
c     [tca1998]  
c     M.G. Martin, J.I. Siepmann,
c     `Calculating Gibbs free energies of transfer from
c     Gibbs ensemble Monte Carlo simulations'
c     Theor. Chem. Acc. 99, 347-350 (1998)
c     [TraPPE-2] 
c     M.G. Martin, J.I. Siepmann, 
c     `Novel configurational-bias Monte Carlo method for
c     branched molecules. Transferable potentials for 
c     phase equilibria. 2. united-atom description of
c     branched alkanes'
c     J. Phys. Chem. B 103, 4508-4517 (1999)
c     [TraPPE-4] 
c     C.D. Wick, M.G. Martin, J.I. Siepmann,
c     `Transferable potentials for phase equilibria. 4.
c     United-atom description of linear and branched 
c     alkenes and alkylbenzenes'
c     J. Phys. Chem. B 104, 8008-8016 (2000)
c     [TraPPE-5] 
c     B. Chen, J. Potoff, J.I. Siepmann
c     `Monte Carlo calculations for alcohols and their 
c     mixtures with alkanes.  Transferable potentials 
c     for phase equilibria. 5. United-atom description
c     of primary, secondary, and tertiary alcohols'
c     J. Phys. Chem. B 105, 3090-3104 (2001)
c     [TraPPE-6]
c     J.M. Stubbs; J.J. Potoff; J.I. Siepmann;
c     "Transferable Potentials for Phase Equilbiria. 6. United-Atom
c     Description for Ethers, Glycols, Ketones, and Aldehydes",
c     J. Phys. Chem. B 108, 17596-17605 (2004)
c     [jpcb2004]
c     G. Kamath; F. Cao; J.J. Potoff;
c     "An Improved Force Field for the Prediction of the Vapor-Liquid
c     Equilibria for Carboxylic Acids", J. Phys. Chem. B 108 
c     14130-14136 (2004)
c     [TraPPE-7]
c     C.D. Wick; J.M. Stubbs; N. Rai; J.I. Siepmann;
c     "Transferable Potentials for Phase Equilibria.  7. Primary,
c     Secondary, and Tertiary Amines, Nitroalkanes and Nitrobenzene,
c     Nitriles, Amides, Pyridine, and Pyrimidine"
c     J. Phys. Chem. B 109, 18974-18982 (2005).
c     [TraPPE-8]
c     N. Lubna; G. Kamath; J.J. Potoff; N. Rai; J.I. Siepmann;
c     "Transferable Potentials for Phase Equilibria. 8. United-Atom
c     Description for Thiols, Sulfides, Disulfides, and Thiophene"
c     J. Phys. Chem. B 109, 24100-24107 (2005).

      type = 0
c     *** Argon ***
      
c     (Ar) Argon [Ttca1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.39d0
      nbclocal(2) = 116.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 39.948d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  'Ar'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ar'
      cbnlocal = 'none'
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Carbon (plus bonded hydrogens) ***

c     *** CH4 molecules ***
      
c     (CH4) Methane CH4 (sp3)  [TraPPE-1]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.73d0
      nbclocal(2) = 148.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 4.0d0*1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH4'
      cbnlocal = 'none'
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** CH3 united-atoms ***
      
c     CH3*(sp3) Alkane methyl group (sp3) bonded to anything 
c     -------- [TraPPE-1]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = 98.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 3.0d0*1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3*(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     *** CH2 united-atoms ***

c     CH2**(sp3) Alkane [CH2] (sp3) bonded to anything [TraPPE-1]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.95d0
      nbclocal(2) = 46.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 2.0d0*1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2**(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (CH2=*2) [CH2] alkene C (sp2) bonded to anything [TraPPE-4]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.675d0
      nbclocal(2) = 85.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0269d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2=*(sp2)'
      cbnlocal = 'Csp2'
      canlocal = 'CH2sp2'
      ctnlocal = 'CH2sp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** CH united-atoms ***
      
c     CHccc(sp3) Alkane [CH] (sp3) bonded only to Carbons [TraPPE-2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.68d0
      nbclocal(2) = 10.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0 + 1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHccc(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     CHcco(sp3) Alcohol [CH] (sp3) bonded to 2 Carbons and 1 Oxygen 
c     -------- [TraPPE-5]
c     also bonded to ether [TraPPE-6]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.33d0
      nbclocal(2) = 10.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0191d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHcco(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     CHccs(sp3) thiol or sulfide [CH] (sp3) bonded to 2 Carbons and
c     1 Sulphur
c     -------- [TraPPE-8]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.68d0
      nbclocal(2) = 10.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0191d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHccs(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     CHc=c(sp2) Alkene [CH] (sp2) bonded to carbons [TraPPE-4]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.73d0
      nbclocal(2) = 47.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0191d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHc=c(sp2)'
      cbnlocal = 'Csp2'
      canlocal = 'CHsp2'
      ctnlocal = 'CHsp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     CHc=o(sp2) CH (sp2) in an aldehyde single bonded to C and double
c     bonded to O [TraPPE-6]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.52d0
      nbclocal(2) = 54.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0191d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHc=o(sp2)'
      cbnlocal = 'CHc=o(sp2)'
      canlocal = 'CHc=o(sp2)'
      ctnlocal = 'CHsp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     CHaa(aro) Aromatic [CH] group bonded to two aromatic carbons 
c     -------- [TraPPE-4]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.695d0
      nbclocal(2) = 50.5d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0191d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CHaa(aro)'
      cbnlocal = 'Caro'
      canlocal = 'CHaro'
      ctnlocal = 'CHaro'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     *** C atoms ***
      
c     Ccccc(sp3) Alkane [C] (sp3) bonded only to Carbons [TraPPE-2]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 6.4d0
      nbclocal(2) = 0.5d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ccccc(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = 'Csp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Cccco(sp3) Alcohol [C] (sp3) bonded to 3 Carbons and 1 Oxygen 
c     -------- [TraPPE-5]
c     also ether [TraPPE-6]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 5.8d0
      nbclocal(2) = 0.5d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cccco(sp3)'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = 'Csp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     Ccc=c(sp2) Alkene [C] (sp2) single bonded to 2 C and double bonded
c     to 1 C [TraPPE-4]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.85d0
      nbclocal(2) = 20.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ccc=c(sp2)'
      cbnlocal = 'Csp2'
      canlocal = 'Csp2'
      ctnlocal = 'Csp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Ccc=o(sp2) Alkene [C] (sp2) single bonded to 2 C and double
c     bonded to 1 O [TraPPE-6]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.82d0
      nbclocal(2) = 40.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ccc=o(sp2)'
      cbnlocal = 'Ccc=o(sp2)'
      canlocal = 'Ccc=o(sp2)'
      ctnlocal = 'Csp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Caac(aro) Aromatic [C] (ar) bonded to two aromatic carbons 
c     -------- [TraPPE-4]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.88d0
      nbclocal(2) = 21.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Caac(aro)'
      cbnlocal = 'Caro'
      canlocal = 'Caro'
      ctnlocal = 'Caro'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Caaa(aro) Aromatic [C] (ar) bonded to 3 aromatic carbons 
c     -------- [TraPPE-4]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.70d0
      nbclocal(2) = 30.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Caaa(aro)'
      cbnlocal = 'Caro'
      canlocal = 'Caro'
      ctnlocal = 'Caro'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Cco=o(sp2) C in a carboxylic acid single bonded to 1 C, 1 O
c     and double bonded to 1 O
c     -------- [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.90d0
      nbclocal(2) = 41.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cco=o(sp2)'
      cbnlocal = 'Cco=o(sp2)'
      canlocal = 'Cco=o(sp2)'
      ctnlocal = 'Cco=o(sp2)'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     *** helium ***

c     (He) Helium [tca1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.11d0
      nbclocal(2) = 4.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 4.0026d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  'He'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'He'
      cbnlocal = 'none'
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Hydrogen ***

c     (Ho) Alcohol [H] bonded to 1 Oxygen [TraPPE-5]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ho'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (Hoacid) Carboxylic acid [H] bonded to 1 Oxygen [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Hoacid'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (Hs) thiol hydrogen bonded to sulpher [TraPPE-8]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 1.0079d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Hs'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** Oxygen ***

c     Och(sp3) Alcohol [O] bonded to 1 C and 1 H [TraPPE-5]
c     Och(sp3) carboxylic acid O bonded to 1 C and 1 H [jpbc2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.020d0
      nbclocal(2) = 93.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Och(sp3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Occ(sp3) ether [O] single bonded to 2 C [TraPPE-6]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.80d0
      nbclocal(2) = 55.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Occ(sp3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     O=c(sp2) aldehyde O (sp2) double bonded to C [TraPPE-6]
c     O=c(sp2) ketone O (sp2) double bonded to C [TraPPE-6]
c     O=c(sp2) carboxylic acid O (sp2) double bonded to C [jpcb2004]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.05d0
      nbclocal(2) = 79.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal = ' O'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O=c(sp2)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     Sulphur in a thiol [TraPPE-8]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.62d0
      nbclocal(2) = 232.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 32.065d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Sch(sp3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     Sulphur in sulfide [TraPPE-8]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.58d0
      nbclocal(2) = 199.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 32.065d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Scc(sp3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     Sulphur in disulfide
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.72d0
      nbclocal(2) = 165.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 32.065d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Scs(sp3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     Sulphur in thiophene
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.60d0
      nbclocal(2) = 180.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 32.065d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' S'
      bplocal = 'arom'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Scc(aro)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE bond lengths                                            *
c     ******************************************************************

c     TraPPE fixed bond length [TraPPE-1]
c     [CHx](sp3) -- [CHy](sp3,sp2,aro) single bonds
c     --- from table 1.6 of McMurry Organic Chemistry 3rd edition
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.54d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     Csp3 - Csp3 [TraPPE-1]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Csp3'
c     Csp3 - Csp2 [TraPPE-4]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Csp2'
c     Csp3 - Caro [TraPPE-4]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Caro'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- TraPPE fixed bond length [TraPPE-4]
c     --- [CHx](sp2) == [CHy](sp2) double bonds 
c     --- from table 1.6 of McMurry Organic Chemistry 3rd edition
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.33d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     Csp2 - Csp2 [TraPPE-4]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp2'
      vnlocal(isame,2) = 'Csp2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- TraPPE fixed bond length [TraPPE-4]
c     --- [CHx](aro) @@ [CHy](aro) aromatic 1.5 bonds
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.40d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     Caro - Caro [TraPPE-4]
      isame = isame + 1
      vnlocal(isame,1) = 'Caro'
      vnlocal(isame,2) = 'Caro'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- fixed bond length for CHx(sp3)-O in alcohol [TraPPE-5]
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.43d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     Csp3 - Och(sp3) [TraPPE-5]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Och(sp3)'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- fixed bond length for O-H in alcohol [TraPPE-5]
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.945d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     O-H [TraPPE-5]
      vnlocal(1,1) = 'Och(sp3)'
      vnlocal(1,2) = 'Ho'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.41d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     CHx - O ether [TraPPE-6]
      vnlocal(1,1) = 'Occ(sp3)'
      vnlocal(1,2) = 'Csp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.229d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     C=O ketone [TraPPE-6]
      vnlocal(1,1) = 'Ccc=o(sp2)'
      vnlocal(1,2) = 'O=c(sp2)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.217d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     CH=O aldehyde [TraPPE-6]
      vnlocal(1,1) = 'CHc=o(sp2)'
      vnlocal(1,2) = 'O=c(sp2)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.52d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CHx - C aldehyde[TraPPE-6]
      isame = isame + 1
      vnlocal(isame,1) = 'CHc=o(sp2)'
      vnlocal(isame,2) = 'Csp3'
c     CHx - C ketone[TraPPE-6]
      isame = isame + 1
      vnlocal(isame,1) = 'Ccc=o(sp2)'
      vnlocal(isame,2) = 'Csp3'
c     CHx - C carboxylic acid [jpcb2004]
      isame = isame + 1
      vnlocal(isame,1) = 'Cco=o(sp2)'
      vnlocal(isame,2) = 'Csp3'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.214d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     C=O carboxylic acid [jpcb2004]
      vnlocal(1,1) = 'Cco=o(sp2)'
      vnlocal(1,2) = 'O=c(sp2)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.364d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     C-O carboxylic acid [jpcb2004]
      vnlocal(1,1) = 'Cco=o(sp2)'
      vnlocal(1,2) = 'Och(sp3)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.970d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     H-O carboxylic acid [jpcb2004]
      vnlocal(1,1) = 'Hoacid'
      vnlocal(1,2) = 'Och(sp3)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.34d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     S-H thiol [TraPPE-8]
      isame = isame + 1
      vnlocal(isame,1) = 'Sch(sp3)'
      vnlocal(isame,2) = 'Hs'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.82d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CHx - S(aliphatic) thiol [TraPPE-8]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Sch(sp3)'
c     CHx - S(aliphatic) monosulfide [TraPPE-8]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Scc(sp3)'
c     CHx - S(aliphatic) disulfide [TraPPE-8]
      isame = isame + 1
      vnlocal(isame,1) = 'Csp3'
      vnlocal(isame,2) = 'Scs(sp3)'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.71d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     CH-S (aromatic) [TraPPE-8]
      isame = isame + 1
      vnlocal(isame,1) = 'Caro'
      vnlocal(isame,2) = 'Scc(aro)'
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE-UA angles - simple harmonic                             *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2(sp3) - CHy [TraPPE-1]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CH2sp3'
      bnlocal(isame,3) = 'Cwild'
c     CHy - CHx - S thiol [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Cwild'
      bnlocal(isame,3) = 'Sch(sp3)'
c     CHy - CHx - S monosulfide [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Cwild'
      bnlocal(isame,3) = 'Scc(sp3)'
c     CHy - CHx - S disulfide [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Cwild'
      bnlocal(isame,3) = 'Scs(sp3)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.0d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH(sp3) - CHy [TraPPE-2]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHsp3'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - C(sp3) - CHy [TraPPE-2]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Csp3'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.70d0
      bclocal(1) = (70420.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx(sp2) = CH(sp2) - CHy [TraPPE-4]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHsp2'
      bnlocal(isame,3) = 'Cwild'
c     CHx(sp2) = C(sp2) - CHy(sp3) [TraPPE-4]
c     CHx - C(sp2) - CHy [TraPPE-4]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Csp2'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- aromatic centered bond angles are rigid in TraPPE-UA
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 120.00d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx(aro) - CHy(aro) - CHz [TraPPE-4]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Caro'
      bnlocal(isame,3) = 'Cwild'
c     CHx(aro) - CHy(aro) - CHz [TraPPE-4]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHaro'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (50400.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CHy(sp3) - O alcohol [TraPPE-5]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CH2sp3'
      bnlocal(isame,3) = 'Och(sp3)'
c     CHx - CHy(sp3) - O alcohol [TraPPE-5]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHsp3'
      bnlocal(isame,3) = 'Och(sp3)'
c     CHx - CHy(sp3) - O alcohol [TraPPE-5]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Csp3'
      bnlocal(isame,3) = 'Och(sp3)'
c     CHx - CHy(sp3) - O ether [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CH2sp3'
      bnlocal(isame,3) = 'Occ(sp3)'
c     CHx - CHy(sp3) - O ether [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHsp3'
      bnlocal(isame,3) = 'Occ(sp3)'
c     CHx - CHy(sp3) - O ether [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Csp3'
      bnlocal(isame,3) = 'Occ(sp3)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.5d0
      bclocal(1) = (55400.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - Och(sp3) - Ho [TraPPE-5]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Och(sp3)'
      bnlocal(isame,3) = 'Ho'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.0d0
      bclocal(1) = (60400.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - O(sp3) - CHy [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Occ(sp3)'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.2d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - C(=O) - CHy [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Ccc=o(sp2)'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.4d0
      bclocal(1) = (62500.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - C = O aldehyde [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'CHc=o(sp2)'
      bnlocal(isame,3) = 'O=c(sp2)'
c     CHx - C = O ketone [TraPPE-6]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Ccc=o(sp2)'
      bnlocal(isame,3) = 'O=c(sp2)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.0d0
      bclocal(1) = (17600.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H - O - C carboxylic acid [jpcb2004]
      isame = isame + 1
      bnlocal(isame,1) = 'Hoacid'
      bnlocal(isame,2) = 'Och(sp3)'
      bnlocal(isame,3) = 'Cco=o(sp2)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.0d0
      bclocal(1) = (40300.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O - C = O carboxylic acid [jpcb2004]
      isame = isame + 1
      bnlocal(isame,1) = 'Och(sp3)'
      bnlocal(isame,2) = 'Cco=o(sp2)'
      bnlocal(isame,3) = 'O=c(sp2)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (35300.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O - C - CH3 carboxylic acid [jpcb2004]
      isame = isame + 1
      bnlocal(isame,1) = 'Och(sp3)'
      bnlocal(isame,2) = 'Cco=o(sp2)'
      bnlocal(isame,3) = 'CH3sp3'
c     O - C - CH2 carboxylic acid (assumed from [jpcb2004])
      isame = isame + 1
      bnlocal(isame,1) = 'Och(sp3)'
      bnlocal(isame,2) = 'Cco=o(sp2)'
      bnlocal(isame,3) = 'CH2sp3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.0d0
      bclocal(1) = (40300.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O = C - CH3 carboxylic acid [jpcb2004]
      isame = isame + 1
      bnlocal(isame,1) = 'O=c(sp2)'
      bnlocal(isame,2) = 'Cco=o(sp2)'
      bnlocal(isame,3) = 'CH3sp3'
c     O = C - CH2 carboxylic acid (assumed from [jpcb2004])
      isame = isame + 1
      bnlocal(isame,1) = 'O=c(sp2)'
      bnlocal(isame,2) = 'Cco=o(sp2)'
      bnlocal(isame,3) = 'CH2sp3'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 96.0d0
      bclocal(1) = (33.830d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHx - S - H thiol [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Sch(sp3)'
      bnlocal(isame,3) = 'Hs'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 99.0d0
      bclocal(1) = (45550.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHy - S - CHx general case monosulfide [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Scc(sp3)'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.0d0
      bclocal(1) = (45550.0d0)*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CHy - S - S general case disulfide [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Scs(sp3)'
      bnlocal(isame,3) = 'Scs(sp3)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 111.5d0
      bclocal(1) = 1.0d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH(aro) - CH(aro) - S (thiophene) [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Caro'
      bnlocal(isame,2) = 'Caro'
      bnlocal(isame,3) = 'Scc(aro)'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 92.2d0
      bclocal(1) = 1.0d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     CH - S - CH (thiophene) [TraPPE-8]
      isame = isame + 1
      bnlocal(isame,1) = 'Caro'
      bnlocal(isame,2) = 'Scc(aro)'
      bnlocal(isame,3) = 'Caro'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE torsions using standard OPLS cosine series              *
c     * k0 + k1*(1 + Cos(phi)) + k2*(1 - Cos(2*phi))                   *
c     * + k3*(1 + Cos(3*phi)                                           *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 355.03d0
      tclocal(2) = -68.19d0
      tclocal(3) = 791.32d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2(sp3) - CH2(sp3) - CHy [TraPPE-1]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Cwild'
c     surprisingly this is also used in carboxylic acids
c     CHx - CH2(sp3) = Cco=o(sp2) - Och(sp3)
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Cco=o(sp2)'
      tnlocal(isame,4) = 'Och(sp3)'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     -- CH split OPLS torsion Siepmann et al Mol Phys 90 687-693 (1997)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = -251.06d0
      tclocal(1) =  428.73d0
      tclocal(2) = -111.85d0
      tclocal(3) =  441.27d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2(sp3) - CH(sp3) - CHy [TraPPE 2]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CHsp3'
      tnlocal(isame,4) = 'Cwild'
c     CHx - CH(sp3) - CH(sp3) - CHy [TraPPE 2]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHsp3'
      tnlocal(isame,3) = 'CHsp3'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- Mundy et al, Faraday Disc 104, 17-36 (1996)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) =    0.00d0
      tclocal(1) =    0.00d0
      tclocal(2) =    0.00d0
      tclocal(3) =  461.29d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- torsional parameters for segment containing a quaternary 
c     --- carbon as a mid-segment 
c     CHx - CH2(sp3) - C(sp3) - CHy [TraPPE 2]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Csp3'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- from jorgensen JACS 106 22 6638-6646 (1984)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) =    688.5d0
      tclocal(1) =     86.36d0
      tclocal(2) =   -109.77d0
      tclocal(3) =   -282.24d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx -- CH2(sp3) -- CH(sp2) == CHy(sp2) [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'CH2sp2'
c     CHx -- CH2(sp3) -- CH(sp2) == CHy(sp2) [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'CHsp2'
c     CHx -- CH2(sp3) -- CH(sp2) == CHy(sp2) [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Csp2'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- alcohols [TraPPE-5]
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.00d0
      tclocal(1) = 176.62d0
      tclocal(2) = -53.34d0
      tclocal(3) = 769.93d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CHx - CH2(sp3) - CH2(sp3) - OH [TraPPE 5]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Och(sp3)'
c     CHx - CH2 - CH2 - O ether [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Occ(sp3)'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- J. Phys. Chem. B 105 3093-3104 (2001)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.00d0
      tclocal(1) = 209.82d0
      tclocal(2) = -29.17d0
      tclocal(3) = 187.93d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CHx - CH2(sp3) - O - H [TraPPE 5]
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'Och(sp3)'
      tnlocal(1,4) = 'Ho'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- J. Phys. Chem. B 105 3093-3104 (2001)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 215.96d0
      tclocal(1) = 197.33d0
      tclocal(2) = 31.46d0
      tclocal(3) = -173.92d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CHx - CH(sp3) - O - H [TraPPE 5]
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CHsp3'
      tnlocal(1,3) = 'Och(sp3)'
      tnlocal(1,4) = 'Ho'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- J. Phys. Chem. B 105 3093-3104 (2001)
      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 0.0d0
      tclocal(2) = 0.0d0
      tclocal(3) = 163.56d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- CHx - C(sp3) - O - H [TraPPE 5]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Csp3'
      tnlocal(isame,3) = 'Och(sp3)'
      tnlocal(isame,4) = 'Ho'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 725.35d0
      tclocal(2) = -163.75d0
      tclocal(3) = 558.20d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CHy - O - CHz [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Csp3'
      tnlocal(isame,3) = 'Occ(sp3)'
      tnlocal(isame,4) = 'Cwild'
c     CHx - CHy - O - CHz [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHsp3'
      tnlocal(isame,3) = 'Occ(sp3)'
      tnlocal(isame,4) = 'Cwild'
c     CHx - CHy - O - CHz [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Occ(sp3)'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 503.24d0
      tclocal(1) = 0.0d0
      tclocal(2) = -251.62d0
      tclocal(3) = 1006.47d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     O - CH2 - CH2 - O ethers [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Occ(sp3)'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Occ(sp3)'
c     O - CH2 - CH2 - O ether/alcohol [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Occ(sp3)'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Och(sp3)'
c     O - CH2 - CH2 - O diol [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Och(sp3)'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Och(sp3)'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 2035.58d0
      tclocal(1) = -736.90d0
      tclocal(2) = 57.84d0
      tclocal(3) = -293.23d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CHy - C = O ketones [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Cwild'
      tnlocal(isame,3) = 'Csp2'
      tnlocal(isame,4) = 'O=c(sp2)'
c     CHx - CHy - C = O aldehydes [TraPPE-6]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Cwild'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'O=c(sp2)'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 52.88d0
      tclocal(1) = -52.05d0
      tclocal(2) = 18.46d0
      tclocal(3) = 306.36d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2 - S - H [TraPPE-8]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Sch(sp3)'
      tnlocal(isame,4) = 'Hs'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 367.6d0
      tclocal(2) = -270.18d0
      tclocal(3) = 581.64d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHy - CHx - S - CHx monosulfide [TraPPE-8]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Cwild'
      tnlocal(isame,3) = 'Scc(sp3)'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 1763.9d0
      tclocal(1) = 398.0d0
      tclocal(2) = -2181.6d0
      tclocal(3) = 505.25d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - S - S - CHx disulfide [TraPPE-8]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Scs(sp3)'
      tnlocal(isame,3) = 'Scs(sp3)'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = -251.06d0
      tclocal(1) = 428.73d0
      tclocal(2) = -111.85d0
      tclocal(3) = 441.27d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2 - CH - S [TraPPE-8]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'CHsp3'
      tnlocal(isame,4) = 'Swild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 0.0d0
      tclocal(2) = 0.0d0
      tclocal(3) = 461.29d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2 - C - S [TraPPE-8]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Csp3'
      tnlocal(isame,4) = 'Swild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 0.0d0
      tclocal(1) = 0.0d0
      tclocal(2) = 0.0d0
      tclocal(3) = 400.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CH2 - C - S [TraPPE-8]
      isame = isame + 1
      tnlocal(isame,1) = 'CH3sp3'
      tnlocal(isame,2) = 'CHsp3'
      tnlocal(isame,3) = 'Sch(sp3)'
      tnlocal(isame,4) = 'Hs'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     ******************************************************************
c     * TraPPE harmonic torsions                                       *
c     * torcoeff0* ( phi - torcoeff1 )^2                               *
c     ******************************************************************
c     --- torsional parameters for TRANS alkene X-(CH sp2)-(CH sp2)-Y
c     --- used in harmonic potential E = torcoeff0 
c     --- * (theta - torcoeff1)**2
c     --- J. Phys. Chem. B 104, 8008-8016 (2000) [TraPPE 4]

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_torstyle(GLB_SET,type,1)
      tclocal(0) =  13400.0d0
      tclocal(1) =  onepi
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ordertor = 'trans'
      call twh_torsionorder(GLB_SET,type,ordertor)
      isame = 0
c     CHx - C(sp2) = C(sp2) - CHy trans
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHsp2'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Cwild'
c     CHx - C(sp2) = C(sp2) - CHy trans
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Csp2'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Cwild'
c     CHx - C(sp2) = C(sp2) - CHy trans
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Csp2'
      tnlocal(isame,3) = 'Csp2'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     --- torsional parameters for CIS alkene X-(CH sp2)-(CH sp2)-Y
c     --- used in harmonic potential E = torcoeff0 
*     --- (theta - torcoeff1)**2
c     --- J. Phys. Chem. B 104, 8008-8016 (2000) [TraPPE 4]
      type = type + 1
      call twh_torstyle(GLB_SET,type,1)
      tclocal(0) =  12400.0d0
      tclocal(1) =  0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ordertor = 'cis'
      call twh_torsionorder(GLB_SET,type,ordertor)
      isame = 0
c     CHx - C(sp2) = C(sp2) - CHy cis
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHsp2'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Cwild'
c     CHx - C(sp2) = C(sp2) - CHy cis
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Csp2'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Cwild'
c     CHx - C(sp2) = C(sp2) - CHy cis
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Csp2'
      tnlocal(isame,3) = 'Csp2'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     ******************************************************************
c     * TraPPE simple cosine term                                      *
c     * vtorsion = torcoeff(x,0)* [1 - cos(2*phi + torcoeff(x,1))]     *
c     ******************************************************************
c     --- C.D. Wick; M.G. Martin; J.I. Siepmann;
c     --- J. Phys. Chem. B 104, 8008-8016 (2000)

      type = type + 1
      call twh_torstyle(GLB_SET,type,7)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 131.0d0
      tclocal(1) = onepi
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx(aro) - C(aro) - CH2(sp3) -- CHx [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'CHaro'
      tnlocal(isame,2) = 'Caro'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Cwild'
c     CHx(aro) - C(aro) - CH2(sp3) -- CHx [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Caro'
      tnlocal(isame,2) = 'Caro'
      tnlocal(isame,3) = 'CH2sp3'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,7)
      tclocal(-1) = 0.5d0
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(0) = 167.0d0
      tclocal(1) = (5.0/3.0)*onepi
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx(aro) - C(aro) - CH(sp2) = CHy [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'CHaro'
      tnlocal(isame,2) = 'Caro'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Cwild'
c     CHx(aro) - C(aro) - CH(sp2) = CHy [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Caro'
      tnlocal(isame,2) = 'Caro'
      tnlocal(isame,3) = 'CHsp2'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     ******************************************************************
c     * multiple rigid torsions                                        *
c     * used for the rigid torsions in the aromatics                   *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.true.)
      call twh_ntorloop(GLB_SET,type,2)
c     --- one-four scaling
      tclocal(-1) = 0.5d0
c     --- tolerance for the rigid dihedral
      tclocal(0) = 1.0d-5
c     --- rigid dihedrals
      tclocal(1) = onepi
      tclocal(2) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CHx - CHx(aro) - CHx(aro) - CHx [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHaro'
      tnlocal(isame,3) = 'CHaro'
      tnlocal(isame,4) = 'Cwild'
c     CHx - CHx(aro) - CHx(aro) - CHx [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'CHaro'
      tnlocal(isame,3) = 'Caro'
      tnlocal(isame,4) = 'Cwild'
c     CHx - CHx(aro) - CHx(aro) - CHx [TraPPE 4]
      isame = isame + 1
      tnlocal(isame,1) = 'Cwild'
      tnlocal(isame,2) = 'Caro'
      tnlocal(isame,3) = 'Caro'
      tnlocal(isame,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

c     ******************************************************************
c     * Potoff torsions are listed in the paper as                     *
c     * c1[1 + Cos(phi + f) + c2[1 - Cos(2 phi)]                       *
c     * However, this is a typo and they actually use                  *
c     * c1[1 + Cos(phi + f) + c2[1 - Cos^2( phi)]                      *
c     * they did not use any 1-4 interactions                          *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,17)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(1) = 630.0d0
      tclocal(2) = 1562.4d0
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     O-C-O-H in carboxlic acids [jcb2004]
      isame = isame + 1
      tnlocal(isame,1) = 'O=c(sp2)'
      tnlocal(isame,2) = 'Cco=o(sp2)'
      tnlocal(isame,3) = 'Och(sp3)'
      tnlocal(isame,4) = 'Hoacid'
c     this is also used for the CHx - CH2 - C = O interaction
c     [personal comm. between J. Potoff and M.G. Martin 2-15-2004]
      isame = isame + 1
      tnlocal(isame,1) = 'CH3sp3'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Cco=o(sp2)'
      tnlocal(isame,4) = 'O=c(sp2)'
c     this is also used for the CHx - CH2 - C = O interaction
c     [personal comm. between J. Potoff and M.G. Martin 2-15-2004]
      isame = isame + 1
      tnlocal(isame,1) = 'CH2sp3'
      tnlocal(isame,2) = 'CH2sp3'
      tnlocal(isame,3) = 'Cco=o(sp2)'
      tnlocal(isame,4) = 'O=c(sp2)'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,17)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(1) = 630.0d0
      tclocal(2) = 1562.4d0
      tclocal(3) = 0.0d0*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CH3-C-O-H in carboxlic acids [jcb2004]
      isame = isame + 1
      tnlocal(isame,1) = 'CH3sp3'
      tnlocal(isame,2) = 'Cco=o(sp2)'
      tnlocal(isame,3) = 'Och(sp3)'
      tnlocal(isame,4) = 'Hoacid'
c     CH2-C-O-H in carboxlic acids [jcb2004]
      isame = isame + 1
      tnlocal(isame,1) = 'CH2sp3'
      tnlocal(isame,2) = 'Cco=o(sp2)'
      tnlocal(isame,3) = 'Och(sp3)'
      tnlocal(isame,4) = 'Hoacid'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * TraPPE-UA bond increments                                      *
c     * inferred from the charge examples in the various papers        *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [TraPPE-1]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Csp3'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [TraPPE-4]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Csp2'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [TraPPE-4]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Caro'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [TraPPE-4]
      qbinlocal(1) = 'Csp2'
      qbinlocal(2) = 'Csp2'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [TraPPE-4]
      qbinlocal(1) = 'Caro'
      qbinlocal(2) = 'Caro'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.435d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     O - H [TraPPE-5]
      qbinlocal(1) = 'Ho'
      qbinlocal(2) = 'Och(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.265d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     CHx - O alcohol [TraPPE-5]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Och(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.25d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     CHx - O ether [TraPPE-6]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Occ(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.043d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     CH(=O) - CHx aldehyde [TraPPE-6]
      qbinlocal(1) = 'CHc=o(sp2)'
      qbinlocal(2) = 'Csp3'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.482d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     CH(=O) = O aldehyde [TraPPE-6]
      qbinlocal(1) = 'CHc=o(sp2)'
      qbinlocal(2) = 'O=c(sp2)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.424d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C(=O) = O ketone [TraPPE-6]
      qbinlocal(1) = 'Ccc=o(sp2)'
      qbinlocal(2) = 'O=c(sp2)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.45d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     O=(C) carboxylic acid [jpcb2004]
      qbinlocal(1) = 'O=c(sp2)'
      qbinlocal(2) = 'Cco=o(sp2)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.12d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     CHx - (C) carboxylic acid [jpcb2004]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Cco=o(sp2)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.09d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     O(sp3) - (C) carboxylic acid [jpcb2004]
      qbinlocal(1) = 'Och(sp3)'
      qbinlocal(2) = 'Cco=o(sp2)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.37d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - O carboxylic acid [jpcb2004]
      qbinlocal(1) = 'Hoacid'
      qbinlocal(2) = 'Och(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.171d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - S thiol [TraPPE-8]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Sch(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.15d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - S mono-sulfide [TraPPE-8]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Scc(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.177d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - S disulfide [TraPPE-8]
      qbinlocal(1) = 'Csp3'
      qbinlocal(2) = 'Scs(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.206d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - S thiol [TraPPE-8]
      qbinlocal(1) = 'Hs'
      qbinlocal(2) = 'Sch(sp3)'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - S thiophene [TraPPE-8]
      qbinlocal(1) = 'Caro'
      qbinlocal(2) = 'Scc(aro)'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

c     ******************************************************************
c     * special 1-5 interactions                                       *
c     ******************************************************************

      type = type + 1
      call twh_ofstyle(GLB_SET,type,2)
      call twh_offf(GLB_SET,type,ffname)
c     --- a/r^12 power repulsion
      call twh_ofcoeff(GLB_SET,type,1,4.0d7)
c     different nonbonded parameters for a 1-5 interaction between
c     an explicit alcohol hydrogen and an ether oxygen
      ofnlocal(1) = 'Ho'
      ofnlocal(2) = 'wild'
      ofnlocal(3) = 'wild'
      ofnlocal(4) = 'wild'
      ofnlocal(5) = 'Occ(sp3)'
      call twh_set_ofnames(type,ofnlocal)

      type = type + 1
      call twh_ofstyle(GLB_SET,type,2)
      call twh_offf(GLB_SET,type,ffname)
c     --- a/r^12 power repulsion
      call twh_ofcoeff(GLB_SET,type,1,7.5d7)
c     different nonbonded parameters for a 1-5 interaction between
c     an explicit alcohol hydrogen and an alcohol oxygen
      ofnlocal(1) = 'Ho'
      ofnlocal(2) = 'wild'
      ofnlocal(3) = 'wild'
      ofnlocal(4) = 'wild'
      ofnlocal(5) = 'Och(sp3)'
      call twh_set_ofnames(type,ofnlocal)

      noftypes = type
      call twh_ffcheck(lfailure,7,type)
      if ( lfailure ) return

      return
      end
