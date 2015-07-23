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
c     ******************************************************************
      subroutine twh_ffoplsua(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the OPLS-ua      *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 05-28-2004 by M.G. Martin                   *
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
      double precision aaa,ccc
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

      ffname = 'OPLS-ua'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Geometric')
      type = 0

c     ******************************************************************
c     * OPLS-UA force field                                            *
c     ******************************************************************
c     Reference list
c     [jacs_1984_CSJ] 'Energy Component Analysis for Dilute Aqueous 
c     Solutions fo Li+, Na+, F-, and Cl- Ions', J. Chandrasekhar;
c     D.C. Spellmeyer; W.L. Jorgensen; J. Am. Chem. Soc. 1984, 106,
c     903-901.
c     [jacs_1984_JMS] 'Optimized Intermolecular Potential Functions for
c     Liquid Hydrocarbons', W.L. Jorgensen; J.D. Madura; C.J Swenson
c     J. Am. Chem. Soc. 1984, 106, 6638-6646.
c     [jacs_1984_CJ] 'Solvent Effects on the Relative Energies of 
c     Carbonium Ions.  Solvation and Internal Rotation for the Allyl
c     Cation in Liquid Hydrogen Fluoride', M.E. Cournoyer; 
c     W.L. Jorgensen; J. Am. Chem. Soc. 1984, 106, 5104-5112.
c     [jacs_1985_a] 'Optimized Intermolecular Potential Functions for
c     Amides and Peptides.  Structure and Properties of Liquid Amides'
c     J. Am. Chem. Soc. 1985, 107, 569-578.
c     [jpc_1986_a] 'Optimized Intermolecular Potential Functions for 
c     Liquid Alcohols', W.L. Jorgensen; J. Phys. Chem. 1986, 90,
c     1276-1284.
c     [jpc_1986_b] 'Intermolecular Potential Functions and Monte Carlo
c     Simulations for Liquid Sulfur Compounds' W.L. Jorgensen; 
c     J. Phys. Chem. 1986, 90, 6379-6392.
c     [mp_1988_JB] 'Monte Carlo simulations of liquid acetonitrile
c     with a three-site model' W.J. Jorgensen; J.M. Briggs; Mol. Phys.
c     63, 547-558 (1988)
c     [jpc_1990_JBC] 'Relative Partition Coefficients for Organic
c     Solutes from Fluid Simulations', W.L. Jorgensen; J.M. Briggs;
c     M.L. Contreras; J. Phys. Chem. 94, 1683-1686 (1990).

c     ** Carbon (united atoms including hydrogens) **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.730d0
      nbclocal(2) = (0.294d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 16.0426d0
      elelocal = ' C'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH4 methane [jacs_1984_JMS]
      nbnlocal = 'CH4'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = (0.207d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH3 (C1) ethane [jacs_1984_JMS]
      nbnlocal = 'CH3-C1'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.175d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal =  ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH3 (C2) n-butane [jacs_1984_JMS]
      nbnlocal = 'CH3-C2'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.910d0
      nbclocal(2) = (0.160d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal =  ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH3 (C3) isobutane [jacs_1984_JMS]
c     CH3 bonded to a C=O. ex. NMA [jacs_1985_a]
c     CH3 bonded to a C=O ex. acetone [jpc_1990_JBC]
      nbnlocal = 'CH3-C3'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.960d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal =  ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH3 (C4) neopentane [jacs_1984_JMS]
      nbnlocal = 'CH3-C4'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.80d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal =  ' C'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH3 bonded to a sulfur that is bonded to an atom other 
c     than hydrogen [jpc_1986_b]
c     CH3 bonded to (N secondary) NMA [jacs_1985_a]
c     CH3 bonded to (N tertiary) DMF [jacs_1985_a]
      nbnlocal = 'CH3-X'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Carbon CH2 **
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0268d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH2 (sp3) n-butane [jacs_1984_JMS]
      nbnlocal = 'CH2sp3'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.80d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0268d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH2 (sp3) bonded to a sulfur that is bonded to an atom other 
c     than hydrogen [jpc_1986_b]
      nbnlocal = 'CH2sp3-SR'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.80d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0268d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH2 (sp3) bonded to a N(secondar) ex. C alpha in Glycine
c     CH2 (sp3) bonded to a N(tertiary) ex proline
      nbnlocal = 'CH2sp3-N'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.850d0
      nbclocal(2) = (0.140d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.0268d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH2 (sp2) 1-butene [jacs_1984_JMS]
      nbnlocal = 'CH2sp2'
      cbnlocal = 'Csp2'
      canlocal = 'CH2sp2'
      ctnlocal = 'CH2sp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Carbon CH **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.850d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0189d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH (sp3) isobutane [jacs_1984_JMS]
      nbnlocal = 'CHsp3'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.80d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0189d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH (sp3) bonded to N(secondary) ex. Calpha in Ala [jacs_1985_a]
c     CH (sp3) bonded to N(tertiary) ex. Calpha in Pro [jacs_1985_a]
      nbnlocal = 'CHsp3-N'
      cbnlocal = 'Csp3'
      canlocal = 'CHsp3'
      ctnlocal = 'CHsp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.115d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0189d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH (sp2) 2-butenes [jacs_1984_JMS]
      nbnlocal = 'CHsp2'
      cbnlocal = 'Csp2'
      canlocal = 'CHsp2'
      ctnlocal = 'CHsp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.115d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0189d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH (sp2) CH=O formamide [jacs_1985_1]
c     same as above, but different bonded parameters
      nbnlocal = 'CHsp2-N'
      cbnlocal = 'Csp2-N'
      canlocal = 'CHsp2-N'
      ctnlocal = 'CHsp2-N'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.110d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 13.0189d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH (arom) benzene [jacs_1984_JMS]
      nbnlocal = 'CHarom'
      cbnlocal = 'Carom'
      canlocal = 'Carom'
      ctnlocal = 'Carom'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Carbon (no H) **
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C (sp3) neopentane [jacs_1984_JMS]
      nbnlocal = 'Csp3'
      cbnlocal = 'Csp3'
      canlocal = 'Csp3'
      ctnlocal = 'Csp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C (sp2) isobutene [jacs_1984_JMS]
c     C (sp2) acetone, acetamide, acetic acid [jpc_1990_JBC]
      nbnlocal = 'Csp2'
      cbnlocal = 'Csp2'
      canlocal = 'Csp2'
      ctnlocal = 'Csp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C (sp2) double bonded to O, NMA, Ala [jacs_1985_a]
c     same nonbond parameters as above, but some different bonded terms
      nbnlocal = 'Csp2-N'
      cbnlocal = 'Csp2-N'
      canlocal = 'Csp2-N'
      ctnlocal = 'Csp2-N'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.650d0
c     this paper used units of kj instead of the normal kcal
      nbclocal(2) = (0.628d0)*(CON_KJTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      elelocal =  ' C'
      bplocal = 'sp'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C with a single bond and a triple bond in acetonitrile
      nbnlocal = 'Csp'
      cbnlocal = 'Csp'
      canlocal = 'Csp'
      ctnlocal = 'Csp'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Chlorine **

      type = type + 1
      ii = twh_getnbtype(type,type)
      aaa = (26000.0d0)*CON_KCALTOK
      ccc = (3500.0d0)*CON_KCALTOK
      nbclocal(1) = (aaa/ccc)**(1.0d0/6.0d0)
      nbclocal(2) = (ccc*ccc)/(4.0d0*aaa)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 35.453d0
      elelocal =  'Cl'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     Cl- ion [jacs_1984_CSJ]
      nbnlocal = 'Cl-'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Fluorine **

      type = type + 1
      ii = twh_getnbtype(type,type)
      aaa = (0.5d0)*CON_KCALTOK
      ccc = (1200.0d0)*CON_KCALTOK
      nbclocal(1) = (aaa/ccc)**(1.0d0/6.0d0)
      nbclocal(2) = (ccc*ccc)/(4.0d0*aaa)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 18.998d0
      elelocal =  ' F'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     F- ion [jacs_1984_CSJ]
      nbnlocal = 'F-'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      aaa = (0.33819d0)*CON_KCALTOK
c     this was listed at -0.42189 in [jacs_1984_CJ], but I assumed that
c     must have been a typo as I have never seen a negative value on 
c     such a parameter before
      ccc = (0.42189d0)*CON_KCALTOK
      nbclocal(1) = (aaa/ccc)**(1.0d0/6.0d0)
      nbclocal(2) = (ccc*ccc)/(4.0d0*aaa)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 18.998d0
      elelocal =  ' F'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     F in HF [jacs_1984_CJ]
      nbnlocal = 'F'
      cbnlocal = 'F'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Hydrogen **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     H (COH) [jpc_1986_a]
c     also H in HF [jacs_1984_CJ]
      nbnlocal = 'H'
      cbnlocal = 'H'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Lithium **

      type = type + 1
      ii = twh_getnbtype(type,type)
      aaa = (0.4d0)*CON_KCALTOK
      ccc = (100.0d0)*CON_KCALTOK
      nbclocal(1) = (aaa/ccc)**(1.0d0/6.0d0)
      nbclocal(2) = (ccc*ccc)/(4.0d0*aaa)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 6.941d0
      elelocal =  'Li'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     Li+ ion [jacs_1984_CSJ]
      nbnlocal = 'Li+'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     ** Nitrogen **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     N(primary) formamide, Asn [jacs_1985_a]
c     N(secondary) NMA, Ala [jacs_1985_a]
c     N(primary) DMF, Pro [jacs_1985_a]
      nbnlocal = 'Nsp2'
      cbnlocal = 'Nsp2'
      canlocal = 'Nsp2'
      ctnlocal = 'Nsp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.200d0
c     parameters expressed in kj/mol in this paper unlike most
      nbclocal(2) = (0.711d0)*(CON_KJTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 14.007d0
      elelocal =  ' N'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     triple bonded Nitrogen (acetonitrile) [mp_1988_JB]
      nbnlocal = 'Ns'
      cbnlocal = 'Ns'
      canlocal = 'Ns'
      ctnlocal = 'Ns'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     ** Oxygen **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.070d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     O (COH) [jpc_1986_a]
      nbnlocal = 'Osp3'
      cbnlocal = 'Osp3'
      canlocal = 'Osp3'
      ctnlocal = 'Osp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.000d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     protonated oxygen in a carboxylic acid [jpc_1990_JBC]
c     oxygen bonded to 2 carbons in an ester [jpc_1990_JBC]
      nbnlocal = 'Osp3-COO'
      cbnlocal = 'Osp3-COO'
      canlocal = 'Osp3-COO'
      ctnlocal = 'Osp3-COO'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     O (C=O) amides [jacs_1985_a}
c     O (C=O) carboxylic acids [jpc_1990_JBC]
      nbnlocal = 'Osp2'
      cbnlocal = 'Osp2'
      canlocal = 'Osp2'
      ctnlocal = 'Osp2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Sodium **

      type = type + 1
      ii = twh_getnbtype(type,type)
      aaa = (14.0d0)*CON_KCALTOK
      ccc = (300.0d0)*CON_KCALTOK
      nbclocal(1) = (aaa/ccc)**(1.0d0/6.0d0)
      nbclocal(2) = (ccc*ccc)/(4.0d0*aaa)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 22.990d0
      elelocal =  'Na'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     Na+ ion [jacs_1984_CSJ]
      nbnlocal = 'Na+'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** Sulfur **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.70d0
      nbclocal(2) = (0.25d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 32.066d0
      elelocal =  ' S'
      bplocal = 'p'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     S (H2S) [jpc_1986_b]
      nbnlocal = 'S-HSH'
      cbnlocal = 'S'
      canlocal = 'S'
      ctnlocal = 'S'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 32.066d0
      elelocal =  ' S'
      bplocal = 'p'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     S (H2S) [jpc_1986_b]
      nbnlocal = 'S-RSR'
      cbnlocal = 'S'
      canlocal = 'S'
      ctnlocal = 'S'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed bond lengths                                             *
c     ******************************************************************

c     to try and make this easy to search the bond lengths are 
c     listed in alphabetical order 


c     ** Carbon - Carbon+ bond lengths

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.53d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     --- C4-C4 [jacs_1984_JMS] where this means a C sp3 - C sp3 bond
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Csp3'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.50d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     C4-C3 where this case is Csp3 - Csp2 [jacs_1984_JMS]
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Csp2'
c     C4-C3 where this case is Csp3 - Carom [jacs_1984_JMS]
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Carom'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.458d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     Csp3 - Csp from acetonitrile [mp_1988_JB]
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Csp'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.34d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     --- C3=C3 [jacs_1984_JMS] where this means a C sp2 - C sp2 
c     --- nonaromatic bond
      ii = ii + 1
      vnlocal(ii,1) = 'Csp2'
      vnlocal(ii,2) = 'Csp2'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.522d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     Csp2-N - Csp3 general in an amide/amino acid [jacs_1985_JS]
      ii = ii + 1
      vnlocal(ii,1) = 'Csp2-N'
      vnlocal(ii,2) = 'Csp3'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.40d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     --- C3=C3b [jacs_1984_JMS] where this means a C sp2 - C sp2 
c     --- aromatic bond
      ii = ii + 1
      vnlocal(ii,1) = 'Carom'
      vnlocal(ii,2) = 'Carom'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.335d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     Csp2 - Nsp2 [jacs_1985_JS] in amides and peptides
      ii = ii + 1
      vnlocal(ii,1) = 'Csp2-N'
      vnlocal(ii,2) = 'Nsp2'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.449d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     Csp3 - Nsp2 [jacs_1985_JS] in amides and peptides
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Nsp2'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.157d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     Csp - Ns bond in acetonitrile [mp_1988_JB]
      ii = ii + 1
      vnlocal(ii,1) = 'Csp'
      vnlocal(ii,2) = 'Ns'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.229d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     C=Osp2 [jacs_1985_JS] in amides
      ii = ii + 1
      vnlocal(ii,1) = 'Csp2-N'
      vnlocal(ii,2) = 'Osp2'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.430d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     C-Osp3 [jpc_1986_a] in aliphatic alcohols (primary)
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'Osp3'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.82d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     C-S [jpc_1986_b]
      ii = ii + 1
      vnlocal(ii,1) = 'Csp3'
      vnlocal(ii,2) = 'S'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ** Hydrogen - Hydrogen+ bond lengths **

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.917d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     H-F [jacs_1984_CJ]
      ii = ii + 1
      vnlocal(ii,1) = 'H'
      vnlocal(ii,2) = 'F'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.96d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     H-Nsp2 [jacs_1985_JS] in an amide
      ii = ii + 1
      vnlocal(ii,1) = 'H'
      vnlocal(ii,2) = 'Nsp2'
c     --- set nvsame 
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.945d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     H-O [jpc_1986_a] in an alcohol
      ii = ii + 1
      vnlocal(ii,1) = 'H'
      vnlocal(ii,2) = 'Osp3'
c     --- set nvsame 
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.34d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     H-S [jpc_1986_b]
      ii = ii + 1
      vnlocal(ii,1) = 'H'
      vnlocal(ii,2) = 'S'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ** Sulfur - Sulfur+ bond lengths

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 2.04d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     S-S [jpc_1986_b]
      ii = ii + 1
      vnlocal(ii,1) = 'S'
      vnlocal(ii,2) = 'S'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed angles                                                   *
c     ******************************************************************

c     ** arranged alphabetically by central atom **

c     ** Carbon central atom **

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 112.00d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     --- C-C4-C [jacs_1984_JMS] with a CH2sp3 as the center
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Cwild'
c     --- C-C4-C [jacs_1984_JMS] with CHsp3 as the center
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'Cwild'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 109.47d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     --- C-C4-C [jacs_1984_JMS] with Csp3 as the center, inferred from
c     --- geometrical considerations (must be tetrahedral).
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'Cwild'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 124.00d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     --- C-C3-C [jacs_1984_JMS] with CHsp2 as the center
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'CHsp2'
      bnlocal(ii,3) = 'Cwild'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 120.00d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     --- C-C3-C with Csp2 as the center from geometrical
c     --- considerations (must be 120.0)
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Csp2'
      bnlocal(ii,3) = 'Cwild'
c     --- C-C3-C with Carom as the center must be 120
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Carom'
      bnlocal(ii,3) = 'Cwild'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 108.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CCO [jpc_1986_a] aliphatic alcohols (primary)
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'Osp3'
c     CCO [jpc_1986_a] aliphatic alcohols (secondary)
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'CHsp3'
      bnlocal(ii,3) = 'Osp3'
c     CCO [jpc_1986_a] aliphatic alcohols (tertiary)
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Csp3'
      bnlocal(ii,3) = 'Osp3'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 116.6d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     C-C-N in an amide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Csp2-N'
      bnlocal(ii,3) = 'Nsp2'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 120.5d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     C-C-O in an amide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Csp2-N'
      bnlocal(ii,3) = 'Osp2'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 122.9d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     N-C-O in an amide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'Nsp2'
      bnlocal(ii,2) = 'Csp2-N'
      bnlocal(ii,3) = 'Osp2'
c     N-C-O in formamide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'Nsp2'
      bnlocal(ii,2) = 'CHsp2-N'
      bnlocal(ii,3) = 'Osp2'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 180.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     C-C#N in acetonitrile [mp_1988_JB]
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'Csp'
      bnlocal(ii,3) = 'Ns'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ** Nitrogen central atom **

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 121.9d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Csp2-N-C in an amide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'Csp2-N'
      bnlocal(ii,2) = 'Nsp2'
      bnlocal(ii,3) = 'Cwild'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 119.8d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Csp2-N-H in an amide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'Csp2-N'
      bnlocal(ii,2) = 'Nsp2'
      bnlocal(ii,3) = 'H'
c     CHsp2-N-H in formamide [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp2-N'
      bnlocal(ii,2) = 'Nsp2'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 118.3d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Csp3-N-H in an amide [jacs_1985_JS]
c     CH2sp3
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'Nsp2'
      bnlocal(ii,3) = 'H'
c     Csp3-N-H in an amide [jacs_1985_JS]
c     CHsp3
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'Nsp2'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 120.4d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     H-N-H in an amide inferred from [jacs_1985_JS]
      ii = ii + 1
      bnlocal(ii,1) = 'H'
      bnlocal(ii,2) = 'Nsp2'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ** Oxygen central atom **

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 108.5d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      ii = 0
c     COH [jpc_1986_a] aliphatic alcohols (primary)
      ii = ii + 1
      bnlocal(ii,1) = 'CH2sp3'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'H'
c     COH [jpc_1986_a] aliphatic alcohols (secondary)
      ii = ii + 1
      bnlocal(ii,1) = 'CHsp3'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'H'
c     COH [jpc_1986_a] aliphatic alcohols (tertiary)
      ii = ii + 1
      bnlocal(ii,1) = 'Csp3'
      bnlocal(ii,2) = 'Osp3'
      bnlocal(ii,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 114.4d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CCS [jpc_1986_b] (methyl groups)
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'C32sp3'
      bnlocal(ii,3) = 'S'
c     CCS [jpc_1986_b] (methylene groups)
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'CH2sp3'
      bnlocal(ii,3) = 'S'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ** sulfur central atom **

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 92.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     HSH [jpc_1986_b]
      ii = ii + 1
      bnlocal(ii,1) = 'H'
      bnlocal(ii,2) = 'S'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 96.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CSH [jpc_1986_b]
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'S'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 98.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     SSH [jpc_1986_b]
      ii = ii + 1
      bnlocal(ii,1) = 'S'
      bnlocal(ii,2) = 'S'
      bnlocal(ii,3) = 'H'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 99.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CSC [jpc_1986_b]
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'S'
      bnlocal(ii,3) = 'Cwild'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 103.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     CSS [jpc_1986_b]
      ii = ii + 1
      bnlocal(ii,1) = 'Cwild'
      bnlocal(ii,2) = 'S'
      bnlocal(ii,3) = 'S'
c     --- set nbname
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * OPLS-AA torsions using standard OPLS cosine series             *
c     * k0 + k1*(1 + Cos(phi)) + k2*(1 - Cos(2*phi))                   *
c     * + k3*(1 + Cos(3*phi)                                           *
c     ******************************************************************
c     OPLS has a list of hydrocarbons that give the atoms and an example
c     molecule.  Thus, I have had to take these examples and translate
c     them into appropriate atom types.  Hopefully this all worked out
c     properly.  The OPLS notation is at the top of each torsion, along
c     with the example system.  M.G. Martin (07-05-2001)

c     ** Torsions using the OPLS-UA torsional style **
c     ** This is the majority of the torsions **

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (1.363d0)*CON_KCALTOK
      tclocal(1) =  (0.343d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-0.436d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (-1.121d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- 1-butene [jacs_1984_JMS]
      tnlocal(1,1) = 'CH2sp2'
      tnlocal(1,2) = 'CHsp2'
      tnlocal(1,3) = 'CH2sp3'
      tnlocal(1,4) = 'CH3sp3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (2.713d0)*CON_KCALTOK
      tclocal(1) =  (1.526d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (0.533d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (-3.453d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- isopentane [jacs_1984_JMS]
      tnlocal(1,1) = 'CH3sp3'
      tnlocal(1,2) = 'CHsp3'
      tnlocal(1,3) = 'CH2sp3'
      tnlocal(1,4) = 'CH3sp3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.0d0)*CON_KCALTOK
      tclocal(1) =  (1.522d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-0.315d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (3.207d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- n-butane [jacs_1984_JMS]
      tnlocal(1,1) = 'CH3sp3'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'CH2sp3'
      tnlocal(1,4) = 'CH3sp3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.0d0)*CON_KCALTOK
      tclocal(1) =  (1.411d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-0.271d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (3.145d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- other n-alkanes [jacs_1984_JMS] including one CH3 end
      tnlocal(1,1) = 'CH3sp3'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'CH2sp3'
      tnlocal(1,4) = 'CH2sp3'
c     --- other n-alkanes [jacs_1984_JMS] including no CH3 ends
      tnlocal(2,1) = 'CH2sp3'
      tnlocal(2,2) = 'CH2sp3'
      tnlocal(2,3) = 'CH2sp3'
      tnlocal(2,4) = 'CH2sp3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.0d0)*CON_KCALTOK
      tclocal(1) =  (0.834d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-0.116d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (0.747d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     ethanol/propanol C-O bond central [jpc_1986_a]
c     assuming this holds for longer primary alkanols as well
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'CH2sp3'
      tnlocal(ii,3) = 'Osp3'
      tnlocal(ii,4) = 'H'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.0d0)*CON_KCALTOK
      tclocal(1) =  (0.702d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-0.212d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (3.060d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     1-propanol C1-C2 bond central [jpc_1986_a]
c     assuming this is true for longer primary alkanols as well
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'CH2sp3'
      tnlocal(1,4) = 'Osp3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.429d0)*CON_KCALTOK
      tclocal(1) =  (0.784d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (0.125d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (-0.691d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     2-propanol C-O bond central [jpc_1986_a]
c     assuming this is true for longer secondary alkanols as well
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CHsp3'
      tnlocal(1,3) = 'Osp3'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.0d0)*CON_KCALTOK
      tclocal(1) =  (0.0d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (0.0d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (0.650d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     2-methyl-2-propanol C-O bond central [jpc_1986_a]
c     --- assuming this is true for longer tertiary alkanols as well
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'Csp3'
      tnlocal(1,3) = 'Osp3'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.2100d0)*CON_KCALTOK
      tclocal(1) =  (-0.2067d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-0.0733d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (1.2167d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CH3-CH2-S-H, C-S bond[jpc_1986_b]
c     assuming this is true for longer alkyl chains as well
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (0.00d0)*CON_KCALTOK
      tclocal(1) =  (1.4600d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-1.0730d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (2.3100d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CH3-CH2-S-CH3; C-S bond [jpc_1986_b]
c     CH3-CH2-S-CH2; C-S bond [jpc_1986_b]
c     assuming this is true for longer alkyl tails as well
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CH2sp3'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'Cwild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,11)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = (7.0100d0)*CON_KCALTOK
      tclocal(1) =  (1.5820d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (-8.6700d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) =  (2.0080d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CH3-S-S-CH3; S-S bond [jpc_1986_b]
c     assuming this is true for longer alkyl tails as well
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'Cwild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     ** Torsions that use the 2 term cosine series  **
c     ** this is much less common for this forcefield **

      type = type + 1
      call twh_torstyle(GLB_SET,type,13)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(1) =  (2.8d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (21.2d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     Cwild - Co - N - Cwild bond in amides [jacs_1985_JS]
      ii = ii + 1
      tnlocal(ii,1) = 'Cwild'
      tnlocal(ii,2) = 'Csp2-N'
      tnlocal(ii,3) = 'Nsp2'
      tnlocal(ii,4) = 'Cwild'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,13)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(1) =  (0.0d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) =  (20.0d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ii = 0
c     wild - Co - N - H bond in amides [jacs_1985_JS]
      ii = ii + 1
      tnlocal(ii,1) = 'wild'
      tnlocal(ii,2) = 'Csp2-N'
      tnlocal(ii,3) = 'Nsp2'
      tnlocal(ii,4) = 'H'
c     wild - Co - N - H bond in formamide [jacs_1985_JS]
      ii = ii + 1
      tnlocal(ii,1) = 'wild'
      tnlocal(ii,2) = 'CHsp2-N'
      tnlocal(ii,3) = 'Nsp2'
      tnlocal(ii,4) = 'H'
c     --- set ntsame
      call twh_set_tornames(type,ii,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
