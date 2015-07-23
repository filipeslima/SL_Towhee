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
      subroutine twh_ffoplsaa(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the oplsaa       *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
c     * Added additional OPLS improper torsion names AAC 9-27-05       *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      integer ii,type
      character*2 c2temp
      character*3 c3temp
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
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)

      ffname = 'OPLS-aa'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Geometric')
      type = 0

c     **************************************************************
c     * OPLS-AA force field                                        *
c     **************************************************************
c     Reference list
c     [jacs1984] J. Chandrasekhar, D. Spellmeyer, W.L. Jorgensen
c     `Energy component analysis for dilute aqueous 
c     solutions of Li+, Na+, F- and Cl- Ions'
c     J. Am. Chem. Soc. 106, 903-910 (1984)
c     [jacs1991] J. Pranata, S.G. Wierschke, W.L. Jorgensen
c     `OPLS potential functions for nucleotide bases.
c     Relatvie association constants of hydrogen-bonded 
c     base pairs in chloroform' J. Am. Chem. Soc.
c     113, 2810-2819 (1991)         
c     [jpc1994]  G. Kaminski, E.M. Duffy, T. Matsui, W.L. Jorgensen
c     `Free energies of hydration and pure liquid 
c     properties of hydrocarbons from the OPLS all-atom 
c     model', J. Phys. Chem. 98, 13077-13082 (1994)
c     [jacs1996] W.L. Jorgensen, D.S. Maxwell, J. Tirado-Rives 
c     `Development and testing of the OPLS all-atom force
c     field on conformational energetics and properties 
c     of organic liquids' J. Am. Chem. Soc. 
c     118, 11225-11236 (1996) {see the sup. info.}
c     [wlj1996]  `OPLS all-atom parameters for organic molecules,
c     ions, & nucleic acids 12/96' available by sending 
c     email to W.L. Jorgensen (bill@adrik.chem.yale.edu)
c     [jcc1997]  W. Damm, A. Frontera, J. Tirado-Rives, 
c     W.L. Jorgensen
c     `OPLS all-atom force field for carbohydrates'
c     [jms1998]  W.L. Jorgensen, N.A. McDonald `Development of an 
c     all-atom force field for heterocycles.  Properties 
c     of liquid pyridine and diazenes' J. Mol. Structure 
c     Theochem, 424, 145-155 (1998)
c     [jpcb1998] N.A. McDonald, W.L. Jorgensen `Development of an
c     all-atom force field for heterocycles.  Properties
c     of liquid pyrrole, furan, diazoles, and oxazoles'
c     J. Phys. Chem. B 102, 8049-8059 (1998)
c     [jacs1998] N.A. McDonald, E.M. Duffy, W.L. Jorgensen
c     `Monte Carlo investigations of selective anion
c     complexation by a 
c     bis(phenylurea)p-tert-butylcalix[4]arene'
c     J. Am. Chem. Soc. 120, 5104-5111 (1998)
c     [jacs1999] R.C. Rizzo, W.L. Jorgensen `OPLS all-atom model for
c     amines: resolution of the amine hydration problem'
c     J. Am. Chem. Soc. 121, 4827-4836 (1999)
c     [jcp2000]  M.W. Mahoney, W.L. Jorgensen
c     `A five-site model for liquid water and the 
c     reproduction of the density anomoly by rigid,
c     nonpolarizable potential functions' J. Chem. Phys.
c     112, 8910-8922 (2000)
c     [wlj2001]  `OPLS all-atom parameters for organic molecules,
c     ions, & nucleic acids 05/01' available by sending 
c     email to W.L. Jorgensen (bill@adrik.chem.yale.edu)
c     OPLS does not generally specify parameter names that are
c     --- useful for nonbond, bond, etc.  The code I use here is the
c     --- Amber types (as opls often uses the amber bonded terms)
c     --- followed by a -# if there are different nonbonded terms 
c     --- that have the same amber type.  OPLS has defined a few 
c     --- new atom types so where that is done I use those instead.
c     Be careful with the codes as the first two letters are used in
c     --- the ffbond and ffangle routines to determine types, while
c     --- the first three characters are used to determine torsions
c     --- Where a third character is used the following applies
c     --- ~ is for a neutral atom
c     --- + is for a positively charged ion
c     --- x is for an atom which is special for rna/dna
c     --- f is for atoms which are bonded to F
c     --- m is for O in an mono-alcohol
c     --- d is for O in a diol
c     --- p is for O in a polyol
c     --- a is for O in an acid
c     --- k is for a C in a ketone
c     --- j is for a CA at the junction of rings
c     To assign the proper charges for a molecule you will need 
c     --- to look in the appropriate OPLS references.


c     (Br-) Bromide ion
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.5950d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0900d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 79.904d0
      elelocal = 'Br'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Br-'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C k) Carbon in C=O not bonded to N [jacs1996]
      nbnlocal ='C k'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C a) Carbon in O-C=O acid [jacs1996]
      nbnlocal ='C a'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C n) Carbon in N-C=O [jacs1996] 
      nbnlocal ='C n'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CA~) neutral aromatic carbon [jms1998]
      nbnlocal ='CA~'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CN) sp2 aromatic CE in TRP [wlj2001]
      nbnlocal = 'CN'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CQ) sp2 C in 6-mem ring between lone pair N's [jms1998]
      nbnlocal ='CQ'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CW) sp2 aromatic C in 5-mem ring next to C & NH [jpcb1998]
      nbnlocal ='CW'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CS) arom C which does not have 120 bond angles with all
c     ---- of its neighbors.  This is a newly defined type and
c     ---- does not appear in amber. [jpcb1998]
      nbnlocal ='CS'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CU) arom C which abuts NB in a hteroatom-heteroatom bond
c     ---- new type for OPLS, not in amber [jpcb1998]
      nbnlocal ='CU'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CR) arom C in 5-mem ring next to two Nitrogens [jpcb1998]
      nbnlocal ='CR'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CV) arom C in 5-mem ring next to C and lone pair N [jpcb1998]
      nbnlocal ='CV'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (C*) arom C in 5-mem ring next to 2 carbons [jacs1996]
      nbnlocal ='C*'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CAj) neutral aromatic carbon at the junction of a 5-mem and 
c     ----- 6-mem rings [mgm]
      nbnlocal = 'CAj'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CA+) aromatic carbon in guanidinium C+ [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.250
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CA+'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Carbon parameters for the nucleotide bases
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CAxna) aromatic C in DAP C2,C3,C4 pyridine, 
c     ------ Cytosine C4,C5,C6
      nbnlocal ='CAxna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CBxna) Adenine C4,C5,C6 : Guanine C2,C4,C5
      nbnlocal ='CBxna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CKxna) Adenine C8 Guanine
      nbnlocal ='CKxna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CTxna) Thymine C-C5, 9-Me A or G C-N9
      nbnlocal ='CTxna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CMxna) Carbon in Uracil C5,C6
      nbnlocal ='CMxna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CM) alkene sp2 carbon (non-aromatic) [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.076d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CM'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CT) aliphatic sp3 hybrid carbon [jpc1994]
      nbnlocal ='CT'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (CO) sp3 anomeric carbon (bonded to ether O and alcohol O)
c     ---- this type shows up in carbohydrates [jcc1997]
      nbnlocal ='CO'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CTf) aliphatic sp3 hybrid carbon bonded to F [wlj1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.062d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CTf'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CB) arom C at junction of 5-mem and 6-mem rings [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CB'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Cl-) Chlorine (ionic -1) [jacs1984]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.41724
      nbclocal(2) = (0.11779d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 35.453d0
      elelocal = 'Cl'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl-'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Cl) Chlorine (nonionic, bonded) [wlj1996]
c     --- note next to it says chloroalkene, might be some newer
c     --- version of this in BOSS, but not yet generally available
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400
      nbclocal(2) = (0.300d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 35.453d0
      elelocal = 'Cl'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (F) Fluorine nonionic [wlj1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.94
      nbclocal(2) = (0.061d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (F-) Fluorine (ionic -1) [jacs1984]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.73295
      nbclocal(2) = (0.720d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F-'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (H) H attached to N [jpcb1998]
      nbnlocal ='H'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HO) H attached to O in alcohol [jacs1996]
      nbnlocal ='HO'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HS) H attached to S [jacs1996]
      nbnlocal ='HS'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (HC~1) hydrgen attached to an sp3 carbon in most cases
c     --- see special cases of hydrogen below [jpc1994]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HC~1'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.42
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HA) H bonded to an aromatic ring [jms1998]
      nbnlocal ='HC~2'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.42
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HC~2) H bonded to a non-aromatic sp2 carbon which is not 
c     ----- bonded to a N, not double bonded to O, and is not bonded
c     ----- to C=O [jacs1996]
      nbnlocal ='HA'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HC~3) H bonded to a C=O or bonded to a C which is bonded to
c     ------ a C=O [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.42
      nbclocal(2) = (0.015d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HC~3'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HC~4) H bonded to CT bonded to NT [jacs1999]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HC~4'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     Hydrogens bonded to Carbon on the nucleotide bases [jacs1991]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (HCxna) DAP H3,H4 : Uracil H-C5,H-C6 : Thymine H-CC5
c     --- Cytosine H-C5,H-C6 : Adenine H-C2,H-C8 Guanine
      nbnlocal ='HCxna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (H2xna) DAP H-amine, Cytosine, Adenine
      nbnlocal ='H2xna'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (I-) Iodine (ionic -1) [jacs1998]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 5.400
      nbclocal(2) = (0.0700d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal =  126.90d0
      elelocal = ' I'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='I-'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Li+) Lithium (ionic +1) [jacs1984]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.25992
      nbclocal(2) = (6.25d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 6.941d0
      elelocal = 'Li'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Li+'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NA) sp2 aromatic N with H attached [jpcb1998]
      nbnlocal ='NA'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NB) sp2 N in 5-mem ring with lone pair [jpcb1998]
      nbnlocal ='NB'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NC) sp2 N in 6-mem ring with lone pair [jms1998] 
      nbnlocal ='NC'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (N2) sp2 N of aromatic amines and guanidinium ions [jacs1991]
      nbnlocal ='N2'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (N) N in an amide [wlj1996]
      nbnlocal ='N'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (N3~) sp3 neutral N of amines [jacs1996] NOTE the N3 
c     ------- parameters are obsoleted by NT parameters [jacs1999]
      nbnlocal ='N3~'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (N3+) sp3 +1 N of ammonium ions [jacs1996] NOTE the N3 
c     ------- parameters are obsoleted by NT parameters [jacs1999]
      nbnlocal ='N3+'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.42
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NT~1) sp3 N in ammonia 
      nbnlocal ='NT~1'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.42
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NT+1) sp3 N in ammonium ion.  
c     ------ This is new OPLS type that replaces
c     ------ the old N3 in Amber.  Note that N3 parameters were 
c     ------ given in [jacs1996] and these are now obsolete with the
c     ------ publication of [jacs1999]
      nbnlocal ='NT+1'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.30
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NT~2) sp3 N in primary,secondary, tertiary amine 
      nbnlocal ='NT~2'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.30
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (NT+2) sp3 N in primary,secondary, tertiary ammonium ion
c     ------  This is a new OPLS type that replaces
c     ------ the old N3 in Amber.  Note that N3 parameters were 
c     ------ given in [jacs1996] and these are now obsolete with the
c     ------ publication of [jacs1999]         
      nbnlocal ='NT+2'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25
      nbclocal(2) = (0.120d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     --- NO (N in nitro R-NO2) [wlj2001]
      nbnlocal ='NO'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Na+) Sodium (ionic +1) [jacs1984]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.89744
      nbclocal(2) = (1.60714d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 22.990d0
      elelocal = 'Na'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na+'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (O) Oxygen in C=O [jacs1996]
      nbnlocal ='O'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (O2) Oxygen in COO- or COOH [jacs1996]
      nbnlocal ='O2'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OHm) Oxygen in a mono-alcohol [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.120
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OHm'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OHp) O in polyols or phenol [jacs1996]  
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.07
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OHp'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OHa) O in RCOOH [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.000
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OHa'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     --- (ON) O in nitro R-NO2 [wlj2001]
      nbnlocal ='ON'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (OS) sp3 O in an ether [jpcb1998] or acetal [jacs1996]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.90
      nbclocal(2) = (0.140d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OS'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (SH) S in thiols [jacs1996]
      nbnlocal ='SH'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (S) Sulfur in sulfide, disulfide [jacs1996]
      nbnlocal ='S'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Amber param96 bond potentials                                  *
c     ******************************************************************
c     -- W.D. Cornell et al. J. Am. Chem. Soc. (1995), 117, 5179-5197

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.409d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - CA [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.419d0
      vclocal(1) = (447.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C - CB [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CB'
c     --- CB - CN [wlj1996]
      vnlocal(2,1) = 'CB'
      vnlocal(2,2) = 'CN'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.444d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - CM [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - CT [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (490.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - N [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'N'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.383d0
      vclocal(1) = (424.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - N* [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.388d0
      vclocal(1) = (418.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - NA [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.358d0
      vclocal(1) = (457.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - NC [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.229d0
      vclocal(1) = (570.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - O [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.250d0
      vclocal(1) = (656.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - O2 [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.364d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C - OH [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'OH'
c     --- CA - OH from GROMACS OPLS-aa AAC
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'OH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.459d0
      vclocal(1) = (388.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C* - CB [wlj1996]
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'CB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.495d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C* - CT [wlj1996]
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.352d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C* - CC [wlj1996]
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'CC'
c     --- C* - CW [wlj1996]
      vnlocal(2,1) = 'C*'
      vnlocal(2,2) = 'CW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.080d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 7
c     --- CA - HA [jms1998]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'HA'
c     --- CQ - HA  [jms1998]
      vnlocal(2,1) = 'CQ'
      vnlocal(2,2) = 'HA'
c     --- HA - CW [jpcb1998]
      vnlocal(3,1) = 'HA'
      vnlocal(3,2) = 'CW'
c     --- HA - CS [jpcb1998]
      vnlocal(4,1) = 'HA'
      vnlocal(4,2) = 'CS'
c     --- HA - CU [jpcb1998]
      vnlocal(5,1) = 'HA'
      vnlocal(5,2) = 'CU'
c     --- HA - CR [jpcb1998]
      vnlocal(6,1) = 'HA'
      vnlocal(6,2) = 'CR'
c     --- HA - CV [jpcb1998]
      vnlocal(7,1) = 'HA'
      vnlocal(7,2) = 'CV'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.400d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 8
c     --- C - CA [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CA'
c     --- C - CD [wlj1996]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'CD'
c     --- CA - CA [wlj1996]
      vnlocal(3,1) = 'CA'
      vnlocal(3,2) = 'CA'
c     --- CA - CD [wlj1996]
      vnlocal(4,1) = 'CA'
      vnlocal(4,2) = 'CD'
c     --- CA - CN [wlj1996]
      vnlocal(5,1) = 'CA'
      vnlocal(5,2) = 'CN'
c     --- CB - CD [wlj1996]
      vnlocal(6,1) = 'CB'
      vnlocal(6,2) = 'CD'
c     --- CA - CC [wlj1996]
      vnlocal(7,1) = 'CA'
      vnlocal(7,2) = 'CC'
c     --- CD - CC [wlj1996]
      vnlocal(8,1) = 'CD'
      vnlocal(8,2) = 'CC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.404d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - CB
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.433d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - CM [wlj1996]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.510d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CA - CT [wlj1996]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CT'
c     --- CM - CT [wlj1996]
      vnlocal(2,1) = 'CM'
      vnlocal(2,2) = 'CT'
c     --- CB - CT [wlj1996]
      vnlocal(3,1) = 'CB'
      vnlocal(3,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (481.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CA - NT [jacs1999]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'NT'
c     --- CA - N2 [wlj1996]
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'N2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- NA - CW [jpc1998]
      vnlocal(1,1) = 'NA'
      vnlocal(1,2) = 'CW'
c     --- CA - NA [wlj1996]
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'NA'
c     --- CA - N [wlj1996]
      vnlocal(3,1) = 'CA'
      vnlocal(3,2) = 'N'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.339d0
      vclocal(1) = (483.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - NC [wlj1996]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.370d0
      vclocal(1) = (520.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CV - CW [jpcb1998]
      vnlocal(1,1) = 'CV'
      vnlocal(1,2) = 'CW'
c     --- CB - CB [wlj1996]
      vnlocal(2,1) = 'CB'
      vnlocal(2,2) = 'CB'
c     --- CC - CB [wlj1996]
      vnlocal(3,1) = 'CC'
      vnlocal(3,2) = 'CB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.374d0
      vclocal(1) = (436.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CB - N* [wlj1996]
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'N*'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.391d0
      vclocal(1) = (414.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CB - NB [wlj1996]
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'NB'
c     --- CA - NB [wlj1996]
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.354d0
      vclocal(1) = (461.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CB - NC [wlj1996]
      vnlocal(1,1) = 'CB'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.504d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CC - CT [wlj1996]
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CT'
c     --- CV - CT from GROMACS opls-aa AAC
      vnlocal(2,1) = 'CV'
      vnlocal(2,2) = 'CT'
c     --- CW - CT from TINKER -AAC
      vnlocal(3,1) = 'CW'
      vnlocal(3,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (512.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CC - CV [wlj1996]
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CV'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (518.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CC - CW [wlj1996]
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CW'
c     --- CW - CW from TINKER - aac
      vnlocal(2,1) = 'CW'
      vnlocal(2,2) = 'CW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.385d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CC - NA [wlj1996]
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.394d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- NB - CV [jpcb1998]
      vnlocal(1,1) = 'NB'
      vnlocal(1,2) = 'CV'
c     --- CC - NB [jpcb1998]
      vnlocal(2,1) = 'CC'
      vnlocal(2,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CK - N*
      vnlocal(1,1) = 'CK'
      vnlocal(1,2) = 'N*'
c     --- CK - NA
      vnlocal(2,1) = 'CK'
      vnlocal(2,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.304d0
      vclocal(1) = (529.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CK - NB [wlj1996]
      vnlocal(1,1) = 'CK'
      vnlocal(1,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (448.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CM - N* [wlj1996]
      vnlocal(1,1) = 'CM'
      vnlocal(1,2) = 'N*'
c     --- CM - NA [wlj1996]
      vnlocal(2,1) = 'CM'
      vnlocal(2,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.380d0
      vclocal(1) = (428.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CN - NA [wlj1996]
      vnlocal(1,1) = 'CN'
      vnlocal(1,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.324d0
      vclocal(1) = (502.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CQ - NC [jms1998]
      vnlocal(1,1) = 'CQ'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (477.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- NA - CR [jpcb1998]
      vnlocal(1,1) = 'NA'
      vnlocal(1,2) = 'CR'
c     --- CP - NA [wlj1996]
      vnlocal(2,1) = 'CP'
      vnlocal(2,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (488.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- NB - CR [jpcb1998]
      vnlocal(1,1) = 'NB'
      vnlocal(1,2) = 'CR'
c     --- CP - NB [wlj1996]
      vnlocal(2,1) = 'CP'
      vnlocal(2,2) = 'NB'
c     --- CR - NB [wlj1996]
      vnlocal(3,1) = 'CR'
      vnlocal(3,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.090d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CT - HC [jacs1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'HC'
c     --- CO - HC [jcc1997]
      vnlocal(2,1) = 'CO'
      vnlocal(2,2) = 'HC'
c     --- HC - C [jcc1997]
      vnlocal(3,1) = 'HC'
      vnlocal(3,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - N [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT - N* [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N*'
c     --- CT - NA [wlj1996]
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'NA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.463d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - N2
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - N3 [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'N3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.410d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT - OH [jcc1997]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'OH'
c     --- CT - OS [jcc1997]
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.010d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 4
c     --- H - NA [jpc1998]
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'NA'
c     --- H - NT [jacs1999]
      vnlocal(2,1) = 'H'
      vnlocal(2,2) = 'NT'
c     --- H - N [wlj1996]
      vnlocal(3,1) = 'H'
      vnlocal(3,2) = 'N'
c     --- H - N2 [wlj1996]
      vnlocal(4,1) = 'H'
      vnlocal(4,2) = 'N2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.336d0
      vclocal(1) = (274.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HS - SH [wlj1996]
      vnlocal(1,1) = 'HS'
      vnlocal(1,2) = 'SH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.480d0
      vclocal(1) = (525.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- O2 - P [wlj1996]
      vnlocal(1,1) = 'O2'
      vnlocal(1,2) = 'P'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.610d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- OH - P [wlj1996]
      vnlocal(1,1) = 'OH'
      vnlocal(1,2) = 'P'
c     --- OS - P [wlj1996]
      vnlocal(2,1) = 'OS'
      vnlocal(2,2) = 'P'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.038d0
      vclocal(1) = (166.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- S - S [wlj1996]
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.529d0
      vclocal(1) = (268.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT - CT [jacs1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'CT'
c     --- CT - CO [jcc1997]
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'CO'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.320d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- NC - NC [jms1998]
      vnlocal(1,1) = 'NC'
      vnlocal(1,2) = 'NC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.367d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CS - CW [jpcb1998]
      vnlocal(1,1) = 'CS'
      vnlocal(1,2) = 'CW'
c     --- CA - CW [jpcb1998]
      vnlocal(2,1) = 'CA'
      vnlocal(2,2) = 'CW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.424d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CS - CS [jpcb1998]
      vnlocal(1,1) = 'CS'
      vnlocal(1,2) = 'CS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.349d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- NA - NB [jpcb1998]
      vnlocal(1,1) = 'NA'
      vnlocal(1,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.320d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- NB - CU [jpcb1998]
      vnlocal(1,1) = 'NB'
      vnlocal(1,2) = 'CU'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.421d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CU - CA [jpcb1998]
      vnlocal(1,1) = 'CU'
      vnlocal(1,2) = 'CA'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.360d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CW - OS [jpcb1998]
      vnlocal(1,1) = 'CW'
      vnlocal(1,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.399d0
      vclocal(1) = (462.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OS - NB [jpcb1998]
      vnlocal(1,1) = 'OS'
      vnlocal(1,2) = 'NB'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.357d0
      vclocal(1) = (462.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OS - CR [jpcb1998]
      vnlocal(1,1) = 'OS'
      vnlocal(1,2) = 'CR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.448d0
      vclocal(1) = (382.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - NT [jacs1999]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'NT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.380d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CO - OS  [jcc1997]
      vnlocal(1,1) = 'CO'
      vnlocal(1,2) = 'OS'
c     --- CO - OH  [jcc1997]
      vnlocal(2,1) = 'CO'
      vnlocal(2,2) = 'OH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.945d0
      vclocal(1) = (553.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- HO - OH  [jcc1997]
      vnlocal(1,1) = 'HO'
      vnlocal(1,2) = 'OH'
c     --- HO - OS  [wlj1996]
      vnlocal(2,1) = 'HO'
      vnlocal(2,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)
       
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.080d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- C* - HC  [wlj1996]
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'HC'
c     --- CK - HA  [wlj1996]
      vnlocal(2,1) = 'CK'
      vnlocal(2,2) = 'HA'
c     --- CM - HC  [wlj1996]
      vnlocal(3,1) = 'CM'
      vnlocal(3,2) = 'HC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.327d0
      vclocal(1) = (214.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - OS  [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.382d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - CW  [wlj1996]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CM - CM  [wlj1996]
      vnlocal(1,1) = 'CM'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.810d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CT - S  [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'S'
c     --- CT - SH  [wlj1996]
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'SH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.781d0
      vclocal(1) = (245.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - Cl  [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'Cl'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.843d0
      vclocal(1) = (212.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - P  [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'P'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.332d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - F  [wlj1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'F'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.225d0
      vclocal(1) = (550.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- NO - ON [wlj2001]
      vnlocal(1,1) = 'NO'
      vnlocal(1,2) = 'ON'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.490d0
      vclocal(1) = (375.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - NO [wlj2001]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'NO'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.460d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CA - NO [wlj2001]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'NO'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CA - CA - CA [jms1998]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CA'
c     --- CA - CA - CB [wlj1996]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 19
c     --- CA - CA - HA [jms1998]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'HA'
c     --- H - NA - CW [jpcb1998]
      bnlocal(2,1) = 'H'
      bnlocal(2,2) = 'NA'
      bnlocal(2,3) = 'CW'
c     --- H - NA - CR [jpcb1998]
      bnlocal(3,1) = 'H'
      bnlocal(3,2) = 'NA'
      bnlocal(3,3) = 'CR'
c     --- HA - CR - NA [jpcb1998]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CR'
      bnlocal(4,3) = 'NA'
c     --- HA - CR - NB [jpcb1998]
      bnlocal(5,1) = 'HA'
      bnlocal(5,2) = 'CR'
      bnlocal(5,3) = 'NB'
c     --- HA - CV - NB [jpcb1998]
      bnlocal(6,1) = 'HA'
      bnlocal(6,2) = 'CV'
      bnlocal(6,3) = 'NB'
c     --- C - CA - HA [wlj1996]
      bnlocal(7,1) = 'C'
      bnlocal(7,2) = 'CA'
      bnlocal(7,3) = 'HA'
c     --- CB - CA - HA [wlj1996]
      bnlocal(8,1) = 'CB'
      bnlocal(8,2) = 'CA'
      bnlocal(8,3) = 'HA'
c     --- CN - CA - HA [wlj1996]
      bnlocal(9,1) = 'CN'
      bnlocal(9,2) = 'CA'
      bnlocal(9,3) = 'HA'
c     --- CM - CM - HC [wlj1996]
      bnlocal(10,1) = 'CM'
      bnlocal(10,2) = 'CM'
      bnlocal(10,3) = 'HC'
c     --- O - C - HC [wlj1996]
      bnlocal(11,1) = 'O'
      bnlocal(11,2) = 'C'
      bnlocal(11,3) = 'HC'
c     --- HA - CK - NA [wlj1996]
      bnlocal(12,1) = 'HA'
      bnlocal(12,2) = 'CK'
      bnlocal(12,3) = 'NA'
c     --- HA - CK - NB [wlj1996]
      bnlocal(13,1) = 'HA'
      bnlocal(13,2) = 'CK'
      bnlocal(13,3) = 'NB'
c     --- C* - CW - HA [wlj1996]
      bnlocal(14,1) = 'C*'
      bnlocal(14,2) = 'CW'
      bnlocal(14,3) = 'HA'
c     --- HA - CA - NA [wlj1996]
      bnlocal(15,1) = 'HA'
      bnlocal(15,2) = 'CA'
      bnlocal(15,3) = 'NA'
c     --- H - N - H [wlj1996]
      bnlocal(16,1) = 'H'
      bnlocal(16,2) = 'N'
      bnlocal(16,3) = 'H'
c     --- CA - N2 - H [wlj1996]
      bnlocal(17,1) = 'CA'
      bnlocal(17,2) = 'N2'
      bnlocal(17,3) = 'H'
c     --- CA - NA - H [wlj1996]
      bnlocal(18,1) = 'CA'
      bnlocal(18,2) = 'NA'
      bnlocal(18,3) = 'H'
c     --- H - N2 - H from GROMACS AAC
      bnlocal(19,1) = 'H'
      bnlocal(19,2) = 'N2'
      bnlocal(19,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - CB - NB [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NB'
c     --- CT - CM - CT [wlj1996]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CM - CT [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.70d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CM - HC [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.10d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - CT - CT [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
c     --- CM - CT - CT [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'CT'
c     --- CT - CT - NO (wlj nitro) [wlj2001]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'NO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 9
c     --- CT - CT - OH [jcc1997]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'OH'
c     --- CO - CT - OH [jcc1997]
      bnlocal(2,1) = 'CO'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'OH'
c     --- CT - CT - OS [jcc1997]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'OS'
c     --- CT - CO - OS [jcc1997]
      bnlocal(4,1) = 'CT'
      bnlocal(4,2) = 'CO'
      bnlocal(4,3) = 'OS'
c     --- CT - CO - OH [jcc1997]
      bnlocal(5,1) = 'CT'
      bnlocal(5,2) = 'CO'
      bnlocal(5,3) = 'OH'
c     --- CA - NT - CT [jacs1999]
      bnlocal(6,1) = 'CA'
      bnlocal(6,2) = 'NT'
      bnlocal(6,3) = 'CT'
c     --- CT - CT - N* [wlj1996]
      bnlocal(7,1) = 'CT'
      bnlocal(7,2) = 'CT'
      bnlocal(7,3) = 'N*'
c     --- N* - CT - OS [wlj1996]
      bnlocal(8,1) = 'N*'
      bnlocal(8,2) = 'CT'
      bnlocal(8,3) = 'OS'
c     --- CT - CT - F [wlj1996]
      bnlocal(9,1) = 'CT'
      bnlocal(9,2) = 'CT'
      bnlocal(9,3) = 'F'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.10d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CT - N [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.20d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 6
c     --- C - CT - N2 [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N2'
c     --- C - CT - N3 [wlj1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'N3'
c     --- C - CT - NT [wlj1996]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'NT'
c     --- CT - CT - N3 [wlj1996]
      bnlocal(4,1) = 'CT'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'N3'
c     --- CT - CT - NT [wlj1996]
      bnlocal(5,1) = 'CT'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'NT'
c     --- CT - CT - N2 from GROMACS AAC [JCP 76 1439]
      bnlocal(6,1) = 'CT'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'N2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.90d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - N - CT [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CT - CW - NA [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'NA'
c     --- C - N* - CM [wlj1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'N*'
      bnlocal(2,3) = 'CM'
c     --- C - NA - CM [wlj1996]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'NA'
      bnlocal(3,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - N* - CT [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NA - C [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - NA - CA [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'CA'
c     --- CA - NA - CA [wlj1996]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'NA'
      bnlocal(2,3) = 'CA'
c     --- CA - NB - CA [wlj1996]
      bnlocal(3,1) = 'CA'
      bnlocal(3,2) = 'NB'
      bnlocal(3,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NC - CA [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - OH - HO [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'HO'
c     --- CA - OH - HO from GROMACS - AAC
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'OH'
      bnlocal(2,3) = 'HO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.60d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C* - CT - CT [wlj1996]
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- CA - CA - NA [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NA'
c     --- C* - CW - NA [wlj1996]
      bnlocal(2,1) = 'C*'
      bnlocal(2,2) = 'CW'
      bnlocal(2,3) = 'NA'
c     --- CA - CA - NB [wlj1996]
      bnlocal(3,1) = 'CA'
      bnlocal(3,2) = 'CA'
      bnlocal(3,3) = 'NB'
c     --- CB - CA - NA [wlj1996]
      bnlocal(4,1) = 'CB'
      bnlocal(4,2) = 'CA'
      bnlocal(4,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 14
c     --- NA - CR - NB [jpcb1998]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'NB'
c     --- CA - C - OH [wlj19996]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'OH'
c     --- CD - C - OH [wlj19996]
      bnlocal(3,1) = 'CD'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'OH'
c     --- CA - CA - CT [wlj19996]
      bnlocal(4,1) = 'CA'
      bnlocal(4,2) = 'CA'
      bnlocal(4,3) = 'CT'
c     --- N2 - CA - N2 [wlj19996]
      bnlocal(5,1) = 'N2'
      bnlocal(5,2) = 'CA'
      bnlocal(5,3) = 'N2'
c     --- NA - CR - NA [wlj19996]
      bnlocal(6,1) = 'NA'
      bnlocal(6,2) = 'CR'
      bnlocal(6,3) = 'NA'
c     --- NA - CK - NB [wlj19996]
      bnlocal(7,1) = 'NA'
      bnlocal(7,2) = 'CK'
      bnlocal(7,3) = 'NB'
c     --- NA - CR - CT [wlj19996]
      bnlocal(8,1) = 'NA'
      bnlocal(8,2) = 'CR'
      bnlocal(8,3) = 'CT'
c     --- NB - CR - CT [wlj19996]
      bnlocal(9,1) = 'NB'
      bnlocal(9,2) = 'CR'
      bnlocal(9,3) = 'CT'
c     --- CT - CA - NA [wlj19996]
      bnlocal(10,1) = 'CT'
      bnlocal(10,2) = 'CA'
      bnlocal(10,3) = 'NA'
c     --- CA - CA - N [wlj19996]
      bnlocal(11,1) = 'CA'
      bnlocal(11,2) = 'CA'
      bnlocal(11,3) = 'N'
c     --- CA - CB - CT [wlj19996]
      bnlocal(12,1) = 'CA'
      bnlocal(12,2) = 'CB'
      bnlocal(12,3) = 'CT'
c     --- CB - CB - CT [wlj19996]
      bnlocal(13,1) = 'CB'
      bnlocal(13,2) = 'CB'
      bnlocal(13,3) = 'CT'
c     --- CA - CA - OH from GROMACS - AAC
      bnlocal(14,1) = 'CA'
      bnlocal(14,2) = 'CA'
      bnlocal(14,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - CB - NB [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.30d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - CM - HC [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CA - CN - NA [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CN'
      bnlocal(1,3) = 'NA'
c     --- CD - CN - NA [wlj1996]
      bnlocal(2,1) = 'CD'
      bnlocal(2,2) = 'CN'
      bnlocal(2,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.00d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- CA - CT - CT [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
c     --- CB - CT - CT [wlj1996]
      bnlocal(2,1) = 'CB'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'CT'
c     --- CT - CT - CV from GROMACS HIE AAC
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'CV'
c     --- CT - CT - CW from TINKER HIE/HID AAC
      bnlocal(4,1) = 'CT'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.20d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - N2 - CT [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'N2'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - NC - CB [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.6d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CA - NC - CQ [jms1997]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CQ'
c     --- NA - C - NA [wlj1996]
      bnlocal(2,1) = 'NA'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - C - NA [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.80d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - C - O [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CB - C* - CT [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CT'
c     --- CB - CA - CT [wlj1996]
      bnlocal(2,1) = 'CB'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.40d0
      bclocal(1) = (63.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CB - C* - CW [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CW'
c     --- CB - CA - CW [wlj1996]
      bnlocal(2,1) = 'CB'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CA - N2 [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'N2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CA - NC [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1 
c     --- CB - CB - N* [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'N*'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CB - NB [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CB - NC [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - CN - NA [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'CN'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.40d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - N* - CK [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CK'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CB - N* - CT [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CU - CA - CW [jpcb1998]
      bnlocal(1,1) = 'CU'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CW'
c     --- CB - NB - CK [jpcb1998]
      bnlocal(2,1) = 'CB'
      bnlocal(2,2) = 'NB'
      bnlocal(2,3) = 'CK'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- NB - CV - CW [jpcb1998]
      bnlocal(1,1) = 'NB'
      bnlocal(1,2) = 'CV'
      bnlocal(1,3) = 'CW'
c     --- CA - CV - NB [wlj1996]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'CV'
      bnlocal(2,3) = 'NB'
c     --- CB - NC - CQ [wlj1996]
      bnlocal(3,1) = 'CB'
      bnlocal(3,2) = 'NC'
      bnlocal(3,3) = 'CQ'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CA - NC - CA [jms1998]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NC'
      bnlocal(1,3) = 'CA'
c     --- CA - NC - NC [jms1998]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'NC'
      bnlocal(2,3) = 'NC'
c     --- CT - C - O2 [wlj1996]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CK - N* - CT [wlj1996]
      bnlocal(1,1) = 'CK'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CM - C - NA [wlj1996]
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.30d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CM - C - O [wlj1996]
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CA - CA - NT [jacs1999]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NT'
c     --- CM - CA - N2 [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'N2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CM - CA - NC [wlj1996]
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
c     --- CM - CA - NA [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CM - CM - N* [wlj1996]
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'N*'
c     --- CM - CM - NA [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'NA'
c     --- CM - N* - CT [wlj1996]
      bnlocal(3,1) = 'CM'
      bnlocal(3,2) = 'N*'
      bnlocal(3,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CN - NA - CW [wlj1996]
      bnlocal(1,1) = 'CN'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CT - C - N [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
c     --- C - C - N [wlj1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.40d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - C - O [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - C* - CW [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (40.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CA - CT - CA [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.70d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CT - CT - N [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'N'
c     --- CT - CT - NT [wlj1996]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'NT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.70d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - CT - S [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.60d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - CT - SH [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'SH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (50.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1 
c     --- CT - N - CT [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.50d0
      bclocal(1) = (55.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CT - OH - HO [jcc1997]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'HO'
c     --- CO - OH - HO [jcc1997]
      bnlocal(2,1) = 'CO'
      bnlocal(2,2) = 'OH'
      bnlocal(2,3) = 'HO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (60.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CT - OS - CO
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'CO'
c     --- CT - OS - CT
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'OS'
      bnlocal(2,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.50d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CT - OS - P [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'P'
c     --- P - OS - P [wlj1996]
      bnlocal(2,1) = 'P'
      bnlocal(2,2) = 'OS'
      bnlocal(2,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 98.90d0
      bclocal(1) = (62.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - S - CT [WLJ1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 103.70d0
      bclocal(1) = (68.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - S - S [WLJ1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.10d0
      bclocal(1) = (77.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- F - CT - F [wlj1996]
      bnlocal(1,1) = 'F'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'F'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 20
c     --- HC - CT - OS [jcc1997]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'OS'
c     --- HC - CT - OH [jcc1997]
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'OH'
c     --- HC - CO - OS [jcc1997]
      bnlocal(3,1) = 'HC'
      bnlocal(3,2) = 'CO'
      bnlocal(3,3) = 'OS'
c     --- HC - CO - OH [jcc1997]
      bnlocal(4,1) = 'HC'
      bnlocal(4,2) = 'CO'
      bnlocal(4,3) = 'OH'
c     --- CT - NT - H [jacs1999]
      bnlocal(5,1) = 'CT'
      bnlocal(5,2) = 'NT'
      bnlocal(5,3) = 'H'
c     --- HC - CT - NT [jacs1999]
      bnlocal(6,1) = 'HC'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'NT'
c     --- C - CT - HC [wlj1996]
      bnlocal(7,1) = 'C'
      bnlocal(7,2) = 'CT'
      bnlocal(7,3) = 'HC'
c     --- C* - CT - HC [wlj1996]
      bnlocal(8,1) = 'C*'
      bnlocal(8,2) = 'CT'
      bnlocal(8,3) = 'HC'
c     --- CA - CT - HC [wlj1996]
      bnlocal(9,1) = 'CA'
      bnlocal(9,2) = 'CT'
      bnlocal(9,3) = 'HC'
c     --- CM - CT - HC [wlj1996]
      bnlocal(10,1) = 'CM'
      bnlocal(10,2) = 'CT'
      bnlocal(10,3) = 'HC'
c     --- HC - CT - N [wlj1996]
      bnlocal(11,1) = 'HC'
      bnlocal(11,2) = 'CT'
      bnlocal(11,3) = 'N'
c     --- HC - CT - N* [wlj1996]
      bnlocal(12,1) = 'HC'
      bnlocal(12,2) = 'CT'
      bnlocal(12,3) = 'N*'
c     --- HC - CT - NA [wlj1996]
      bnlocal(13,1) = 'HC'
      bnlocal(13,2) = 'CT'
      bnlocal(13,3) = 'NA'
c     --- HC - CT - N2 [wlj1996]
      bnlocal(14,1) = 'HC'
      bnlocal(14,2) = 'CT'
      bnlocal(14,3) = 'N2'
c     --- HC - CT - N3 [wlj1996]
      bnlocal(15,1) = 'HC'
      bnlocal(15,2) = 'CT'
      bnlocal(15,3) = 'N3'
c     --- HC - CT - S [wlj1996]
      bnlocal(16,1) = 'HC'
      bnlocal(16,2) = 'CT'
      bnlocal(16,3) = 'S'
c     --- HC - CT - SH [wlj1996]
      bnlocal(17,1) = 'HC'
      bnlocal(17,2) = 'CT'
      bnlocal(17,3) = 'SH'
c     --- CB - CT - HC [wlj1996]
      bnlocal(18,1) = 'CB'
      bnlocal(18,2) = 'CT'
      bnlocal(18,3) = 'HC'
c     --- CV - CT - HC from GROMACS - AAC
      bnlocal(19,1) = 'CV'
      bnlocal(19,2) = 'CT'
      bnlocal(19,3) = 'HC'
c     --- CW - CT - HC from TINKER for HIE -AAC
      bnlocal(20,1) = 'CW'
      bnlocal(20,2) = 'CT'
      bnlocal(20,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.10d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- HC - CM - N* [wlj1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'N*'
c     --- HC - CM - NA [wlj1996]
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'NA'
c     --- HA - CA - NB [wlj1996]
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CA'
      bnlocal(3,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.05d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- HA - CK - N* [wlj1996]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CK'
      bnlocal(1,3) = 'N*'
c     --- HA - CK - NB [wlj1996]
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CK'
      bnlocal(2,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.45d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NC - CQ - HA [jms1998]
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'CQ'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.50d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HO - OH - P [wlj1996]
      bnlocal(1,1) = 'HO'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 92.07d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HS - SH - HS [wlj1996]
      bnlocal(1,1) = 'HS'
      bnlocal(1,2) = 'SH'
      bnlocal(1,3) = 'HS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.90d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N - C - O [wlj1996]
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.20d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N* - CB - NC [wlj1996]
      bnlocal(1,1) = 'N*'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.90d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N* - CK - NB [wlj1996]
      bnlocal(1,1) = 'N*'
      bnlocal(1,2) = 'CK'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CT - C - CT [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CT'
c     --- NC - CA - CT [wlj1996]
      bnlocal(2,1) = 'NC'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'CT'
c     --- N2 - CA - NA [wlj1996]
      bnlocal(3,1) = 'N2'
      bnlocal(3,2) = 'CA'
      bnlocal(3,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N2 - CA - NC [wlj1996]
      bnlocal(1,1) = 'N2'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.60d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NA - C - O [wlj1996]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- NA - CA - NC [wlj1996]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
c     --- NA - CA - NB [wlj1996]
      bnlocal(2,1) = 'NA'
      bnlocal(2,2) = 'CA'
      bnlocal(2,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.50d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NC - C - O [wlj1996]
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 129.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NC - CQ - NC [jms1998]
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'CQ'
      bnlocal(1,3) = 'NC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.00d0
      bclocal(1) = (80.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- O - C - O [wlj1996]
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
c     --- O2 - C - O2 [wlj1996]
      bnlocal(2,1) = 'O2'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O2'
c     --- O - C - O2 [wlj1996]
      bnlocal(3,1) = 'O'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.90d0
      bclocal(1) = (140.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 - P - O2 [wlj1996]
      bnlocal(1,1) = 'O2'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.23d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 - P - OH [wlj1996]
      bnlocal(1,1) = 'O2'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.23d0
      bclocal(1) = (100.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 - P - OS [wlj1996]
      bnlocal(1,1) = 'O2'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 102.60d0
      bclocal(1) = (45.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- OH - P - OS [wlj1996]
      bnlocal(1,1) = 'OH'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OS'
c     --- OS - P - OS [wlj1996]
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'P'
      bnlocal(2,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.80d0
      bclocal(1) = (33.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HC - CT - HC [jacs1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.70d0
      bclocal(1) = (37.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- HC - CT - CT [jacs1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
C     --- CO - CT - HC [jcc1997]
      bnlocal(2,1) = 'CO'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'HC'
C     --- CT - CO - HC [jcc1997]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'CO'
      bnlocal(3,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.70d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT - CT - CT [jacs1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
C     --- CT - CT - CO [jcc1997]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'CO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.55d0
      bclocal(1) = (92.60d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- OS - CO - OS [jcc1997]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'CO'
      bnlocal(1,3) = 'OS'
C     --- OS - CO - OH [jcc1997]
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'CO'
      bnlocal(2,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CA - CA - NC [jms1998]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
C     --- CM - CM - CT [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.00d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NC - CA - HA [jms1998]
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'HA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.70d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NA - CW - CS [jpcb1998]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.80d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- CW - NA - CW [jpcb1998]
      bnlocal(1,1) = 'CW'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'CW'
C     --- CR - NA - CW [jpcb1998]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'NA'
      bnlocal(2,3) = 'CW'
C     --- CA - NA - CD [wlj1996]
      bnlocal(3,1) = 'CA'
      bnlocal(3,2) = 'NA'
      bnlocal(3,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CW - CS - CS [jpcb1998]
      bnlocal(1,1) = 'CW'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'CS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.60d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CW - NA [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'NA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.70d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CS - CW [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 132.10d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CW - CS [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.50d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CS - CS [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CS'
      bnlocal(1,3) = 'CS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.10d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NA - NB  - CU [jpcb1998]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'NB'
      bnlocal(1,3) = 'CU'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.10d0
      bclocal(1) = (56.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NB - NA - CW [jpcb1998]
      bnlocal(1,1) = 'NB'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NA - CW - CA [jpcb1998]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.40d0
      bclocal(1) = (56.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- H - NA - NB [jpcb1998]
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.90d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CU - NB [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CU'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.90d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CA - CW [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.90d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NB - CU - CA [jpcb1998]
      bnlocal(1,1) = 'NB'
      bnlocal(1,2) = 'CU'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.60d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CU - CA [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CU'
      bnlocal(1,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.70d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- HA - CW - CA [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CA'
C     --- HA - CW - CV [jpcb1998]
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CW'
      bnlocal(2,3) = 'CV'
C     --- HA - CW - CW from TINKER for HIE - AAC
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CW'
      bnlocal(3,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.20d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- HA - CA - CU [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CU'
C     --- HA - CV - CW [jpcb1998]
      bnlocal(2,1) = 'HA'
      bnlocal(2,2) = 'CV'
      bnlocal(2,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- NA - CW - CV [jpcb1998]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CV'
C     --- NA - CW - CW from TINKER for HID/HIE -AAC
      bnlocal(2,1) = 'NA'
      bnlocal(2,2) = 'CW'
      bnlocal(2,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CR - NB - CV [jpcb1998]
      bnlocal(1,1) = 'CR'
      bnlocal(1,2) = 'NB'
      bnlocal(1,3) = 'CV'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.60d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- OS - CW - CS [jpcb1998]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CS'
C     --- CA - CW - OS [jpcb1998]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'CW'
      bnlocal(2,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.50d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CW - OS - CW [jpcb1998]
      bnlocal(1,1) = 'CW'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.40d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CW - OS [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.30d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- OS - NB - CU [jpcb1998]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'NB'
      bnlocal(1,3) = 'CU'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.90d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CW - OS - NB [jpcb1998]
      bnlocal(1,1) = 'CW'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'NB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NB - CR - OS [jpcb1998]
      bnlocal(1,1) = 'NB'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CR - OS - CW [jpcb1998]
      bnlocal(1,1) = 'CR'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CV - CW - OS [jpcb1998]
      bnlocal(1,1) = 'CV'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'OS'
C     --- CT - C - OH [wlj1996]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (35.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- HA - CR - OS [jpcb1998]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'OS'
C     --- CT - CM - HC [wlj1996]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'HC'
C     --- HC - CM - HC [wlj1996]
      bnlocal(3,1) = 'HC'
      bnlocal(3,2) = 'CM'
      bnlocal(3,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.40d0
      bclocal(1) = (43.60d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- H - NT - H [jacs1999]
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NT'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CA - NT - H [jacs1999]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NT'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (56.20d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - CT - NT [jacs1999]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'NT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.20d0
      bclocal(1) = (51.80d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - NT - CT [jacs1999]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'NT'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.80d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- HC - C* - CW [wlj1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CW'
C     --- HC - C* - CB [wlj1996]
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'C*'
      bnlocal(2,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HA - CA - CW [wlj1996]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.00d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 11
C     --- CA - C - CA [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CA'
C     --- CD - C - CD [wlj1996]
      bnlocal(2,1) = 'CD'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'CD'
C     --- C - CA - CA [wlj1996]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'CA'
      bnlocal(3,3) = 'CA'
C     --- CA - CA - CN [wlj1996]
      bnlocal(4,1) = 'CA'
      bnlocal(4,2) = 'CA'
      bnlocal(4,3) = 'CN'
C     --- CD - CA - CD [wlj1996]
      bnlocal(5,1) = 'CD'
      bnlocal(5,2) = 'CA'
      bnlocal(5,3) = 'CD'
C     --- C - CD - CD [wlj1996]
      bnlocal(6,1) = 'C'
      bnlocal(6,2) = 'CD'
      bnlocal(6,3) = 'CD'
C     --- CA - CD - CD [wlj1996]
      bnlocal(7,1) = 'CA'
      bnlocal(7,2) = 'CD'
      bnlocal(7,3) = 'CD'
C     --- CB - CD - CD [wlj1996]
      bnlocal(8,1) = 'CB'
      bnlocal(8,2) = 'CD'
      bnlocal(8,3) = 'CD'
C     --- CD - CD - CD [wlj1996]
      bnlocal(9,1) = 'CD'
      bnlocal(9,2) = 'CD'
      bnlocal(9,3) = 'CD'
C     --- CD - CD - CN [wlj1996]
      bnlocal(10,1) = 'CD'
      bnlocal(10,2) = 'CD'
      bnlocal(10,3) = 'CN'
C     --- CA - CA - NO (wlj nitro) [wlj2001]
      bnlocal(11,1) = 'CA'
      bnlocal(11,2) = 'CA'
      bnlocal(11,3) = 'NO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HC - C - N [wlj1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.80d0
      bclocal(1) = (69.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - CT - Cl [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'Cl'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.70d0
      bclocal(1) = (78.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- Cl - CT - Cl [wlj1996]
      bnlocal(1,1) = 'Cl'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'Cl'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.60d0
      bclocal(1) = (51.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HC - CT - Cl [wlj1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'Cl'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.20d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- N - C - N [wlj1996]
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.00d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- O - C - OH [wlj1996]
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.40d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CB - C* - CW [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'C*'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.40d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CA - CA - CW [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NC - CQ - CT [wlj1996]
      bnlocal(1,1) = 'NC'
      bnlocal(1,2) = 'CQ'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.50d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- NB - CV - CT [wlj1996]
      bnlocal(1,1) = 'NB'
      bnlocal(1,2) = 'CV'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.20d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C - CB - CB [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 134.90d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
C     --- C* - CB - CA [wlj1996]
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CA'
C     --- C* - CB - CD [wlj1996]
      bnlocal(2,1) = 'C*'
      bnlocal(2,2) = 'CB'
      bnlocal(2,3) = 'CD'
C     --- CA - CB - CA [wlj1996]
      bnlocal(3,1) = 'CA'
      bnlocal(3,2) = 'CB'
      bnlocal(3,3) = 'CA'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.80d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C* - CB - CN [wlj1996]
      bnlocal(1,1) = 'C*'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CN'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.30d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CA - CB - CB [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CB'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.20d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CA - CB - CN [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CB'
      bnlocal(1,3) = 'CN'
C     --- CD - CB - CN [wlj1996]
      bnlocal(2,1) = 'CD'
      bnlocal(2,2) = 'CB'
      bnlocal(2,3) = 'CN'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.70d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C - CM - CM [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.00d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CA - CM - CM [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'CM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - C - HC [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.70d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CA - CN - CB [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CN'
      bnlocal(1,3) = 'CB'
C     --- CB - CN - CD [wlj1996]
      bnlocal(2,1) = 'CB'
      bnlocal(2,2) = 'CN'
      bnlocal(2,3) = 'CD'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (41.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HC - CT - P [wlj1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (43.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT - CT - P [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'P'
C     --- CA - CT - P [wlj1996]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'P'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.80d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- C - N - H [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'H'
C     --- CA - N - H [wlj1996]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'N'
      bnlocal(2,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.40d0
      bclocal(1) = (38.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - N - H [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.40d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CT - N2 - H from GROMACS - AAC
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N2'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.20d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C - N* - H [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 127.30d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CB - N* - H [wlj1996]
      bnlocal(1,1) = 'CB'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.20d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CM - N* - H [wlj1996]
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'N*'
      bnlocal(1,3) = 'H'
C     --- CM - NA - H [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'NA'
      bnlocal(2,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - N3 - CT [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'N3'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.80d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C - NA - H [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.00d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CA - NA - H [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.10d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CN - NA - H [wlj1996]
      bnlocal(1,1) = 'CN'
      bnlocal(1,2) = 'NA'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 123.40d0
      bclocal(1) = (83.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- O - C - OS [wlj1996]
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.00d0
      bclocal(1) = (100.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- OS - CT - OS [wlj1996]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'OS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.50d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT - P - OS [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'P'
      bnlocal(1,3) = 'OS'
C     --- CT - P - O2 [wlj1996]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'P'
      bnlocal(2,3) = 'O2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 96.00d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - SH - HS [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'SH'
      bnlocal(1,3) = 'HS'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.20d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CT - C - C [wlj1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.40d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C - C - O [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.0d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- ON - NO - ON (wlj nitro) [wlj2001]
      bnlocal(1,1) = 'ON'
      bnlocal(1,2) = 'NO'
      bnlocal(1,3) = 'ON'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.5d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CT - NO -ON (wlj nitro) [wlj2001]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'NO'
      bnlocal(1,3) = 'ON'
C     --- CA - NO -ON (wlj nitro) [wlj2001]
      bnlocal(2,1) = 'CA'
      bnlocal(2,2) = 'NO'
      bnlocal(2,3) = 'ON'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.0d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- HC - CT - NO (wlj nitro) [wlj2001]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'NO'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.7d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 0
c     --- CT - CV - CW for HIE from GROMACS - AAC
      nbslocal = nbslocal + 1
      bnlocal(nbslocal,1) = 'CT'
      bnlocal(nbslocal,2) = 'CV'
      bnlocal(nbslocal,3) = 'CW'
c     --- CT - CW - CW for HIE from TINKER - AAC
      nbslocal = nbslocal + 1
      bnlocal(nbslocal,1) = 'CT'
      bnlocal(nbslocal,2) = 'CW'
      bnlocal(nbslocal,3) = 'CW'
c     --- CT - CW - CV for HID inspired by TINKER value
      nbslocal = nbslocal + 1
      bnlocal(nbslocal,1) = 'CT'
      bnlocal(nbslocal,2) = 'CW'
      bnlocal(nbslocal,3) = 'CV'
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
c     --- Unlike the use in the ffbond and ffangle subroutines, OPLS 
c     uses a three letter code for determining torsions as some times 
c     the same atom types have different torsional interactions.

c     H-C-C-H alkanes
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.318d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HC - CT - CT - HC [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
c     --- HC - CM - CT - HC [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-C alkanes
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.366d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC - CT - CT - CT [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-C alkanes
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.740d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.157d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.279d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - CT - CT [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C=C alkenes
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.372d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC - CT - CM - CM [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'CM'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- null Torsion with all parameters set to 0.0
c     --- Still needed so that the 1-4 interactions are computed
c     --- properly
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 0
c     H-C-C[ar]-C[ar] ethylbenzene
c     --- HC - CT - CA - CA [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CA~'
      tnlocal(ntslocal,4) = 'CA~'
c     C-C-C[ar]-C[ar] ethylbenzene
c     --- CT - CT - CA - CA [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CA~'
      tnlocal(ntslocal,4) = 'CA~'
c     H-C-C(O)-O carboxylate ion
c     --- HC~ - CT - C_a - O2 [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C a'
      tnlocal(ntslocal,4) = 'O2'
c     C-C-C(O)-H aldehyde
c     --- CT - CT - C_k - HC~ [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C k'
      tnlocal(ntslocal,4) = 'HC~'
c     H-C-C(O)-N acetamide and peptide Psi''
c     --- HC~ - CT - C_n - N [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C n'
      tnlocal(ntslocal,4) = 'N'
c     H-C-N-C tertiary amide
c     --- HC~ - CT - N - CT [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'N'
      tnlocal(ntslocal,4) = 'CT'
c     C(O)-N-C-H  peptide Phi''
c     --- C_n - N - CT - HC~ [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C n'
      tnlocal(ntslocal,2) = 'N'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'HC~'
c     H-N-C-X peptides
c     --- H - N - CT - X [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'H'
      tnlocal(ntslocal,2) = 'N'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'wild'
c     H-C-C3-C 3-ethylindole
c     --- HC~ - CT - CA~ - CAj [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CA~'
      tnlocal(ntslocal,4) = 'CAj'
c     C-C-C3-C 3-ethylindole 
c     --- CT - CT - CA~ - CAj [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CA~'
      tnlocal(ntslocal,4) = 'CAj'
c     H-C-C(O)-O aldehyde 
c     -- HC~ - CT - C_k - O [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C k'
      tnlocal(ntslocal,4) = 'O'
c     H-C-C(O)-O amide
c     -- HC~ - CT - C_n - O [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C n'
      tnlocal(ntslocal,4) = 'O'
c     HC-CT-NO-ON nitro compounds
c     -- HC~ - CT - NO - ON [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'NO'
      tnlocal(ntslocal,4) = 'ON'
c     --NT+  - CT -- Cn  - O for N-term - from Gromacs -AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'NT+'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C n'
      tnlocal(ntslocal,4) = 'O'
c     -- N - CT - C_n - O for peptide bond from Gromacs - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C n'
      tnlocal(ntslocal,4) = 'O'
c     --- HC~ - CT - N2 - H for guanidinium & ARG-Gromacs&Tinker-AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'N2'
      tnlocal(ntslocal,4) = 'H'
c     --- CT - CT - CV - CW for HIS from GROMACS&Tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CV'
      tnlocal(ntslocal,4) = 'CW'
c     --- HC~ - CT - CV - CW for HIS from GROMACS&Tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CV'
      tnlocal(ntslocal,4) = 'CW'
c     --- N - CT - C_a - O2 for C-term from GROMACS&Tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C a'
      tnlocal(ntslocal,4) = 'O2'
c     --- CT - CT - CW - CW for HIE/HID from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'CW'
c     --- CW - CW - CT - HC~ for HIE/HID from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CW'
      tnlocal(ntslocal,2) = 'CW'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'HC~'
c     --- N - C n - N - H for UREA from gromacs -AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'C n'
      tnlocal(ntslocal,3) = 'N'
      tnlocal(ntslocal,4) = 'H'
c     --- CT-CT-CW-?? for TRP, taken from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'wild'
c     --- HC~-CT-CW-?? taken from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'wild'
c     --- CT-CT-C*-CB taken from [wlj2001] for TRP
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C*'
      tnlocal(ntslocal,4) = 'CB'
c     --- HC-CT-C*-?? taken from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C*'
      tnlocal(ntslocal,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
      
c     H-C-C-C[ar] ethylbenzene
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.462d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 0
c     --- HC - CT - CT - CA [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CA~'
c     --- HC - CT - CT - CV aromatic from gromacs&tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CV'
c     --- HC - CT - CT - CW for HID/HIE from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CW'
c     --- HC-CT-CT-C* from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'HC~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C*'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-O-H alcohols
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.450d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HC - CT - OH - HO (mono)[jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
c     --- HC - CT - OH - HO (poly)[jcc1997]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'OHp'
      tnlocal(2,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-O-H alcohols
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.356d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.174d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.492d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - OH - HO [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-O alcohols
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.468d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- HC - CT - CT - OH (mono) [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
c     --- HC - CT - CT - OH (poly) [jcc1997]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OHp'
c     --- HC - CT - CT - OS [jcc1997]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-O alcohols
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.711d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.663d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - CT - OH [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H - O - C[bz] - C[bz] phenol
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.682d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HO - OH - CA - CA [jacs1996]
      tnlocal(1,1) = 'HO'
      tnlocal(1,2) = 'OHp'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-S-H thiol
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.451d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC - CT - SH - HS [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-S-H thiol
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.759d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.282d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.603d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - SH - HS [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-S thiol
c     H-C-C-S sulfide (added by M.G. Martin 07-05-2001)
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.452d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HC - CT - CT - SH [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     --- HC - CT - CT - S [mgm]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-S thiol
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.876d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - CT - SH [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-S-C sulfide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.647d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC - CT - S - CT [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-S sulfide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.619d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.620d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.258d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - CT - S [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-S-C sulfide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.925d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.576d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.677d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - S -CT [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-S-S-C disulfide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-7.414d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.705d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - S - S -CT [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-S-S disulfide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.558d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC - CT - S - S [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-S-S disulfide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.941d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.836d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.935d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - S - S [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.400d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- HC - CT - N3 - H (primary amine) [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N3~'
      tnlocal(1,4) = 'H'
c     --- HC - CT - NT - H (primary amine) [jacs1999]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'NT~'
      tnlocal(2,4) = 'H'
c     --- HC - CT - CT - Cl (alkyl chloride) [wlj1996]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'Cl'
c     --- CT - CT - CT - Cl (alkyl chloride) [wlj1996]
      tnlocal(4,1) = 'CT'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'Cl'
c     --- HC - CT - NT+ - H for N terminus - AAC
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'NT+'
      tnlocal(5,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-N primary amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.013d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.709d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.473d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- HC - CT - CT - N3 [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'N3~'
c     --- HC - CT - CT - NT [jacs1999]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'NT~'
c     --- HC - CT - CT - NT+ for N term from GROMACS - AAC
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'NT+'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-N-H primary amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.190d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.417d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.418d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- CT - CT - N3 - H [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N3~'
      tnlocal(1,4) = 'H'
c     --- CT - CT - NT - H [jacs1999]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'NT~'
      tnlocal(2,4) = 'H'
c     --- C_n - CT - NT+ - H for N term from GROMACS - AAC
      tnlocal(3,1) = 'C n'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'NT+'
      tnlocal(3,4) = 'H'
c     --- CT - CT - NT+ - H from GROMACS for N term - AAC
      tnlocal(4,1) = 'CT'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'NT+'
      tnlocal(4,4) = 'H'
c     --- CT - CT - N2 - H for ARG from GROMACS(not Tinker) - AAC
      tnlocal(5,1) = 'CT'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'N2'
      tnlocal(5,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-N primary amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.392d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.674d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.550d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CT - CT - CT - N3 [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'N3~'
c     --- CT - CT - CT - NT [jacs1999]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'NT~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-N-H ammonium ion
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.261d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - N3+ - H [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N3+'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-N-H ammonium ion
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.347d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - N3+ - H [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N3+'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-N ammonium ion
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.384d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - CT - N3+ [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'N3+'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-N ammonium ion
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.732d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.229d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.485d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- CT - CT - CT - N3+ [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'N3+'
c     --- CT - CT - CT - N2 for Arg - TINKER/GROMACS - AAC
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'N2'
c     --- CT - CT - CT - NT+ for Lys - Tinker - AAC
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'NT+'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-O-C ethers
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.760d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - OS - CT [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-O-C-C ethers
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.650d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.250d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.670d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - OS - CT - CT [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-O-C-O acetals
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
c     these parameters are from the original jacs1996 paper, they were
c     changed in the jcc1997 paper
c      tclocal(1) =  (-0.574d0)*(CON_KCALTOK)*(0.5d0)
c      tclocal(2) =  (-0.997d0)*(CON_KCALTOK)*(0.5d0)
c      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.375d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.358d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.004d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - OS - CT - OS [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C(O)-O carboxylate ion
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.820d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C_a - O2 [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'O2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.225d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     H-C-C-C(O)     carboxylate ion
c     --- HC~ - CT - CT - C_a [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C a'
c     HC-CT-CT-NO nitroethane
c     --- HC~ - CT - CT - NO [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'NO'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-C(O) carboxylate ion
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-3.185d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.825d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.493d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - CT - C_a [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C a'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     O-C-O-H carboxylic acids
c     C-C-O-H carboxylic acids
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.83d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- O2 - C_a - OHa - HO [jacs1996]
      tnlocal(1,1) = 'O2'
      tnlocal(1,2) = 'C a'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
c     --- CT - C_a - OHa - HO [jacs1996]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'C a'
      tnlocal(2,3) = 'OHa'
      tnlocal(2,4) = 'HO'
c     --- CT - CV - NB - CR for HIS from TINKER (uses 4.80?) -AAC
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CV'
      tnlocal(3,3) = 'NB'
      tnlocal(3,4) = 'CR'
c     --- CW - CV - NB - CR for HIS from TINKER (uses 4.80?) -AAC
      tnlocal(4,1) = 'CW'
      tnlocal(4,2) = 'CV'
      tnlocal(4,3) = 'NB'
      tnlocal(4,4) = 'CR'
c     --- HA - CV - NB - CR for HIS from TINKER (uses 4.80?) -AAC
      tnlocal(5,1) = 'HA'
      tnlocal(5,2) = 'CV'
      tnlocal(5,3) = 'NB'
      tnlocal(5,4) = 'CR'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     HC-CT-C-HC aldehyde
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.360d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - C_k - HC~ [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     HC-CT-C(O)-CT ketone
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.275d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - C_k - CT [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     CT-CT-C-O      aldehyde and ketone
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.277d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.228d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.694d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C_k - O [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     CT-CT-C(O)-CT ketone
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.454d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.144d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.775d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C_k - CT [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     HC-CT-CT-C(O) aldehyde and ketone
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.076d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - CT - C_k [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     CT-CT-CT-C(O) aldehyde and ketone 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.697d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.456d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.585d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 0
c     --- CT - CT - CT - C_k [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C k'
c     from GROMACS,Tinker - AAC
c     --- C_n - CT- CT - C_n
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C n'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C n'
c     --- N - CT - CT - CV for aromatic chi1 - Gromacs&tinker -AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CV'
c     --- C_n - CT - CT - CV for HIS from TINKER (!Gromacs) -AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C n'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CV'
c     --- CA~ - CT - CT - C_n for aromatic from TINKER -AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CA~'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C n'
c     --- C_n - CT - CT - C_a for aromatic from TINKER -AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C n'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C a'
c     --- C_n - CT - CT - CW for HID/HIE from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C n'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CW'
c     --- C*-CT-CT-C from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C*'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C n'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C(O)-N-C-H amide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.139d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- C_n - N - CT - HC~ [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C(O)-N-C-C alkyl-amide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.396d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.427d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- C_n - N - CT - CT [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C-H N-ethylformamide, peptide side chains
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.464d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- N - CT - CT - HC~ [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C-C N-propylformamide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.964d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.659d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- N - CT - CT - CT [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C(O)-N propanamide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (3.250d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.402d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.136d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C_n - N [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C(O)-O propanamide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.166d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C_n - O [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-C(O)  propanamide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.100d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - CT - C_n [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C n'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-C(O)    butanamide
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-2.060d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.313d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.315d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - CT - C_n [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C n'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H,C-C(O)-N-H amides
c     O-C(O)-N-H amides 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- HC~ - C_n - N - H [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'H'
c     --- CT - C_n - N - H [jacs1996]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'H'
c     --- O - C_n - N - H [jacs1996]
      tnlocal(3,1) = 'O'
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H,C-C(O)-N-C amides
c     O-C(O)-N-C amides 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.800d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.089d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- HC~ - C_n - N - CT [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
c     --- CT - C_n - N - CT [jacs1996]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'CT'
c     --- O - C_n - N - CT [jacs1996]
      tnlocal(3,1) = 'O'
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C(O)-N-C-C(O) peptides Phi
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-2.365d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.912d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.850d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- C_n - N - CT - C_k  [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C k'
c     --- C_n - N - CT - C_a  [jacs1996]
      tnlocal(2,1) = 'C n'
      tnlocal(2,2) = 'N'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'C a'
c     --- C_n - N - CT - C_n for peptide Phi - AAC
      tnlocal(3,1) = 'C n'
      tnlocal(3,2) = 'N'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'C n'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C(O)-N-C-C peptides Phi'
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.462d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- C_n - N - CT - CT  [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C(O)-N peptides Psi
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.816d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.222d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.581d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- N - CT - C_n - N  [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'N'
c     --- NT+ - CT - C_n - N for N-terminus AAC
      tnlocal(2,1) = 'NT+'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'C n'
      tnlocal(2,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     Cbeta-C-C(O)-N peptides Psi'
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.173d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.189d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-1.200d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - C_n - N  [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C-C peptides Chi-1 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.845d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.962d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.713d0)*(CON_KCALTOK)*(0.5d0)
c      tclocal(3) =  (-1.200d0)*(CON_KCALTOK)*(0.5d0)
c     --- commented out by AAC - duplicate entry?
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 0
c     --- N - CT - CT - CT  [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CT'
c     --- N - CT - CT - C_n from gromacs&tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C n'
c     --- N - CT - CT - CA~ aromatic from tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CA~'
c     --- N - CT - CT - C_a from gromacs&tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C a'
c     --- N - CT - CT - CV aromatic from tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CV'
c     --- N - CT - CT - CW aromatic from tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'CW'
c     --- C_a - CT - CT - C_a from gromacs&tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C a'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C a'
c     --- C_a - CT - CT - NT+ from gromacs&tinker - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'C a'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'NT+'
c     --- N -CT-CT-C* from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'N'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CT'
      tnlocal(ntslocal,4) = 'C*'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C-O Chi for Ser and Thr
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (6.280d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.467d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (2.030d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- N - CT - CT - OH  [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C(O)-C-C-O Chi for Ser and Thr
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-6.180d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- C_n - CT - CT - OH  [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
c     --- C_a - CT - CT - OH  [jacs1996]
      tnlocal(2,1) = 'C a'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C-S Chi for Cys
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.583d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.163d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.141d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- N - CT - CT - SH  [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     --- N - CT - CT - S  [jacs1996]
      tnlocal(2,1) = 'N'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C(O)-C-C-S Chi for Cys
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-4.214d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-2.114d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.969d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- C_n - CT - CT - SH  [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     --- C_n - CT - CT - S  [jacs1996]
      tnlocal(2,1) = 'C n'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'S'
c     --- C_a - CT - CT - SH  [jacs1996]
      tnlocal(3,1) = 'C a'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'SH'
c     --- C_a - CT - CT - S  [jacs1996]
      tnlocal(4,1) = 'C a'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C-N   HID, HIE, HIP
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.420d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HC~ - CT - CA - NA  [wlj1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'NA'
c     --- HC~ - CT - CA - NB  [wlj1996]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CA~'
      tnlocal(2,4) = 'NB'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-N   HID, HIE, HIP 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.349d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.312d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.540d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CT - CT - CA - NA  [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'NA'
c     --- CT - CT - CA - NB  [jacs1996]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CA~'
      tnlocal(2,4) = 'NB'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     H-C-C3-C2 3-ethylindole
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.480d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC~ - CT - CA - CA  [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C3-C2 3-ethylindole
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.714d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 0
c     --- CT - CT - CA - CA  [jacs1996]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'CA~'
      tnlocal(ntslocal,4) = 'CA~'
c     --- CT-CT-C*-CW for TRP from [wlj2001]
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CT'
      tnlocal(ntslocal,3) = 'C*'
      tnlocal(ntslocal,4) = 'CW'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     O-C-O-H carbohydrates
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.257d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.806d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.003d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- OS - CT - OHp - HO  [jcc1997]
      tnlocal(1,1) = 'OS'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHp'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-C-O carbohydrates, ethers
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.336d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CT - CT - CT - OHp  [jcc1997]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHp'
c     --- CT - CT - CT - OS  [jcc1997]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     C-C-O-H polyols, carbohydrates
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.674d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-2.883d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.026d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - CT - OHp - HO  [jcc1997]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHp'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     OH-C-C-OS carbohydrates 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.319d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- OHp - CT - CT - OS  [jcc1997]
      tnlocal(1,1) = 'OHp'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     OH-C-C-OH carbohydrates, polyols
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (9.066d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- OHp - CT - CT - OHp  [jcc1997]
      tnlocal(1,1) = 'OHp'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHp'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     CT-NT-CT-CT aliphatic amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.416d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.128d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.695d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CT - NT - CT - CT  [jacs1999]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'NT~'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     HC-CT-NT-CT  aliphatic amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.560d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- HC - CT - NT - CT  [jacs1999]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     CA-CA-NT-H  aromatic amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.030d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CA~ - CA~ - NT~ - H  [jacs1999]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     CA-CA-NT-CT  aromatic amine
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-7.582d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.431d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (3.198d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CA~ - CA~ - NT~ - CT  [jacs1999]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (7.250d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- X - CA~ - CA~ - Y  [jacs1999]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'wild'
c     --- ?? - NC - CA - ??  [wlj2001]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'NC'
      tnlocal(2,3) = 'CA~'
      tnlocal(2,4) = 'wild'
c     --- ?? - CS - CW - ??  [wlj2001]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CS'
      tnlocal(3,3) = 'CW'
      tnlocal(3,4) = 'wild'
c     --- CT - CV - CW - HA from Tinker&Gromacs aromatic - AAC
      tnlocal(4,1) = 'CT'
      tnlocal(4,2) = 'CV'
      tnlocal(4,3) = 'CW'
      tnlocal(4,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.25d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- Cl - CT - CT - Cl (dichloride) [wlj1996]
      tnlocal(1,1) = 'Cl'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'Cl'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (8.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- NT - CT - CT - OH (2-aminoethanol 6-31G* fit wj) [wlj2001]
      tnlocal(1,1) = 'NT~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.550d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- OS - CT - CT - OS (polyethers, crown ethers) [wlj2001]
      tnlocal(1,1) = 'OS'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (14.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- ?? - CM - CM - ?? (alkene) [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'wild'
c     --- CT - CM - CM - CT (alkene) [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'CM'
      tnlocal(2,4) = 'CT'
c     --- CT - CM - CM - HC (alkene) [wlj2001]
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CM'
      tnlocal(3,3) = 'CM'
      tnlocal(3,4) = 'HC~'
c     --- HC - CM - CM - HC (alkene) [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CM'
      tnlocal(4,3) = 'CM'
      tnlocal(4,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.40d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-NO-ON nitroethane
c     --- CT - CT - NO - ON [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'NO'
      tnlocal(1,4) = 'ON'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.15d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CA-CA-NO-ON nitrobenzene
c     --- CA~ - CA~ - NO - ON [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'NO'
      tnlocal(1,4) = 'ON'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.20d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     --- ?? - NA - CW - ?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'wild'
c     --- C? - NA - CW - C? [wlj2001]
      tnlocal(2,1) = 'Cwild'
      tnlocal(2,2) = 'NA'
      tnlocal(2,3) = 'CW'
      tnlocal(2,4) = 'Cwild'
c     --- N? - NA - CW - C? [wlj2001]
      tnlocal(3,1) = 'Nwild'
      tnlocal(3,2) = 'NA'
      tnlocal(3,3) = 'CW'
      tnlocal(3,4) = 'Cwild'
c     --- C? - NA - CW - H? [wlj2001]
      tnlocal(4,1) = 'Cwild'
      tnlocal(4,2) = 'NA'
      tnlocal(4,3) = 'CW'
      tnlocal(4,4) = 'Hwild'
c     --- N? - NA - CW - H? [wlj2001]
      tnlocal(5,1) = 'Nwild'
      tnlocal(5,2) = 'NA'
      tnlocal(5,3) = 'CW'
      tnlocal(5,4) = 'Hwild'
c     --- H - NA - CW - wild [wlj2001]
      tnlocal(6,1) = 'H'
      tnlocal(6,2) = 'NA'
      tnlocal(6,3) = 'CW'
      tnlocal(6,4) = 'wild'
c     --- H - NA - CR - wild [wlj2001]
      tnlocal(7,1) = 'H'
      tnlocal(7,2) = 'NA'
      tnlocal(7,3) = 'CR'
      tnlocal(7,4) = 'wild'
c     --- H - NA - CB - wild [wlj2001]
      tnlocal(8,1) = 'H'
      tnlocal(8,2) = 'NA'
      tnlocal(8,3) = 'CB'
      tnlocal(8,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.753d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.734d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Pro CD-N-CA-CB from TINKER&GROMACS - AAC
c     --- CT - N - CT - CT
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.737d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.251d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-3.501d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Pro CD-N-CA-C from TINKER&GROMACS - AAC
c     --- C_n - CT - N - CT
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.829d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.243d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.498d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Guanidinium (ARG?) from TINKER&GROMACS - AAC
c     --- CT - CT - N2 - CA+
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N2'
      tnlocal(1,4) = 'CA+'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (7.936d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Guanidinium (ARG) from TINKER&GROMACS - AAC
c     --- CT - N2 - CA+ - N2
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'N2'
      tnlocal(1,3) = 'CA+'
      tnlocal(1,4) = 'N2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.582d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Guanidinium (ARG) from GROMACS,wlj1996,not_tinker - AAC
c     --- HC~ - CT - CT - N2
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'N2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.90d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Guanidinium (ARG) from GROMACS,wlj1996,Tinker - AAC
c     --- H - N2 - CA+ - N2
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'N2'
      tnlocal(1,3) = 'CA+'
      tnlocal(1,4) = 'N2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.177d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Guanidinium (ARG) from GROMACS,Tinker - AAC
c     --- HC~ - CT- N2 - CA+
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N2'
      tnlocal(1,4) = 'CA+'
      call twh_set_tornames(type,ntslocal,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.366d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.262d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.505d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CT - CT- CV - NB for HID/HIE from GROMACS&TINKER - AAC
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CV'
      tnlocal(1,4) = 'NB'
c     --- CT - CT - CW - NA for HIE/HID from TINKER - AAC
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CW'
      tnlocal(2,4) = 'NA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (10.75d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 0
c     --- CT - CV- CW - NA for HID/HIE from GROMACS&TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CV'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'NA'
c     --- NB - CV - CW - NA for aromatic from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'NB'
      tnlocal(ntslocal,2) = 'CV'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'NA'
c     --- NB - CV - CW - HA for aromatic from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'NB'
      tnlocal(ntslocal,2) = 'CV'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'HA'
c     --- CT - CW - CW - NA for HID/HIE from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CW'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'NA'
c     --- CT - CW - CW - HA for HID/HIE from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'CT'
      tnlocal(ntslocal,2) = 'CW'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'HA'
c     --- NA - CW - CW - NA for HIE/HID from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'NA'
      tnlocal(ntslocal,2) = 'CW'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'NA'
c     --- NA - CW - CW - HA for HIE/HID from TINKER - AAC
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'NA'
      tnlocal(ntslocal,2) = 'CW'
      tnlocal(ntslocal,3) = 'CW'
      tnlocal(ntslocal,4) = 'HA'
c     --- ??-CW-CV-?? taken from [wlj2001] for HIS
      ntslocal = ntslocal + 1
      tnlocal(ntslocal,1) = 'wild'
      tnlocal(ntslocal,2) = 'CW'
      tnlocal(ntslocal,3) = 'CV'
      tnlocal(ntslocal,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (10.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- NA - CR - NB - CV for aromatic from GROMACS&TINKER - AAC
      tnlocal(1,1) = 'NA'
      tnlocal(1,2) = 'CR'
      tnlocal(1,3) = 'NB'
      tnlocal(1,4) = 'CV'
c     --- HA - CR - NB - CV for aromatic from TINKER - AAC
      tnlocal(2,1) = 'HA'
      tnlocal(2,2) = 'CR'
      tnlocal(2,3) = 'NB'
      tnlocal(2,4) = 'CV'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.419d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- HC~ - CT- CW - NA for aromatic from TINKER&GROMACS - AAC
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'NA'
c     --- HC~ - CT- CV - NB for aromatic from TINKER&GROMACS - AAC
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CV'
      tnlocal(2,4) = 'NB'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.65d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- CW - NA - CR - NB for aromatic from GROMACS&TINKER - AAC
      tnlocal(1,1) = 'CW'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CR'
      tnlocal(1,4) = 'NB'
c     --- CW - NA - CR - NA for aromatic from TINKER - AAC
      tnlocal(2,1) = 'CW'
      tnlocal(2,2) = 'NA'
      tnlocal(2,3) = 'CR'
      tnlocal(2,4) = 'NA'
c     --- CW - NA - CR - HA for aromatic from TINKER - AAC
      tnlocal(3,1) = 'CW'
      tnlocal(3,2) = 'NA'
      tnlocal(3,3) = 'CR'
      tnlocal(3,4) = 'HA'
      call twh_set_tornames(type,ntslocal,tnlocal)
    
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (13.05d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- ??-C*-CW-?? from [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C*'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
    
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.35d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- ??-C*-CB-?? from [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C*'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
    
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- ??-CB-CN-?? from [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'CN'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
    
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (7.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- ??-CA-CB-?? from [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
    
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.05d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- ??-NA-CN-?? from [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CN'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
    
      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * OPLS-aa improper torsions                                      *
c     ******************************************************************
c     According to personal communication with Julian Tirado-Rives I 
c     now think they actually are using the sterocenter version - the 
c     same that Amber uses.  This would make sense as they share
c     many bonded parameters with the Amber family of force fields.
c     The listing of i-j-k-l has j in the center, bonded to i,k, and l

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,2)
      iclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(2) =  (2.200d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 5
c     --- Z - CA~ - X - Y [wlj2001]
c     ---any aromatic carbon i.e. CW, CA+, CA~, CAj, CB', CR, C*, CV
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'CA~'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'wild'
c     --- Z - CW - X - Y
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'CW'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'wild'
c     --- Z - CR - X - Y
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'CR'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'wild'
c     --- Z - CV - X - Y
      ninlocal(4,1) = 'wild'
      ninlocal(4,2) = 'CV'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'wild'
c     --- Z - C* - X - Y
      ninlocal(5,1) = 'wild'
      ninlocal(5,2) = 'C*'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'wild'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,2)
      iclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(2) =  (21.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 4
c     --- Z - C - X - Y [wlj2001]
c     --- generically (O/Z)-C-X-Y for non aromatic, trigonal planar carbons
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'C a'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'wild'
c     --- Z - C - X - Y [wlj2001]
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'C k'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'wild'
c     --- Z - C - X - Y [wlj2001]
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'C n'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'wild'
c     --- Z - CA+ - X - Y for arginine CZ AAC 9/05
      ninlocal(4,1) = 'wild'
      ninlocal(4,2) = 'CA+'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'wild'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,2)
      iclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(2) =  (2.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 3
c     --- Z - N - X - Y [wlj2001]
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'N'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'wild'
c     --- Z - N2 - X - Y  for ARG - AAC
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'N2'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'wild'
c     --- Z - NA - X - Y  for HIS - AAC
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'NA'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'wild'
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

      return
      end
