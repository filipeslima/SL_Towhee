#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2003-2014 Marcus G. Martin                       *
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
      subroutine twh_ffopls2001(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the opls-aa force*
c     * field.                                                         *
c     * This forcefield is based upon the 2001 reparameterization of   *
c     * OPLS as documented in                                          *
c     * Kaminski et al J. Phys. Chem. B 105, 6474-6487 (2001)          *
c     * suplemented with some files I recieved from W.L. Jorgensen     *
c     * sent to me in 2001 and also 2004                               *
c     * Note that the OPLS force field only uses the first 2 character *
c     * for the bond and angle assignment and 3 characters for torsion *
c     * this is handled using the c2temp and c3temp variables          *
c     *                                                                *
c     * rewritten from oplsaa.F 09-23-2003 by M.G. Martin              *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      integer ii,type,isame
      character*1 c1temp
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

      ffname = 'OPLS-2001'
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
      nbnlocal ='CAj'
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
c     (CM) alkene sp2 carbon (non-aromatic) [jacs1996]
      nbnlocal ='CM'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

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
c     (C=) interior sp2 carbons of a diene (ex. C=C-C=C center atoms)
c     [wlj2001]
      nbnlocal ='C='
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
c     (CY) aliphatic sp3 carbon in a small ring (cyclopropane) 
c     [wlj2001]
      nbnlocal ='CY'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (CTf) aliphatic sp3 hybrid carbon bonded to F  
c     in monoalkyl fluorides and perfluoroalkanes (excludes CF4) 
c     [wlj2001]
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
      nbnlocal ='CTf'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (CTcf4) aliphatic sp3 hybrid carbon bonded to F in CF4 
c     [wlj2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.097d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CTcf4'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      ctnlocal = 'none'
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
      
c     (F) Fluorine nonionic [wlj2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.94
      nbclocal(2) = (0.061d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
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

c     (F) Fluorine nonionic in perfluoroalkanes[wlj2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.95
      nbclocal(2) = (0.053d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fpf'
      c1temp = 'F'
      cbnlocal = c1temp
      canlocal = c1temp
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
c     (N) N in an amide [wlj1996]
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
c     (NA) sp2 aromatic N with H attached [jpcb1998]
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
c     (NB) sp2 N in 5-mem ring with lone pair [jpcb1998]
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
c     (NC) sp2 N in 6-mem ring with lone pair [jms1998] 
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
c     (N2) sp2 N of aromatic amines and guanidinium ions [jacs1991]
      nbnlocal ='N'
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
c     (NT+1) sp3 N in ammonium ion.[jacs1999]
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
c     (NT+2) sp3 N in primary,secondary,tertiary ammonium ion [jacs1999]
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
c     (O2) Oxygen in COO- or COOH [jacs1996]
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
c     (O) Oxygen in C=O [jacs1996]
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
c     (OHp, OHd, OHt) O in phenol,diols,triols or polyols [wlj2001] 
      nbnlocal ='OHp'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

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
c     (OHp, OHd, OHt) O in phenol,diols,triols or polyols [wlj2001]  
      nbnlocal ='OHd'
      c2temp = nbnlocal(1:2)
      cbnlocal = c2temp
      canlocal = c2temp
      c3temp = nbnlocal(1:3)
      ctnlocal = c3temp
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

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
c     (OHp, OHd, OHt) O in phenol,diols,triols or polyols [wlj2001]   
      nbnlocal ='OHt'
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
      
c     (P ) dimethylphosphate anion [wlj2001]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.740d0
      nbclocal(2) = (0.200d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 30.974d0
      elelocal = ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P'
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
c     * bond potentials                                                *
c     ******************************************************************

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
      nvslocal = 1
c     --- C - OH [wlj1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'OH'
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
      nvslocal = 4
c     --- CA - CT [wlj1996]
      vnlocal(1,1) = 'CA'
      vnlocal(1,2) = 'CT'
c     --- CM - CT [wlj1996]
      vnlocal(2,1) = 'CM'
      vnlocal(2,2) = 'CT'
c     --- CB - CT [wlj1996]
      vnlocal(3,1) = 'CB'
      vnlocal(3,2) = 'CT'
c     --- C= - CT [wlj1996]
      vnlocal(4,1) = 'C='
      vnlocal(4,2) = 'CT'
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
      nvslocal = 1
c     --- CC - CT [wlj1996]
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CT'
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
      nvslocal = 1
c     --- CC - CW [wlj1996]
      vnlocal(1,1) = 'CC'
      vnlocal(1,2) = 'CW'
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
      nvslocal = 3
c     --- CT - CT [jacs1996]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'CT'
c     --- CT - CO [jcc1997]
      vnlocal(2,1) = 'CT'
      vnlocal(2,2) = 'CO'
c     --- CO - CO added for diols and triols
      vnlocal(3,1) = 'CO'
      vnlocal(3,2) = 'CO'
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
      nvslocal = 4
c     --- C* - HC  [wlj1996]
      vnlocal(1,1) = 'C*'
      vnlocal(1,2) = 'HC'
c     --- CK - HA  [wlj1996]
      vnlocal(2,1) = 'CK'
      vnlocal(2,2) = 'HA'
c     --- CM - HC  [wlj1996]
      vnlocal(3,1) = 'CM'
      vnlocal(3,2) = 'HC'
c     --- C= - HC  [wlj1996]
      vnlocal(4,1) = 'C='
      vnlocal(4,2) = 'HC'
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
      nvslocal = 2
c     --- CM - CM  [wlj1996]
      vnlocal(1,1) = 'CM'
      vnlocal(1,2) = 'CM'
c     --- CM - C=  [wlj2001]
      vnlocal(2,1) = 'CM'
      vnlocal(2,2) = 'C='
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

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CW - C= [wlj2001]
      vnlocal(1,1) = 'CW'
      vnlocal(1,2) = 'C='
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.460d0
      vclocal(1) = (385.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C= - C= [wlj2001]
      vnlocal(1,1) = 'C='
      vnlocal(1,2) = 'C='
c     --- C= - C  [wlj2001]
      vnlocal(2,1) = 'C='
      vnlocal(2,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.445d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- NT-NT [wlj2001]
      isame = isame + 1
      vnlocal(isame,1) = 'NT'
      vnlocal(isame,2) = 'NT'
c     --- set nvsame
      nvslocal = isame
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
      nbslocal = 21
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
c     --- CM - C= - HC [wlj2001]
      bnlocal(19,1) = 'CM'
      bnlocal(19,2) = 'C='
      bnlocal(19,3) = 'HC'
c     --- C= - CM - HC [wlj2001]
      bnlocal(20,1) = 'C='
      bnlocal(20,2) = 'CM'
      bnlocal(20,3) = 'HC'
c     --- C= - C= - HC [wlj2001]
      bnlocal(21,1) = 'C='
      bnlocal(21,2) = 'C='
      bnlocal(21,3) = 'HC'
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
      nbslocal = 10
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
c     --- CO - CO - OH [jcc1997]
      bnlocal(6,1) = 'CO'
      bnlocal(6,2) = 'CO'
      bnlocal(6,3) = 'OH'
c     --- CA - NT - CT [jacs1999]
      bnlocal(7,1) = 'CA'
      bnlocal(7,2) = 'NT'
      bnlocal(7,3) = 'CT'
c     --- CT - CT - N* [wlj1996]
      bnlocal(8,1) = 'CT'
      bnlocal(8,2) = 'CT'
      bnlocal(8,3) = 'N*'
c     --- N* - CT - OS [wlj1996]
      bnlocal(9,1) = 'N*'
      bnlocal(9,2) = 'CT'
      bnlocal(9,3) = 'OS'
c     --- CT - CT - F [wlj1996]
      bnlocal(10,1) = 'CT'
      bnlocal(10,2) = 'CT'
      bnlocal(10,3) = 'F'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.0d0
      bclocal(1) = (40.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HC - CT - F [wlj2001]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'F'
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
      ii = 0
c     --- C - CT - N2 [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'C'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'N2'
c     --- C - CT - NT [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'C'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'NT'
c     --- CT - CT - NT [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'CT'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'NT'
c     --- set nbsame
      nbslocal = ii
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
      nbslocal = 1
c     --- C - OH - HO [wlj1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'HO'
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
      nbslocal = 13
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
      nbslocal = 2
c     --- CA - CT - CT [wlj1996]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
c     --- CB - CT - CT [wlj1996]
      bnlocal(2,1) = 'CB'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'CT'
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
      ii = 0 
c     --- HC - CT - OS [jcc1997]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'OS'
c     --- HC - CT - OH [jcc1997]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'OH'
c     --- HC - CO - OS [jcc1997]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CO'
      bnlocal(ii,3) = 'OS'
c     --- HC - CO - OH [jcc1997]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CO'
      bnlocal(ii,3) = 'OH'
c     --- CT - NT - H [jacs1999]
      ii = ii + 1
      bnlocal(ii,1) = 'CT'
      bnlocal(ii,2) = 'NT'
      bnlocal(ii,3) = 'H'
c     --- NT - NT - H [wlj pers. comm. 11-23-2004 by M.G. Martin]
      ii = ii + 1
      bnlocal(ii,1) = 'NT'
      bnlocal(ii,2) = 'NT'
      bnlocal(ii,3) = 'H'
c     --- HC - CT - NT [jacs1999]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'NT'
c     --- C - CT - HC [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'C'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'HC'
c     --- C* - CT - HC [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'C*'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'HC'
c     --- CA - CT - HC [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'CA'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'HC'
c     --- CM - CT - HC [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'CM'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'HC'
c     --- HC - CT - N [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'N'
c     --- HC - CT - N* [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'N*'
c     --- HC - CT - NA [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'NA'
c     --- HC - CT - N2 [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'N2'
c     --- HC - CT - S [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'S'
c     --- HC - CT - SH [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'HC'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'SH'
c     --- CB - CT - HC [wlj1996]
      ii = ii + 1
      bnlocal(ii,1) = 'CB'
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'HC'
c     --- C= - CT - HC [wlj2001]
      ii = ii + 1
      bnlocal(ii,1) = 'C='
      bnlocal(ii,2) = 'CT'
      bnlocal(ii,3) = 'HC'
c     --- set nbname
      nbslocal = ii
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
      nbslocal = 4
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
c     --- C= - C  - CT [wlj1996]
      bnlocal(4,1) = 'C='
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'CT'
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
      nbslocal = 2
C     --- HC - CT - HC [jacs1996] 
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'HC'
C     --- HC - CO - HC [jacs1996]
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'CO'
      bnlocal(2,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.70d0
      bclocal(1) = (37.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
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
C     --- CO - CO - HC [jcc1997]
      bnlocal(4,1) = 'CO'
      bnlocal(4,2) = 'CO'
      bnlocal(4,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.70d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
C     --- CT - CT - CT [jacs1996]
      bnlocal(1,1) = 'CT'
      bnlocal(1,2) = 'CT'
      bnlocal(1,3) = 'CT'
C     --- CT - CT - CO [jcc1997]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'CO'
C     --- CO - CT - CO [jcc1997]
      bnlocal(3,1) = 'CO'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'CO'
C     --- CT - CO - CT [jcc1997]
      bnlocal(4,1) = 'CT'
      bnlocal(4,2) = 'CO'
      bnlocal(4,3) = 'CT'
C     --- CO - CO - CT [jcc1997]
      bnlocal(5,1) = 'CO'
      bnlocal(5,2) = 'CO'
      bnlocal(5,3) = 'CT'
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
      nbslocal = 6
C     --- CA - CA - NC [jms1998]
      bnlocal(1,1) = 'CA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'NC'
C     --- CM - CM - CT [wlj1996]
      bnlocal(2,1) = 'CM'
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'CT'
C     --- CM - C= - C= [wlj2001]
      bnlocal(3,1) = 'CM'
      bnlocal(3,2) = 'C='
      bnlocal(3,3) = 'C='
C     --- CM - C= - CT [wlj2001]
      bnlocal(4,1) = 'CM'
      bnlocal(4,2) = 'C='
      bnlocal(4,3) = 'CT'
C     --- C= - C= - CT [wlj2001]
      bnlocal(5,1) = 'C='
      bnlocal(5,2) = 'C='
      bnlocal(5,3) = 'CT'
C     --- CT - CM - C= [wlj2001]
      bnlocal(6,1) = 'CT'
      bnlocal(6,2) = 'CM'
      bnlocal(6,3) = 'C='
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
C     --- HA - CW - C= [wlj2001]
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CW'
      bnlocal(3,3) = 'C='
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
      nbslocal = 1
C     --- NA - CW - CV [jpcb1998]
      bnlocal(1,1) = 'NA'
      bnlocal(1,2) = 'CW'
      bnlocal(1,3) = 'CV'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.00d0
      bclocal(1) = (70.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- CR - NB - CV [jpcb1998]
      bnlocal(1,1) = 'CR'
      bnlocal(1,2) = 'NB'
      bnlocal(1,3) = 'CV'
C     --- OS - CW - C= furan [wlj2001]
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'CW'
      bnlocal(2,3) = 'C='
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
      nbslocal = 2
C     --- HA - CA - CW [wlj1996]
      bnlocal(1,1) = 'HA'
      bnlocal(1,2) = 'CA'
      bnlocal(1,3) = 'CW'
C     --- HC - C= - CW [wlj2001]
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'C='
      bnlocal(2,3) = 'CW'
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
      bclocal(0) = 123.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
C     --- C= - CM - OS [wlj2001]
      bnlocal(1,1) = 'C='
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'OS'
C     --- C= - CM - OH [wlj2001]
      bnlocal(2,1) = 'C='
      bnlocal(2,2) = 'CM'
      bnlocal(2,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.0d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C= - C  - O [wlj2001]
      bnlocal(1,1) = 'C='
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.0d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- C= - C  - HC [wlj2001]
      bnlocal(1,1) = 'C='
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.0d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CW - C= - C= [wlj2001]
      bnlocal(1,1) = 'CW'
      bnlocal(1,2) = 'C='
      bnlocal(1,3) = 'C='
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.7d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
C     --- CM - C= - C  [wlj2001]
      bnlocal(1,1) = 'CM'
      bnlocal(1,2) = 'C='
      bnlocal(1,3) = 'C'
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
c     --- It appears many of the torsional terms have been changed in 
c     --- the new version of OPLS (05/01).  Below are the values from 
c     --- the publically avaiable version of the parameters at that date

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.300d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 25
c     HC-CT-CT- HC hydrocarbons [wlj2001]
c     --- HC~ - CT - CT - HC~ 
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
c     HC-C -C -HC [wlj2001]
c     --- HC~ - CM - CM - HC~ [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'CM'
      tnlocal(2,4) = 'HC~'
c     HC-C=-C=-HC [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'C='
      tnlocal(3,3) = 'C='
      tnlocal(3,4) = 'HC~'
c     HC-CY-CY-HC small ring [wlj2001]
c     --- HC~ - CY - CY - HC~
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CY'
      tnlocal(4,3) = 'CY'
      tnlocal(4,4) = 'HC~'
c     --- HC~ - CO - CT - HC~ added for alcohols
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CO'
      tnlocal(5,3) = 'CT'
      tnlocal(5,4) = 'HC~'
c     --- HC~ - CO - CO - HC~ added for alcohols
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CO'
      tnlocal(6,3) = 'CO'
      tnlocal(6,4) = 'HC~'
c     HC-CT-CT-CT hydrocarbon [wlj2001]
c     HC~ - CT - CT - CT
      tnlocal(7,1) = 'HC~'
      tnlocal(7,2) = 'CT'
      tnlocal(7,3) = 'CT'
      tnlocal(7,4) = 'CT'
c     HC-CT-CT-CO acetal [wlj2001]
c     --- HC~ - CT - CT - CO
      tnlocal(8,1) = 'HC~'
      tnlocal(8,2) = 'CT'
      tnlocal(8,3) = 'CT'
      tnlocal(8,4) = 'CO'
c     HC-CO-CT-CT acetal [wlj2001]
      tnlocal(9,1) = 'HC~'
      tnlocal(9,2) = 'CO'
      tnlocal(9,3) = 'CT'
      tnlocal(9,4) = 'CT'
c     HC-CO-CT-CO added for alcohols
      tnlocal(10,1) = 'HC~'
      tnlocal(10,2) = 'CO'
      tnlocal(10,3) = 'CT'
      tnlocal(10,4) = 'CO'
c     HC-CO-CO-CT added for alcohols
      tnlocal(11,1) = 'HC~'
      tnlocal(11,2) = 'CO'
      tnlocal(11,3) = 'CO'
      tnlocal(11,4) = 'CT'
c     HC-CT-CO-CT added for alcohols
      tnlocal(12,1) = 'HC~'
      tnlocal(12,2) = 'CT'
      tnlocal(12,3) = 'CO'
      tnlocal(12,4) = 'CT'
c     HC-CT-CO-CO added for alcohols
      tnlocal(13,1) = 'HC~'
      tnlocal(13,2) = 'CT'
      tnlocal(13,3) = 'CO'
      tnlocal(13,4) = 'CO'
c     HC-CT-CY-CY small ring [wlj2001]
      tnlocal(14,1) = 'HC~'
      tnlocal(14,2) = 'CT'
      tnlocal(14,3) = 'CY'
      tnlocal(14,4) = 'CY'
c     HC-CY-CY-CY small ring [wlj2001]
      tnlocal(15,1) = 'HC~'
      tnlocal(15,2) = 'CY'
      tnlocal(15,3) = 'CY'
      tnlocal(15,4) = 'CY'
c     P -OS-CT-HC phosphonates [wlj2001]
      tnlocal(16,1) = 'P'
      tnlocal(16,2) = 'OS'
      tnlocal(16,3) = 'CT'
      tnlocal(16,4) = 'HC~'
c     CT-NT-OS-CT     generic hydroxylamines [wlj2001]
      tnlocal(17,1) = 'CT'
      tnlocal(17,2) = 'CT'
      tnlocal(17,3) = 'OS'
      tnlocal(17,4) = 'CT'
c     CT-NT-OH-HO     generic [wlj2001]
      tnlocal(18,1) = 'CT'
      tnlocal(18,2) = 'NT~'
      tnlocal(18,3) = 'OHm'
      tnlocal(18,4) = 'HO'
c     H -NT-OS-CT     generic [wlj2001]
      tnlocal(19,1) = 'H'
      tnlocal(19,2) = 'NT~'
      tnlocal(19,3) = 'OS'
      tnlocal(19,4) = 'CT'
c     H -NT-OH-HO     generic [wlj2001]
      tnlocal(20,1) = 'H'
      tnlocal(20,2) = 'NT~'
      tnlocal(20,3) = 'OHm'
      tnlocal(20,4) = 'HO'
c     CT-NT-NT-CT     generic hydrazines [wlj2001]
      tnlocal(21,1) = 'CT'
      tnlocal(21,2) = 'NT~'
      tnlocal(21,3) = 'NT~'
      tnlocal(21,4) = 'CT'
c     CT-NT-NT-H      generic [wlj2001]
      tnlocal(22,1) = 'CT'
      tnlocal(22,2) = 'NT~'
      tnlocal(22,3) = 'NT~'
      tnlocal(22,4) = 'H'
c     H -NT-NT-H      generic [wlj2001]
      tnlocal(23,1) = 'H'
      tnlocal(23,2) = 'NT~'
      tnlocal(23,3) = 'NT~'
      tnlocal(23,4) = 'H'
c     HC-C=-C=-HC [wlj2001]
      tnlocal(24,1) = 'HC~'
      tnlocal(24,2) = 'C='
      tnlocal(24,3) = 'C='
      tnlocal(24,4) = 'HC~'
c     CT-C=-C=-HC [wlj2004]
      tnlocal(25,1) = 'CT'
      tnlocal(25,2) = 'C='
      tnlocal(25,3) = 'C='
      tnlocal(25,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.318d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CM-CT-HC alkene [wlj2001]
c     --- HC~ - CM - CT - HC~
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.366d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     HC-CT-CT-CZ alkyne, nitrile [wlj2001]
c     --- HC~ - CT - CT - CZ
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CZ'
c     HC-CT-CT-CM alkene [wlj2001]
c     --- HC~ - CT - CT - CM
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'CM'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.300d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.050d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.200d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     CT-CT-CT-CT hydrocarbon [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     CT-CT-CT-CO hydrocarbon [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'CO'
c     CO-CT-CT-CO hydrocarbon [wlj2001]
      tnlocal(3,1) = 'CO'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'CO'
c     CO-CT-CO-CT hydrocarbon [wlj2001]
      tnlocal(4,1) = 'CO'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CO'
      tnlocal(4,4) = 'CT'
c     CO-CO-CT-CT hydrocarbon [wlj2001]
      tnlocal(5,1) = 'CO'
      tnlocal(5,2) = 'CO'
      tnlocal(5,3) = 'CT'
      tnlocal(5,4) = 'CT'
c     CO-CO-CT-CO hydrocarbon [wlj2001]
      tnlocal(6,1) = 'CO'
      tnlocal(6,2) = 'CO'
      tnlocal(6,3) = 'CT'
      tnlocal(6,4) = 'CO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.372d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     HC-CT-CM=CM alkenes all-atom
c     --- HC~ - CT - CM - CM 
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'CM'
c     CM-C=-C=-CT diene - generic [wlj2001]
      tnlocal(2,1) = 'CM'
      tnlocal(2,2) = 'C='
      tnlocal(2,3) = 'C='
      tnlocal(2,4) = 'CT'
c     HC-CT-CM-C=     alkenes all-atom [wlj2004]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CM'
      tnlocal(3,4) = 'C='
c     HC-CT-C=-C=     alkenes all-atom [wlj2004]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'C='
      tnlocal(4,4) = 'C='
c     HC-CT-C=-CM     alkenes all-atom [wlj2004]
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'C='
      tnlocal(5,4) = 'CM'
c     CM-C=-C=-HC     diene - generic [wlj2001]
      tnlocal(6,1) = 'CM'
      tnlocal(6,2) = 'C='
      tnlocal(6,3) = 'C='
      tnlocal(6,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.468d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 11
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CT - CO - OHm (mono)
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CO'
      tnlocal(1,4) = 'OHm'
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CT - CO - OHp (poly) 
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CO'
      tnlocal(2,4) = 'OHp'
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CT - CO - OHd (diol) 
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CO'
      tnlocal(3,4) = 'OHd'
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CT - CO - OHt (triol) 
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CO'
      tnlocal(4,4) = 'OHt'
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CO - CO - OHp (poly) 
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CO'
      tnlocal(5,3) = 'CO'
      tnlocal(5,4) = 'OHp'
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CO - CO - OHd (diol) 
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CO'
      tnlocal(6,3) = 'CO'
      tnlocal(6,4) = 'OHd'
c     HC-CT-CT-OH alchohols, ethers [wlj2001]
c     --- HC~ - CO - CO - OHt (triol) 
      tnlocal(7,1) = 'HC~'
      tnlocal(7,2) = 'CO'
      tnlocal(7,3) = 'CO'
      tnlocal(7,4) = 'OHt'
c     HC-CT-CT-OS alcohols, ethers [wlj2001]
c     --- HC~ - CT - CT - OS 
      tnlocal(8,1) = 'HC~'
      tnlocal(8,2) = 'CT'
      tnlocal(8,3) = 'CT'
      tnlocal(8,4) = 'OS'
c     HC-CT-CO-OS alcohols, ethers [wlj2001]
      tnlocal(9,1) = 'HC~'
      tnlocal(9,2) = 'CT'
      tnlocal(9,3) = 'CO'
      tnlocal(9,4) = 'OS'
c     HC-CM-CT-O? alcohols, ethers AA [wlj2001]
      tnlocal(10,1) = 'HC~'
      tnlocal(10,2) = 'CM'
      tnlocal(10,3) = 'CT'
      tnlocal(10,4) = 'OHm'
c     HC-CM-CT-O? alcohols, ethers AA [wlj2001]
      tnlocal(11,1) = 'HC~'
      tnlocal(11,2) = 'CM'
      tnlocal(11,3) = 'CT'
      tnlocal(11,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.711d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.663d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     CT-CT-CT-O? alchohols, ethers [wlj2001]
c     --- CT - CT - CT - OHm
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
c     CT-CT-CT-O? alcohols, ethers [wlj2001]
c     --- CT - CT - CT - OS
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OS'
c     CT-CT-CO-O? alcohols, ethers [wlj2001]
c     --- CT - CT - CO - OHm
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CO'
      tnlocal(3,4) = 'OHm'
c     CT-CT-CO-O? alcohols, ethers [wlj2001]
c     --- CT - CT - CO - OS
      tnlocal(4,1) = 'CT'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CO'
      tnlocal(4,4) = 'OS'
c     CT-CM-CT-O? alcohols, ethers AA [wlj2001]
      tnlocal(5,1) = 'CT'
      tnlocal(5,2) = 'CM'
      tnlocal(5,3) = 'CT'
      tnlocal(5,4) = 'OHm'
c     CT-CM-CT-O? alcohols, ethers AA [wlj2001]
      tnlocal(6,1) = 'CT'
      tnlocal(6,2) = 'CM'
      tnlocal(6,3) = 'CT'
      tnlocal(6,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (9.508d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     OH-CT-CT-OH diols only AA [wlj2001]
      tnlocal(1,1) = 'OHd'
      tnlocal(1,2) = 'CO'
      tnlocal(1,3) = 'CO'
      tnlocal(1,4) = 'OHd'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (12.234d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     OH-CT-CT-OH triols only AA [wlj2001]
      tnlocal(1,1) = 'OHt'
      tnlocal(1,2) = 'CO'
      tnlocal(1,3) = 'CO'
      tnlocal(1,4) = 'OHt'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.760d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     HC-CT-OS-C? ethers AA [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'Cwild'
c     HC-CT-OS-CA ethers AA [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CA~'
c     HC-CT-OS-CM ethers AA [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'OS'
      tnlocal(3,4) = 'CM'
c     HC-CM-OS-CT ethers AA [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CM'
      tnlocal(4,3) = 'OS'
      tnlocal(4,4) = 'CT'
c     HC-CO-OS-CT ethers AA [wlj2001]
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CO'
      tnlocal(5,3) = 'OS'
      tnlocal(5,4) = 'CT'
c     HC-CT-OS-CO ethers AA [wlj2001]
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'OS'
      tnlocal(6,4) = 'CO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.650d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.250d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.670d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     CT-OS-CT-CT ethers AA [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     CT-OS-CM-CT ethers AA [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'OS'
      tnlocal(2,3) = 'CM'
      tnlocal(2,4) = 'CT'
c     CT-OS-CO-CT ethers AA [wlj2001]
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'OS'
      tnlocal(3,3) = 'CO'
      tnlocal(3,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.521d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-2.018d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.996d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-OS-CT-O? acetals AA [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
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
      tclocal(2) =  (2.80d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-CW-NA-?? [wlj2001]
c     --- note that this same wildcard pair occurs three times in the 
c     --- file.   this is the one listed first as OPLS torsion 166
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CW'
      tnlocal(1,3) = 'NA'
      tnlocal(1,4) = 'wild'
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
c     --- ??-NA-CW-?? [wlj2001]
c     --- note that this same wildcard pair occurs three times in the 
c     --- file.   this is the one listed second as OPLS torsion 171
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
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-OS-CA-CA anisole [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
c     ??-CW-NA-?? [wlj2001]
c     --- note that this same wildcard pair occurs three times in the
c     --- file.  this is the one listed third  as OPLS torsion 177
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CW'
      tnlocal(2,3) = 'NA'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.423d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.055d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.858d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     CM-C=-C=-CM diene C=C-C=C [wlj2001]
      isame = isame + 1
      tnlocal(isame,1) = 'CM'
      tnlocal(isame,2) = 'C='
      tnlocal(isame,3) = 'C='
      tnlocal(isame,4) = 'CM'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

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
      ntslocal = 90
c     H-C-C[ar]-C[ar] ethyl benzene, toluene [wlj2001]
c     --- HC - CT - CA - CA 
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
c     HC-CT-CW-?? aromatics [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CW'
      tnlocal(2,4) = 'wild'
c     HC-CT-CV-?? aromatics [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CV'
      tnlocal(3,4) = 'wild'
c     HC-CT-CR-?? aromatics [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CR'
      tnlocal(4,4) = 'wild'
c     HC-CT-CS-?? aromatics [wlj2001]
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'CS'
      tnlocal(5,4) = 'wild'
c     HC-CT-CQ-?? aromatics [wlj2001]
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'CQ'
      tnlocal(6,4) = 'wild'
c     HC-CT-CU-?? aromatics [wlj2001]
      tnlocal(7,1) = 'HC~'
      tnlocal(7,2) = 'CT'
      tnlocal(7,3) = 'CU'
      tnlocal(7,4) = 'wild'
c     HC-CT-CK-?? aromatics [wlj2001]
      tnlocal(8,1) = 'HC~'
      tnlocal(8,2) = 'CT'
      tnlocal(8,3) = 'CK'
      tnlocal(8,4) = 'wild'
c     HC-CT-C*-?? aromatics [wlj2001]
      tnlocal(9,1) = 'HC~'
      tnlocal(9,2) = 'CT'
      tnlocal(9,3) = 'C*'
      tnlocal(9,4) = 'wild'
c     CT-CT-CA-CA ethyl benzene [wlj2001]
      tnlocal(10,1) = 'CT'
      tnlocal(10,2) = 'CT'
      tnlocal(10,3) = 'CA'
      tnlocal(10,4) = 'CA'
c     CT-CT-CW-?? aromatics [wlj2001]
      tnlocal(11,1) = 'CT'
      tnlocal(11,2) = 'CT'
      tnlocal(11,3) = 'CW'
      tnlocal(11,4) = 'wild'
c     CT-CT-CV-?? aromatics [wlj2001]
      tnlocal(12,1) = 'CT'
      tnlocal(12,2) = 'CT'
      tnlocal(12,3) = 'CV'
      tnlocal(12,4) = 'wild'
c     CT-CT-CR-?? aromatics [wlj2001]
      tnlocal(13,1) = 'CT'
      tnlocal(13,2) = 'CT'
      tnlocal(13,3) = 'CR'
      tnlocal(13,4) = 'wild'
c     CT-CT-CS-?? aromatics [wlj2001]
      tnlocal(14,1) = 'CT'
      tnlocal(14,2) = 'CT'
      tnlocal(14,3) = 'CS'
      tnlocal(14,4) = 'wild'
c     CT-CT-CQ-?? aromatics [wlj2001]
      tnlocal(15,1) = 'CT'
      tnlocal(15,2) = 'CT'
      tnlocal(15,3) = 'CQ'
      tnlocal(15,4) = 'wild'
c     CT-CT-CU-?? aromatics [wlj2001]
      tnlocal(16,1) = 'CT'
      tnlocal(16,2) = 'CT'
      tnlocal(16,3) = 'CU'
      tnlocal(16,4) = 'wild'
c     CT-CT-CK-?? aromatics [wlj2001]
      tnlocal(17,1) = 'CT'
      tnlocal(17,2) = 'CT'
      tnlocal(17,3) = 'CK'
      tnlocal(17,4) = 'wild'
c     CT-CT-C*-?? aromatics [wlj2001]
      tnlocal(18,1) = 'CT'
      tnlocal(18,2) = 'CT'
      tnlocal(18,3) = 'C*'
      tnlocal(18,4) = 'wild'
c     F - CT-CW-?? fluoromethyl benzene [wlj2001]
      tnlocal(19,1) = 'F'
      tnlocal(19,2) = 'CT'
      tnlocal(19,3) = 'CA~'
      tnlocal(19,4) = 'CA~'
c     Cl-CT-CA-CA chloromethyl benzene [wlj2001]
      tnlocal(20,1) = 'Cl'
      tnlocal(20,2) = 'CT'
      tnlocal(20,3) = 'CA~'
      tnlocal(20,4) = 'CA~'
c     F -CT-CW-?? fluoromethyl aromatic [wlj2001]
      tnlocal(21,1) = 'F'
      tnlocal(21,2) = 'CT'
      tnlocal(21,3) = 'CW'
      tnlocal(21,4) = 'wild'
c     F -CT-CS-?? fluoromethyl aromatic [wlj2001]
      tnlocal(22,1) = 'F'
      tnlocal(22,2) = 'CT'
      tnlocal(22,3) = 'CS'
      tnlocal(22,4) = 'wild'
c     Cl-CT-CW-?? chloromethyl aromatic [wlj2001]
      tnlocal(23,1) = 'Cl'
      tnlocal(23,2) = 'CT'
      tnlocal(23,3) = 'CW'
      tnlocal(23,4) = 'wild'
c     Cl-CT-CS-?? chloromethyl aromatic [wlj2001]
      tnlocal(24,1) = 'Cl'
      tnlocal(24,2) = 'CT'
      tnlocal(24,3) = 'CS'
      tnlocal(24,4) = 'wild'
c     O?-CT-CA-CA benzyl alcohols & ethers [wlj2001]
      tnlocal(24,1) = 'OHm'
      tnlocal(24,2) = 'CT'
      tnlocal(24,3) = 'CA~'
      tnlocal(24,4) = 'CA~'
c     O?-CT-CA-CA benzyl alcohols & ethers [wlj2001]
      tnlocal(25,1) = 'OS'
      tnlocal(25,2) = 'CT'
      tnlocal(25,3) = 'CA~'
      tnlocal(25,4) = 'CA~'
c     C?-CT-NA-C? heterocycles [wlj2001]
      tnlocal(26,1) = 'Cwild'
      tnlocal(26,2) = 'CT'
      tnlocal(26,3) = 'NA'
      tnlocal(26,4) = 'Cwild'
c     H?-CT-NA-C? heterocycles [wlj2001]
      tnlocal(27,1) = 'Hwild'
      tnlocal(27,2) = 'CT'
      tnlocal(27,3) = 'NA'
      tnlocal(27,4) = 'Cwild'
c     C?-CT-N*-C? heterocycles [wlj2001]
      tnlocal(28,1) = 'Cwild'
      tnlocal(28,2) = 'CT'
      tnlocal(28,3) = 'N*'
      tnlocal(28,4) = 'Cwild'
c     H?-CT-N*-C? heterocycles [wlj2001]
      tnlocal(29,1) = 'Hwild'
      tnlocal(29,2) = 'CT'
      tnlocal(29,3) = 'N*'
      tnlocal(29,4) = 'Cwild'
c     CA-CA-CT-N? aromatics [wlj2001]
      tnlocal(30,1) = 'CA~'
      tnlocal(30,2) = 'CA~'
      tnlocal(30,3) = 'CT'
      tnlocal(30,4) = 'Nwild'
c     C -N -CT-HC Phi" peptides AA [wlj2001]
      tnlocal(31,1) = 'C n'
      tnlocal(31,2) = 'N'
      tnlocal(31,3) = 'CT'
      tnlocal(31,4) = 'HC~'
c     HC-CT-C -N Phi" peptides AA [wlj2001]
      tnlocal(32,1) = 'HC'
      tnlocal(32,2) = 'CT'
      tnlocal(32,3) = 'C n'
      tnlocal(32,4) = 'N'
c     HC-CM-C -N [wlj2001]
      tnlocal(33,1) = 'HC~'
      tnlocal(33,2) = 'CM'
      tnlocal(33,3) = 'C n'
      tnlocal(33,4) = 'N'
c     H -N -CT-?? peptides [wlj2001]
      tnlocal(34,1) = 'H'
      tnlocal(34,2) = 'N'
      tnlocal(34,3) = 'CT'
      tnlocal(34,4) = 'wild'
c     H -N -CT-CT peptides [wlj2001]
      tnlocal(35,1) = 'H'
      tnlocal(35,2) = 'N'
      tnlocal(35,3) = 'CT'
      tnlocal(35,4) = 'CT'
c     H -N -CT-C peptides [wlj2001]
      tnlocal(36,1) = 'H'
      tnlocal(36,2) = 'N'
      tnlocal(36,3) = 'CT'
      tnlocal(36,4) = 'C n'
c     ??-CT-C -O peptides [wlj2001]
      tnlocal(37,1) = 'wild'
      tnlocal(37,2) = 'CT'
      tnlocal(37,3) = 'C n'
      tnlocal(37,4) = 'O'
c     CT-CT-C -O peptides [wlj2001]
      tnlocal(38,1) = 'CT'
      tnlocal(38,2) = 'CT'
      tnlocal(38,3) = 'C n'
      tnlocal(38,4) = 'O'
c     N -CT-C -O peptides [wlj2001]
      tnlocal(39,1) = 'N'
      tnlocal(39,2) = 'CT'
      tnlocal(39,3) = 'C n'
      tnlocal(39,4) = 'O'
c     CT-N -CT-HC tert. amide [wlj2001]
      tnlocal(40,1) = 'CT'
      tnlocal(40,2) = 'N'
      tnlocal(40,3) = 'CT'
      tnlocal(40,4) = 'HC~'
c     CA-N -CT-HC tert. amide [wlj2001]
      tnlocal(41,1) = 'CA'
      tnlocal(41,2) = 'N'
      tnlocal(41,3) = 'CT'
      tnlocal(41,4) = 'HC~'
c     HC-CT-N2-H methylguanidinum [wlj2001]
      tnlocal(42,1) = 'HC~'
      tnlocal(42,2) = 'CT'
      tnlocal(42,3) = 'N2'
      tnlocal(42,4) = 'H'
c     HC-CT-N -H N-methylformamide [wlj2001]
      tnlocal(43,1) = 'HC~'
      tnlocal(43,2) = 'CT'
      tnlocal(43,3) = 'N'
      tnlocal(43,4) = 'H'
c     H -N -CT-CT N-ethylformamide [wlj2001]
      tnlocal(44,1) = 'H'
      tnlocal(44,2) = 'N'
      tnlocal(44,3) = 'CT'
      tnlocal(44,4) = 'CT'
c     HC-CT-C -O all carbonyls [wlj2001]
      tnlocal(45,1) = 'HC~'
      tnlocal(45,2) = 'CT'
      tnlocal(45,3) = 'C k'
      tnlocal(45,4) = 'O'
c     HC-C -C -O [wlj2001]
      tnlocal(46,1) = 'HC~'
      tnlocal(46,2) = 'C k'
      tnlocal(46,3) = 'C k'
      tnlocal(46,4) = 'O'
c     HC-CM-C -O [wlj2001]
      tnlocal(47,1) = 'HC~'
      tnlocal(47,2) = 'CM'
      tnlocal(47,3) = 'C k'
      tnlocal(47,4) = 'O'
c     HC-CT-C*-CB 3-methylindole [wlj2001]
      tnlocal(48,1) = 'HC~'
      tnlocal(48,2) = 'CT'
      tnlocal(48,3) = 'C*'
      tnlocal(48,4) = 'CB'
c     CT-CT-C*-CB 3-ethylindole [wlj2001]
      tnlocal(49,1) = 'CT'
      tnlocal(49,2) = 'CT'
      tnlocal(49,3) = 'C*'
      tnlocal(49,4) = 'CB'
c     HC-CT-NO-ON nitro compounds [wlj2001]
      tnlocal(50,1) = 'HC~'
      tnlocal(50,2) = 'CT'
      tnlocal(50,3) = 'NO'
      tnlocal(50,4) = 'ON'
c     HC-CT-C -O2 caboxylates [wlj2001]
      tnlocal(51,1) = 'HC~'
      tnlocal(51,2) = 'CT'
      tnlocal(51,3) = 'C k'
      tnlocal(51,4) = 'O2'
c     HC-CT-C -OH acids [wlj2001]
      tnlocal(52,1) = 'HC~'
      tnlocal(52,2) = 'CT'
      tnlocal(52,3) = 'C a'
      tnlocal(52,4) = 'OHa'
c     CT-CT-C -F acyl halide [wlj2001]
      tnlocal(53,1) = 'CT'
      tnlocal(53,2) = 'CT'
      tnlocal(53,3) = 'C k'
      tnlocal(53,4) = 'F'
c     CT-CT-C -Cl acyl halide [wlj2001]
      tnlocal(54,1) = 'CT'
      tnlocal(54,2) = 'CT'
      tnlocal(54,3) = 'C k'
      tnlocal(54,4) = 'Cl'
c     CT-CT-C -Br acyl halide [wlj2001]
      tnlocal(55,1) = 'CT'
      tnlocal(55,2) = 'CT'
      tnlocal(55,3) = 'C k'
      tnlocal(55,4) = 'Br'
c     CT-CT-C -HC aldehyde [wlj2001]
      tnlocal(56,1) = 'CT'
      tnlocal(56,2) = 'CT'
      tnlocal(56,3) = 'C k'
      tnlocal(56,4) = 'HC'
c     O2-P -OS-CT phosphonates [wlj2001]
      tnlocal(57,1) = 'O2'
      tnlocal(57,2) = 'P'
      tnlocal(57,3) = 'OS'
      tnlocal(57,4) = 'CT'
c     O -P -OS-CA phosphonates [wlj2001]
      tnlocal(58,1) = 'O'
      tnlocal(58,2) = 'P'
      tnlocal(58,3) = 'OS'
      tnlocal(58,4) = 'CA~'
c     O -P -OH-HO phosphonates [wlj2001]
      tnlocal(59,1) = 'O'
      tnlocal(59,2) = 'P'
      tnlocal(59,3) = 'OH'
      tnlocal(59,4) = 'HO'
c     OH-P -OH-HO phosphonates [wlj2001]
      tnlocal(60,1) = 'OH'
      tnlocal(60,2) = 'P'
      tnlocal(60,3) = 'OH'
      tnlocal(60,4) = 'HO'
c     CA-CA-CT-P phosphonates [wlj2001]
      tnlocal(61,1) = 'CA~'
      tnlocal(61,2) = 'CA~'
      tnlocal(61,3) = 'CT'
      tnlocal(61,4) = 'P'
c     CA-CT-P -O2 phosphonates [wlj2001]
      tnlocal(62,1) = 'CA~'
      tnlocal(62,2) = 'CT'
      tnlocal(62,3) = 'P'
      tnlocal(62,4) = 'O2'
c     ??-??-SY-OY sulfonamide [wlj2001]
      tnlocal(63,1) = 'wild'
      tnlocal(63,2) = 'wild'
      tnlocal(63,3) = 'SY'
      tnlocal(63,4) = 'OY'
c     ??-CA-SY-OY sulfonamide [wlj2001]
      tnlocal(64,1) = 'wild'
      tnlocal(64,2) = 'CA~'
      tnlocal(64,3) = 'SY'
      tnlocal(64,4) = 'OY'
c     ??-N -SY-OY sulfonamide [wlj2001]
      tnlocal(65,1) = 'wild'
      tnlocal(65,2) = 'N'
      tnlocal(65,3) = 'SY'
      tnlocal(65,4) = 'OY'
c     ??-CT-SY-OY sulfonamide [wlj2001]
      tnlocal(66,1) = 'wild'
      tnlocal(66,2) = 'CT'
      tnlocal(66,3) = 'SY'
      tnlocal(66,4) = 'OY'
c     ??-CT-SY-N sulfonamide [wlj2001]
      tnlocal(67,1) = 'wild'
      tnlocal(67,2) = 'CT'
      tnlocal(67,3) = 'SY'
      tnlocal(67,4) = 'N'
c     CT-CT-SY-CT     sulfone [wlj2001]
      tnlocal(68,1) = 'CT'
      tnlocal(68,2) = 'CT'
      tnlocal(68,3) = 'SY'
      tnlocal(68,4) = 'CT'
c     CT-CT-SY-OY     sulfone [wlj2001]
      tnlocal(69,1) = 'CT'
      tnlocal(69,2) = 'CT'
      tnlocal(69,3) = 'SY'
      tnlocal(69,4) = 'OY'
c     ??-CT-CZ-NZ     nitriles [wlj2001]
      tnlocal(70,1) = 'wild'
      tnlocal(70,2) = 'CT'
      tnlocal(70,3) = 'CZ'
      tnlocal(70,4) = 'NZ'
c     ??-CZ-CZ-??'     alkynes [wlj2001]
      tnlocal(71,1) = 'wild'
      tnlocal(71,2) = 'CZ'
      tnlocal(71,3) = 'CZ'
      tnlocal(71,4) = 'wild'
c     ??-CT-CZ-CZ     alkynes [wlj2001]
      tnlocal(72,1) = 'wild'
      tnlocal(72,2) = 'CT'
      tnlocal(72,3) = 'CZ'
      tnlocal(72,4) = 'CZ'
c     CT-CT-CZ-CZ     alkynes [wlj2001]
      tnlocal(73,1) = 'CT'
      tnlocal(73,2) = 'CT'
      tnlocal(73,3) = 'CZ'
      tnlocal(73,4) = 'CZ'
c     HC-CT-CZ-CZ     alkynes [wlj2001]
      tnlocal(74,1) = 'HC~'
      tnlocal(74,2) = 'CT'
      tnlocal(74,3) = 'CZ'
      tnlocal(74,4) = 'CZ'
c     HC-CZ-CZ-CT     alkynes [wlj2001]
      tnlocal(75,1) = 'HC~'
      tnlocal(75,2) = 'CZ'
      tnlocal(75,3) = 'CZ'
      tnlocal(75,4) = 'CT'
c     CT-CZ-CZ-CT     alkynes [wlj2001]
      tnlocal(76,1) = 'CT'
      tnlocal(76,2) = 'CZ'
      tnlocal(76,3) = 'CZ'
      tnlocal(76,4) = 'CT'
c     ??-CY-CY-??     small ring [wlj2001]
      tnlocal(77,1) = 'wild'
      tnlocal(77,2) = 'CY'
      tnlocal(77,3) = 'CY'
      tnlocal(77,4) = 'wild'
c     ??-CY-N$-??'     small ring [wlj2001]
      tnlocal(78,1) = 'wild'
      tnlocal(78,2) = 'CY'
      tnlocal(78,3) = 'N$'
      tnlocal(78,4) = 'wild'
c     ??-CY-C$-??     small ring [wlj2001]
      tnlocal(79,1) = 'wild'
      tnlocal(79,2) = 'CY'
      tnlocal(79,3) = 'C$'
      tnlocal(79,4) = 'wild'
c     CY-CY-C$-N$     small ring [wlj2001]
      tnlocal(80,1) = 'CY'
      tnlocal(80,2) = 'CY'
      tnlocal(80,3) = 'C$'
      tnlocal(80,4) = 'N$'
c     CY-CY-C$-O      small ring [wlj2001]
      tnlocal(81,1) = 'CY'
      tnlocal(81,2) = 'CY'
      tnlocal(81,3) = 'C$'
      tnlocal(81,4) = 'O'
c     ??-CY-S -??     small ring [wlj2001]
      tnlocal(82,1) = 'wild'
      tnlocal(82,2) = 'CY'
      tnlocal(82,3) = 'S'
      tnlocal(82,4) = 'wild'
c     CY-CY-N$-??     small ring [wlj2001]
      tnlocal(83,1) = 'CY'
      tnlocal(83,2) = 'CY'
      tnlocal(83,3) = 'N$'
      tnlocal(83,4) = 'wild'
c     HC-CY-N$-C?     small ring [wlj2001]
      tnlocal(84,1) = 'HC~'
      tnlocal(84,2) = 'CY'
      tnlocal(84,3) = 'N$'
      tnlocal(84,4) = 'Cwild'
c     HC-CY-C$-O      small ring [wlj2001]
      tnlocal(85,1) = 'HC~'
      tnlocal(85,2) = 'CY'
      tnlocal(85,3) = 'C$'
      tnlocal(85,4) = 'O'
c     HC-CY-C$-N$     small ring [wlj2001]
      tnlocal(86,1) = 'HC~'
      tnlocal(86,2) = 'CY'
      tnlocal(86,3) = 'C$'
      tnlocal(86,4) = 'N$'
c     N -CY-CY-??     small ring [wlj2001]
      tnlocal(87,1) = 'N'
      tnlocal(87,2) = 'CY'
      tnlocal(87,3) = 'CY'
      tnlocal(87,4) = 'wild'
c     CY-CY-N -H      small ring [wlj2001]
      tnlocal(88,1) = 'CY'
      tnlocal(88,2) = 'CY'
      tnlocal(88,3) = 'N'
      tnlocal(88,4) = 'H'
c     C$-CY-N -H '     small ring [wlj2001]
      tnlocal(89,1) = 'C$'
      tnlocal(89,2) = 'CY'
      tnlocal(89,3) = 'N'
      tnlocal(89,4) = 'H'
c     HC-CY-N -C '     small ring [wlj2001]
      tnlocal(90,1) = 'HC~'
      tnlocal(90,2) = 'CY'
      tnlocal(90,3) = 'N'
      tnlocal(90,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- more null Torsions with all parameters set to 0.0
c     --- needed for the 1-4 interactions 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 9
c     HC-CY-N -H      small ring [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CY'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'H'
c     C?-N$-CT-C?     small ring [wlj2001]
      tnlocal(2,1) = 'Cwild'
      tnlocal(2,2) = 'N$'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'Cwild'
c     C?-N$-CT-HC     small ring [wlj2001]
      tnlocal(3,1) = 'Cwild'
      tnlocal(3,2) = 'N$'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'HC~'
c     HC-NC-NZ-NZ'     azides  [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'NC'
      tnlocal(4,3) = 'NZ'
      tnlocal(4,4) = 'NZ'
c     CT-NC-NZ-NZ'     azides [wlj2001]
      tnlocal(5,1) = 'CT'
      tnlocal(5,2) = 'NC'
      tnlocal(5,3) = 'NZ'
      tnlocal(5,4) = 'NZ'
c     CT-CT-NC-NZ'     azides [wlj2001]
      tnlocal(6,1) = 'CT'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'NC'
      tnlocal(6,4) = 'NZ'
c     HC-CT-C -S='     thiocarbonyl [wlj2001]
      tnlocal(7,1) = 'HC~'
      tnlocal(7,2) = 'CT'
      tnlocal(7,3) = 'C k'
      tnlocal(7,4) = 'S='
c     N -C -C -O '     dicarbonyls  [wlj2001]
      tnlocal(8,1) = 'N'
      tnlocal(8,2) = 'C k'
      tnlocal(8,3) = 'C k'
      tnlocal(8,4) = 'O'
c     HC-CT-C -O all carbonyls [wlj2001]
      tnlocal(9,1) = 'HC~'
      tnlocal(9,2) = 'CT'
      tnlocal(9,3) = 'C a'
      tnlocal(9,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.462d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 9
c     HC-CT-CT-CA ethyl benzene [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
c     HC-CT-CT-CW aromatics [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'CW'
c     HC-CT-CT-CV aromatics [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'CV'
c     HC-CT-CT-CR aromatics [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'CR'
c     HC-CT-CT-CS aromatics [wlj2001]
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'CT'
      tnlocal(5,4) = 'CS'
c     HC-CT-CT-CQ aromatics [wlj2001]
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'CT'
      tnlocal(6,4) = 'CQ'
c     HC-CT-CT-CU aromatics [wlj2001]
      tnlocal(7,1) = 'HC~'
      tnlocal(7,2) = 'CT'
      tnlocal(7,3) = 'CT'
      tnlocal(7,4) = 'CU'
c     HC-CT-CT-CK aromatcis [wlj2001]
      tnlocal(8,1) = 'HC~'
      tnlocal(8,2) = 'CT'
      tnlocal(8,3) = 'CT'
      tnlocal(8,4) = 'CK'
c     HC-CT-CT-C* aromatics [wlj2001]
      tnlocal(9,1) = 'HC~'
      tnlocal(9,2) = 'CT'
      tnlocal(9,3) = 'CT'
      tnlocal(9,4) = 'C*'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-S -S -CT disulfide all-atom [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     HC-CT-S -S disulfide all-atom [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-S -S disulfide all-atom [wlj2001]
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
      tclocal(3) =  (0.480d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-SH-HS thiol all-atom [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.7590)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.282d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.680d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-SH-HS thiol all-atom [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.452d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     HC-CT-CT-SH thiol all-atom [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     HC-CT-CT-S sulfide all-atom [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'S'
c     HC-CT-CT-SY [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'SY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.262d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.198d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.465d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-CT-CT-SH thiol all-atom [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     CT-CT-CT-SY [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'SY'
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
      ntslocal = 7
c     HC-CT-NT-H amina all-atom [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'H'
c     HC-CT-CT-Cl alkyl chloride [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'Cl'
c     CT-CT-CT-Cl alkyl chloride [wlj2001]
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'Cl'
c     HC-CT-CT-F alkyl fluoride [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CTf'
      tnlocal(4,4) = 'F'
c     HC-CT-CT-F alkyl fluoride [wlj2001]
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CTf'
      tnlocal(5,3) = 'CTf'
      tnlocal(5,4) = 'Fpf'
c     HC-CT-CT-Br alkyl bromide [wlj2001]
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'CT'
      tnlocal(6,4) = 'Br'
c     CT-CT-CT-Br alkyl bromide [wlj2001]
      tnlocal(7,1) = 'CT'
      tnlocal(7,2) = 'CT'
      tnlocal(7,3) = 'CT'
      tnlocal(7,4) = 'Br'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.013d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.709d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.473d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-CT-NT amine all-atom [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'NT~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.190d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.417d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.418d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     CT-CT-NT-H amine all-atom [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'H'
c     CT-CT-N2-H guanidinum [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'N2'
      tnlocal(2,4) = 'H'
c     CT-CT-N2-H3 guanidinium [wlj2001]
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'N2'
      tnlocal(3,4) = 'H3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.392d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.674d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.550d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-CT-NT amine all-atom [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'NT~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.682d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HO-OH-CA-CA phenol all-atom [wlj2001]
      tnlocal(1,1) = 'HO'
      tnlocal(1,2) = 'OHm'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     N-C-C-C peptides Chi-1 
      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.845d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.962d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.713d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-1.200d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     N -CT-CT-CT Chi-1 peptides AA [wlj2001]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     N -CT-CT-C* Chi-1 peptides [wlj2001]
      tnlocal(2,1) = 'N'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'C*'
c     N -CT-CT-CA Chi-1 peptides AA [wlj2001]
      tnlocal(3,1) = 'N'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'CA~'
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
      ntslocal = 2
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
      ntslocal = 1
c     --- N - CT - C_n - N  [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'N'
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

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     CT-C -N -H amides [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'H'
c     CM-C -N -H amides [wlj2001]
      tnlocal(2,1) = 'CM'
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'H'
c     HC-C -N -H amides [wlj2001]
      tnlocal(3,1) = 'HC'
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'H'
c     O -C -N -H amides wlj 6/20/97 [wlj2001]
      tnlocal(4,1) = 'O'
      tnlocal(4,2) = 'C n'
      tnlocal(4,3) = 'N'
      tnlocal(4,4) = 'H'
c     H -N -C -C '     dicarbonyls  [wlj2001]
      tnlocal(5,1) = 'H'
      tnlocal(5,2) = 'N'
      tnlocal(5,3) = 'C k'
      tnlocal(5,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.300d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.089d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 9
c     CT-C -N -CT amides -V1 changed to 2.3 [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
c     CT-C- N -CA [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'CA~'
c     HC-C -N -CT [wlj2001]
      tnlocal(3,1) = 'HC'
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'CT'
c     CY-C$-N$-CT'     small ring amides  [wlj2001]
      tnlocal(4,1) = 'CY'
      tnlocal(4,2) = 'C$'
      tnlocal(4,3) = 'N$'
      tnlocal(4,4) = 'CT'
c     CY-C$-N$-CY'     small ring amides [wlj2001]
      tnlocal(5,1) = 'CY'
      tnlocal(5,2) = 'C$'
      tnlocal(5,3) = 'N$'
      tnlocal(5,4) = 'CY'
c     CY-N -CT-CT'     small ring amides [wlj2001]
      tnlocal(6,1) = 'CY'
      tnlocal(6,2) = 'N'
      tnlocal(6,3) = 'CT'
      tnlocal(6,4) = 'CT'
c     CY-C$-N$-CA     small ring amides [wlj2001]
      tnlocal(7,1) = 'CY'
      tnlocal(7,2) = 'C$'
      tnlocal(7,3) = 'N$'
      tnlocal(7,4) = 'CA~'
c     OS-C -N -CT'     carbamates      [wlj2001]
      tnlocal(8,1) = 'OS'
      tnlocal(8,2) = 'C n'
      tnlocal(8,3) = 'N'
      tnlocal(8,4) = 'CT'
c     OH-C -N -CT'     carbamates [wlj2001]
      tnlocal(9,1) = 'OH'
      tnlocal(9,2) = 'C n'
      tnlocal(9,3) = 'N'
      tnlocal(9,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.089d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 6
c     O -C -N -CT amides [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
c     O -C -N -CA amides [wlj2001]
      tnlocal(2,1) = 'O'
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'CA~'
c     O -C -N -CM amides [wlj2001]
      tnlocal(3,1) = 'O'
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'CM'
c     O -C$-N$-CT'     small ring amides [wlj2001]
      tnlocal(4,1) = 'O'
      tnlocal(4,2) = 'C$'
      tnlocal(4,3) = 'N$'
      tnlocal(4,4) = 'CT'
c     O -C$-N$-CA'     small ring amides [wlj2001]
      tnlocal(5,1) = 'O'
      tnlocal(5,2) = 'C$'
      tnlocal(5,3) = 'N$'
      tnlocal(5,4) = 'CA~'
c     O -C -N -CY'     small ring amides [wlj2001]
      tnlocal(6,1) = 'O'
      tnlocal(6,2) = 'C n'
      tnlocal(6,3) = 'N'
      tnlocal(6,4) = 'CY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.1d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     C?-CA-SH-HS aromatic thiol [wlj2001]
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
c     CA-CA-SH-HS aromatic thiol [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'CA~'
      tnlocal(2,3) = 'SH'
      tnlocal(2,4) = 'HS'
c     CA-CA-C -N      aryl amides [wlj2001]
      tnlocal(3,1) = 'CA~'
      tnlocal(3,2) = 'CA~'
      tnlocal(3,3) = 'C n'
      tnlocal(3,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-3.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     N?-CA-SH-HS aromatic thiol [wlj2001]
      tnlocal(1,1) = 'Nwild'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'SH'
      tnlocal(1,4) = 'HS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.647d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     HC-CT-S -CT sulfide all-atom
c     --- HC - CT - S - CT [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
c     HC-CT-S -CA sulfide all-atom [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'S'
      tnlocal(2,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-CT-S sulfide all-atom
c     --- CT - CT - CT - S [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-S -CT sulfide all-atom
c     --- CT - CT - S -CT [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.9d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     H -N2-CA-N2 guanidinium ion [wlj2001]
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'N2'
      tnlocal(1,3) = 'CA+'
      tnlocal(1,4) = 'N2'
c     H -N2-CA-CA guanidinium ion [wlj2001]
      tnlocal(2,1) = 'H'
      tnlocal(2,2) = 'N2'
      tnlocal(2,3) = 'CA+'
      tnlocal(2,4) = 'CA+'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (7.936d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-N2-CA-N2 methyl guanidinium [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'N2'
      tnlocal(1,3) = 'CA+'
      tnlocal(1,4) = 'N2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.177d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-N2-CA methylguanidinium [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N2'
      tnlocal(1,4) = 'CA+'
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
c     CT-CT-N2-CA ethylguanidinum ion [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N2'
      tnlocal(1,4) = 'CA+'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.582d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-CT-N2 ethylguanidinum ion [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'N2'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     C -N -CT-HC N-methylformamide [wlj2001]
c     --- C_n - N - CT - HC~ [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     N -CT-CT-HC N-ethylformamide
c     --- N - CT - CT - HC~ [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     N -CT-CT-CT N-propylformamide
c     --- N - CT - CT - CT [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CM-CM-C -N vinyl amides [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-C -N  propanamide
c     --- CT - CT - C_n - N [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-C -O propanamide
c     --- CT - CT - C_n - O [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C n'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.100d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     HC-CT-CT-C all carbonyls
c     --- HC~ - CT - CT - C_n [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C n'
c     HC-CT-CT-C all carbonyls
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     C -CT-CT-CT  butanamide
c     --- CT - CT - CT - C_n [jacs1996]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C n'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     N -CT-CT-OH Chi for Ser & Thr
c     --- N - CT - CT - OH  [jacs1996]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     C -CT-CT-OH  Chi for Ser & Thr
c     --- C_n - CT - CT - OH  [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHm'
c     C -CT-CT-OH  Chi for Ser & Thr
c     --- C_a - CT - CT - OH  [jacs1996]
      tnlocal(2,1) = 'C a'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.428d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.086d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.029d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     N -CT-CT-S? Chi for Cys (mod 11/99) [wlj2001]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     N -CT-CT-S? Chi for Cys (mod 11/99) [wlj2001]
      tnlocal(2,1) = 'N'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-4.344d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.714d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     C -CT-CT-S? Chi for Cys (mod 11/99) [wlj2001]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'SH'
c     C -CT-CT-S? Chi for Cys (mod 11/99) [wlj2001]
      tnlocal(2,1) = 'C n'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'S'
c     C -CT-CT-S? Chi for Cys (mod 11/99) [wlj2001]
      tnlocal(3,1) = 'C a'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'SH'
c     C -CT-CT-S? Chi for Cys (mod 11/99) [wlj2001]
      tnlocal(4,1) = 'C a'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     HC-CT-C*-CW 3-methylindole [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C*'
      tnlocal(1,4) = 'CW'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.714d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-C*-CW 3-ethylindole [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C*'
      tnlocal(1,4) = 'CW'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.419d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-CC-N? HID, HIE, HIP, 5-ethylimidazole [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CC'
      tnlocal(1,4) = 'Nwild'
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
      ntslocal = 1
c     CT-CT-CC-N? HID, HIE, HIP, 5-ethylimidazole [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CC'
      tnlocal(1,4) = 'Nwild'
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
      ntslocal = 1
c     HC-CT-CT-NO nitroethane [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'NO'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-CT-C carboxylate ion [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C a'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.546d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-CT-C -O RCOOH acid [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'O'
c     CA-CT-C -O RCOOH acid [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'C a'
      tnlocal(2,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.546d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.450d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-C -OH RCOOH acid [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'OHa'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (5.260d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.820d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     NT-CT-C -OH neutral amino acid [wlj2001]
      tnlocal(1,1) = 'NT~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'OHa'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-2.589d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.123d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.270d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-OH-HO axial cyclohexanol [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.550d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -CT-CT-C dicarboxylic acid [wlj2001]
      tnlocal(1,1) = 'C a'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C a'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.5d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     O -C -OH-HO carboxylic acid - aliphatic [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C a'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
c     S=-C -OH-HO'     thioacids    guess             [wlj2001]
      tnlocal(2,1) = 'S='
      tnlocal(2,2) = 'C a'
      tnlocal(2,3) = 'OHa'
      tnlocal(2,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.5d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.5d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-C -OH-HO carboxylic acid - aliphatic [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C a'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
c     HC-C -OH-HO carboxylic acid - alphatic [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'C a'
      tnlocal(2,3) = 'OHa'
      tnlocal(2,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.074d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -CT-CT-HC dicarboxylic acid [wlj2001]
      tnlocal(1,1) = 'C a'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.750d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.550d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.250d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-C -O dicarboxylic acid [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (1.412d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-C -OH dicarboxylic acid [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'OHa'
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
c     CT-N -CT-C Proline (fit to AM1) [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C n'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.753d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.734d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-N -CT-CT Proline (fit to AM1) [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.859d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.058d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-11.266d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-N -CT Proline (fit to AM1) [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-CT-C+-CT carbocation [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C+'
      tnlocal(1,4) = 'CT'
c     CT-CT-C+-HC carbocation [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'C+'
      tnlocal(2,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     HC-CT-C -CT ketone [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.360d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     HC-CT-C -HC aldehyde [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'HC~'
c     HC-CT-C -F acyl halide [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'C k'
      tnlocal(2,4) = 'F'
c     HC-CT-C -Cl acyl halide [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'C k'
      tnlocal(3,4) = 'Cl'
c     HC-CT-C -Br acyl halide [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'C k'
      tnlocal(4,4) = 'Br'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-C -O aldehyde & ketone [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     HC-CT-CT-CO aldehyde & ketone [wlj2001]
c     --- HC~ - CT - CT - C_k [jacs1996]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.697d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.456d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.585d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     CT-CT-CT-C aldehyde & ketone & acyl halide [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'C k'
c     CA-CT-CT-C aldehyde & ketone [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'C k'
c     C*-CT-CT-C [wlj2001]
      tnlocal(3,1) = 'C*'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-C -CT ketone [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (3.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -C -OH-HO oxalic acid, etc. [wlj2001]
      tnlocal(1,1) = 'C a'
      tnlocal(1,2) = 'C a'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.25d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CA-CT-P -OS phosphonates [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'P'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (3.5d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-3.3d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.50d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-P -OS-CT phosphonates [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'P'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.25d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     HC-CT-P -O2 phosphonates [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'P'
      tnlocal(1,4) = 'O2'
c     HC-CT-P -OS phosphonates [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'P'
      tnlocal(2,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.562d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O2-P -OS-C? MeOPO3 (2-) mll [wlj2001]
      tnlocal(1,1) = 'O2'
      tnlocal(1,2) = 'P'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'Cwild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.990d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CA-CA-OS-P PhOPO3 (2-) mll [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'P'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     NT-CT-CT-Cl 2-chloroethylamines [wlj2001]
      tnlocal(1,1) = 'NT~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'Cl'
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
c     Cl-CT-CT-Cl dichloride [wlj2001]
      tnlocal(1,1) = 'Cl'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'Cl'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.478d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-2.1746d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HO-OH-CT-CT trifluoroethanol [wlj2001]
      tnlocal(1,1) = 'HO'
      tnlocal(1,2) = 'OHm'
      tnlocal(1,3) = 'CTf'
      tnlocal(1,4) = 'CTf'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.476d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-OH-HO trifluoroethanol [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CTf'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.450d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     HC-CT-OH-HO alcohols AA [wlj2001]
c     HC-CO-OHm-HO mono AA
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CO'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
c     HC-CO-OHp-HO phenol and poly AA
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CO'
      tnlocal(2,3) = 'OHp'
      tnlocal(2,4) = 'HO'
c     HC-CO-OHd-HO diol AA
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CO'
      tnlocal(3,3) = 'OHd'
      tnlocal(3,4) = 'HO'
c     HC-CO-OHt-HO triol AA
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CO'
      tnlocal(4,3) = 'OHt'
      tnlocal(4,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.3137d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     F -CT-CT-HC trifluoroethanol [wlj2001]
      tnlocal(1,1) = 'F'
      tnlocal(1,2) = 'CTf'
      tnlocal(1,3) = 'CTf'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.5401d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     F -CT-CT-OH trifluoroethanol [wlj2001]
      tnlocal(1,1) = 'F'
      tnlocal(1,2) = 'CTf'
      tnlocal(1,3) = 'CTf'
      tnlocal(1,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.5401d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     F -CT-CT-OH trifluoroethanol [wlj2001]
      tnlocal(1,1) = 'F'
      tnlocal(1,2) = 'CTf'
      tnlocal(1,3) = 'CTf'
      tnlocal(1,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.671d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-4.901d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.669d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     H -N -SY-CA sulfonamide [wlj2001]
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'SY'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.362d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.457d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.149d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-N -SY sulfonamide [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'SY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.074d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-2.966d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (2.473d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-N -SY-CA sulfonamide [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'SY'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.929d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-2.533d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.497d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-N -SY sulfonamide [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'SY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.656d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.768d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.117d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CA-CA-SY-N sulfonamide [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'SY'
      tnlocal(1,4) = 'N'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.300d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.400d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-CT-CT-F alkyl fluoride [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CTf'
      tnlocal(1,4) = 'F'
c     CT-CT-CT-F alkyl fluoride [wlj2001]
      tnlocal(2,1) = 'CTf'
      tnlocal(2,2) = 'CTf'
      tnlocal(2,3) = 'CTf'
      tnlocal(2,4) = 'Fpf'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-2.5d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.250d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     F -CT-CT-F 1,2-difluoride [wlj2001]
      tnlocal(1,1) = 'F'
      tnlocal(1,2) = 'CTf'
      tnlocal(1,3) = 'CTf'
      tnlocal(1,4) = 'F'
c     F -CT-CT-F 1,2-difluoride [wlj2001]
      tnlocal(2,1) = 'F'
      tnlocal(2,2) = 'CTf'
      tnlocal(2,3) = 'CTf'
      tnlocal(2,4) = 'Fpf'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.552d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 9
c     CT-CT-CT-O? polyols, including diols, triols, phenol [wlj2001]
c     --- CT - CT - CO - OHp 
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CO'
      tnlocal(1,4) = 'OHp'
c     --- CT - CT - CO - OHd 
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CO'
      tnlocal(2,4) = 'OHd'
c     --- CT - CT - CO - OHt 
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'CO'
      tnlocal(3,4) = 'OHt'
c     --- CO - CT - CO - OHp 
      tnlocal(4,1) = 'CO'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CO'
      tnlocal(4,4) = 'OHp'
c     --- CO - CT - CO - OHd 
      tnlocal(5,1) = 'CO'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'CO'
      tnlocal(5,4) = 'OHd'
c     --- CO - CT - CO - OHt 
      tnlocal(6,1) = 'CO'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'CO'
      tnlocal(6,4) = 'OHt'
c     --- CT - CO - CO - OHp 
      tnlocal(7,1) = 'CT'
      tnlocal(7,2) = 'CO'
      tnlocal(7,3) = 'CO'
      tnlocal(7,4) = 'OHp'
c     --- CT - CO - CO - OHd 
      tnlocal(8,1) = 'CT'
      tnlocal(8,2) = 'CO'
      tnlocal(8,3) = 'CO'
      tnlocal(8,4) = 'OHd'
c     --- CT - CO - CO - OHt 
      tnlocal(9,1) = 'CT'
      tnlocal(9,2) = 'CO'
      tnlocal(9,3) = 'CO'
      tnlocal(9,4) = 'OHt'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.552d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 7
c     CT-CT-OH-HO alcohols AA
c     --- CT - CO - OHm - HO   
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CO'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
c     --- CT - CO - OHp - HO  
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CO'
      tnlocal(2,3) = 'OHp'
      tnlocal(2,4) = 'HO'
c     --- CT - CO - OHd - HO  
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CO'
      tnlocal(3,3) = 'OHd'
      tnlocal(3,4) = 'HO'
c     --- CT - CO - OHt - HO  
      tnlocal(4,1) = 'CT'
      tnlocal(4,2) = 'CO'
      tnlocal(4,3) = 'OHt'
      tnlocal(4,4) = 'HO'
c     --- CO - CO - OHp - HO  
      tnlocal(5,1) = 'CO'
      tnlocal(5,2) = 'CO'
      tnlocal(5,3) = 'OHp'
      tnlocal(5,4) = 'HO'
c     --- CO - CO - OHd - HO  
      tnlocal(6,1) = 'CO'
      tnlocal(6,2) = 'CO'
      tnlocal(6,3) = 'OHd'
      tnlocal(6,4) = 'HO'
c     --- CO - CO - OHt - HO  
      tnlocal(7,1) = 'CO'
      tnlocal(7,2) = 'CO'
      tnlocal(7,3) = 'OHt'
      tnlocal(7,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-CT-OH-HO     hexopyranoses [wlj2001]
c     --- CT - CT - OHp - HO  [jcc1997]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHp'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     OH-CT-CT-OH hexopyranoses [wlj2001]
c     --- OHp - CT - CT - OHp  [jcc1997]
      tnlocal(1,1) = 'OHp'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OHp'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.375d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.358d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.004d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-OS-CO-OS hexopyranoses [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CO'
      tnlocal(1,4) = 'OS'
c     CT-OS-CO-OH hexopyranoses [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'OS'
      tnlocal(2,3) = 'CO'
      tnlocal(2,4) = 'OHp'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     OS-CO-OH-HO hexopyranoses [wlj2001]
      tnlocal(1,1) = 'OS'
      tnlocal(1,2) = 'CO'
      tnlocal(1,3) = 'OHp'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     OH-CT-CT-OS hexopyranoses [wlj2001]
      tnlocal(1,1) = 'OHp'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.65d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-NA-CR-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CR'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (10.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 8
c     ??-CR-NB-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CR'
      tnlocal(1,3) = 'NB'
      tnlocal(1,4) = 'wild'
c     ??-CW-NB-?? [wlj2001]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CW'
      tnlocal(2,3) = 'NB'
      tnlocal(2,4) = 'wild'
c     ??-CR-NA-?? [wlj2001]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CR'
      tnlocal(3,3) = 'NA'
      tnlocal(3,4) = 'wild'
c     ??-CR-NC-?? [wlj2001]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CR'
      tnlocal(4,3) = 'NC'
      tnlocal(4,4) = 'wild'
c     ??-CK-NB-?? [wlj2001]
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CK'
      tnlocal(5,3) = 'NB'
      tnlocal(5,4) = 'wild'
c     ??-CK-NA-?? [wlj2001]
      tnlocal(6,1) = 'wild'
      tnlocal(6,2) = 'CK'
      tnlocal(6,3) = 'NA'
      tnlocal(6,4) = 'wild'
c     ??-CK-NC-?? [wlj2001]
      tnlocal(7,1) = 'wild'
      tnlocal(7,2) = 'CK'
      tnlocal(7,3) = 'NC'
      tnlocal(7,4) = 'wild'
c     ??-NA-NB-?? [wlj2001]
      tnlocal(8,1) = 'wild'
      tnlocal(8,2) = 'NA'
      tnlocal(8,3) = 'NB'
      tnlocal(8,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.80d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-NB-CV-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'NB'
      tnlocal(1,3) = 'CV'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (10.75d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     ??-CW-CV-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CW'
      tnlocal(1,3) = 'CV'
      tnlocal(1,4) = 'wild'
c     ??-CW-CW-?? [wlj2001]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CW'
      tnlocal(2,3) = 'CW'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (13.05d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-C*-CW-?? 'wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C*'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.35d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-C*-CB-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C*'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (7.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     ??-CA-CB-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'CB'
      tnlocal(1,4) = 'wild'
c     ??-C -CB-?? [wlj2001]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C k'
      tnlocal(2,3) = 'CB'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-CB-CN-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CB'
      tnlocal(1,3) = 'CN'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.05d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ??-NA-CN-?? [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CN'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.1320d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-CT-C -OS     esters [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.669d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.124d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CT-C -OS-CT     esters [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
c     HC-C -OS-CT     esters [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'C k'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.124d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O -C -OS-CT     esters [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.1980d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -OS-CT-HC     esters [wlj2001]
      tnlocal(1,1) = 'C k'
      tnlocal(1,2) = 'OS'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.5530d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-C -OS     esters [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.220d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.126d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.4220d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-OS-C      esters [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'C k'
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
c     CT-CT-NO-ON nitroethane [wlj2001]
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
c     CA-CA-NO-ON nitrobenzene [wlj2001]
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
      tclocal(2) =  (3.651d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-NY-CA-NZ     neutral arg [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'NY'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'NZ'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CT-NT-CT-CT amine [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'NT~'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (11.035d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.968d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.270d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     NT-CT-CT-NT     amine all-atom [wlj2001]
      tnlocal(1,1) = 'NT~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'NT~'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     HC-CT-NT-CT     amine all-atom [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.700d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.600d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-CW-NA     2-ethyl pyrrole [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'CW'
      tnlocal(1,4) = 'NA'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.230d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.505d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     CT-C=-C=-CM 2-Me-1,3-butadiene [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C='
      tnlocal(1,3) = 'C='
      tnlocal(1,4) = 'CM'
c     CT-C=-C -O? 2-Me-1,3-butadiene-like [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'C='
      tnlocal(2,3) = 'C k'
      tnlocal(2,4) = 'Owild'
c     CT-CM-C -O?'     2-Me-1,3-butadiene-like
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CM'
      tnlocal(3,3) = 'C k'
      tnlocal(3,4) = 'Owild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.800d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-3.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-C -C=-CM     methyl vinyl ketone [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C='
      tnlocal(1,4) = 'CM'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (3.200d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-3.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CM-C=-C -OH     acrylic acid [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'C='
      tnlocal(1,3) = 'C a'
      tnlocal(1,4) = 'OHa'
c     CM-CM-C -OH     acrylic acid-like [wlj2001]
      tnlocal(2,1) = 'CM'
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'C a'
      tnlocal(2,4) = 'OHa'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CM-C=-C -O      acrolein [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'C='
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'OHm'
c     CM-CM-C -O      acrolein-like [wlj2001]
      tnlocal(2,1) = 'CM'
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'C k'
      tnlocal(2,4) = 'OHm'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.20d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CA-CA-C -CT     aryl ketone [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
c     CA-CA-C -HC     aryl aldehyde [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'CA~'
      tnlocal(2,3) = 'C k'
      tnlocal(2,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.10d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 7
c     CA-CA-N -??     N-phenylamide [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'wild'
c     CA-CA-N -H '     N-phenylamide [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'CA~'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'H'
c     CA-CA-N -C      N-phenylamide [wlj2001]
      tnlocal(3,1) = 'CA~'
      tnlocal(3,2) = 'CA~'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'C n'
c     CA-CA-N -CT     N-phenylamide[wlj2001]
      tnlocal(4,1) = 'CA~'
      tnlocal(4,2) = 'CA~'
      tnlocal(4,3) = 'N'
      tnlocal(4,4) = 'CT'
c     CA-CA-C -O      aryl carbonyl [wlj2001]
      tnlocal(5,1) = 'CA~'
      tnlocal(5,2) = 'CA~'
      tnlocal(5,3) = 'C k'
      tnlocal(5,4) = 'O'
c     CA-CA-C -OS     aryl acid, ester [wlj2001]
      tnlocal(6,1) = 'CA~'
      tnlocal(6,2) = 'CA~'
      tnlocal(6,3) = 'C k'
      tnlocal(6,4) = 'OS'
c     CA-CA-C -OH     aryl acid, ester [wlj2001]
      tnlocal(7,1) = 'CA~'
      tnlocal(7,2) = 'CA~'
      tnlocal(7,3) = 'C a'
      tnlocal(7,4) = 'OHa'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.50d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CA-CA-OS-C      phenyl acetate [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'C k'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     O -C -OH-HO     benzoic acids [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C a'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
c     O -C -OS-CT     benzoic esters [wlj2001]
      tnlocal(2,1) = 'O'
      tnlocal(2,2) = 'C k'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT'
c     O -C -OS-CA'     phenyl acetate [wlj2001]
      tnlocal(3,1) = 'O'
      tnlocal(3,2) = 'C k'
      tnlocal(3,3) = 'OS'
      tnlocal(3,4) = 'CA~'
c     S=-C -OS-CT'     thioesters   guess [wlj2001]
      tnlocal(4,1) = 'S='
      tnlocal(4,2) = 'C k'
      tnlocal(4,3) = 'OS'
      tnlocal(4,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CA-C -OH-HO     benzoic acids [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'C a'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
c     CA-C -OS-CT     benzoic esters [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'C a'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-C -OS-CA     phenyl acetate [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.030d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     CA-CA-NT-H      aniline [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'H'
c     ??-CA-N2-H      aniline-like [wlj2001]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CA~'
      tnlocal(2,3) = 'N2'
      tnlocal(2,4) = 'H'
c     ??-CQ-N2-H '     aniline-like [wlj2001]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CQ'
      tnlocal(3,3) = 'N2'
      tnlocal(3,4) = 'H'
c     N?-CA-N2-H      aniline-like [wlj2001]
      tnlocal(4,1) = 'Nwild'
      tnlocal(4,2) = 'CA~'
      tnlocal(4,3) = 'N2'
      tnlocal(4,4) = 'H'
c     CB-CA-N2-H      aniline-like [wlj2001]
      tnlocal(5,1) = 'CB'
      tnlocal(5,2) = 'CA~'
      tnlocal(5,3) = 'N2'
      tnlocal(5,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

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
c     CA-CA-NT-CT     substituted-aniline [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'NT~'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.817d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.169d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.543d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CM-CT-CT     alkenes [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.346d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.405d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.904d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CM-CM-CT-CT     alkenes [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     C=-CM-CT-CT     alkenes [wlj2001]
      tnlocal(2,1) = 'C='
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.241d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.353d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (-0.286d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CM-CM-CA-CA     styrene [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.205d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.531d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CM-CA-CA     1-methylstyrene [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
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
      ntslocal = 19
c     ??-CM-CM-??'     alkene [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'wild'
c     CT-CM-CM-CT'     alkene [wlj2001]
      tnlocal(2,1) = 'CT'
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'CM'
      tnlocal(2,4) = 'CT'
c     CT-CM-C=-CT     alkene [wlj2001]
      tnlocal(3,1) = 'CT'
      tnlocal(3,2) = 'CM'
      tnlocal(3,3) = 'C='
      tnlocal(3,4) = 'CT'
c     CT-CM-C=-HC     alkene [wlj2001]
      tnlocal(4,1) = 'CT'
      tnlocal(4,2) = 'CM'
      tnlocal(4,3) = 'C='
      tnlocal(4,4) = 'HC~'
c     HC-CM-C=-HC     alkene [wlj2001]
      tnlocal(5,1) = 'HC~'
      tnlocal(5,2) = 'CM'
      tnlocal(5,3) = 'C='
      tnlocal(5,4) = 'HC~'
c     HC-CM-C=-CT     alkene [wlj2001]
      tnlocal(6,1) = 'HC~'
      tnlocal(6,2) = 'CM'
      tnlocal(6,3) = 'C='
      tnlocal(6,4) = 'CT'
c     CT-CM-CM-HC     alkene [wlj2001]
      tnlocal(7,1) = 'CT'
      tnlocal(7,2) = 'CM'
      tnlocal(7,3) = 'CM'
      tnlocal(7,4) = 'HC~'
c     C -CM-CM-HC     alkene [wlj2001]
      tnlocal(8,1) = 'C k'
      tnlocal(8,2) = 'CM'
      tnlocal(8,3) = 'CM'
      tnlocal(8,4) = 'HC~'
c     HC-CM-CM-HC     alkene [wlj2001]
      tnlocal(9,1) = 'HC~'
      tnlocal(9,2) = 'CM'
      tnlocal(9,3) = 'CM'
      tnlocal(9,4) = 'HC~'
c     Cl-CM-CM-HC     chloroalkene [wlj2001]
      tnlocal(10,1) = 'Cl'
      tnlocal(10,2) = 'CM'
      tnlocal(10,3) = 'CM'
      tnlocal(10,4) = 'HC~'
c     HC-CM-CM-OS     vinyl ether [wlj2001]
      tnlocal(11,1) = 'HC~'
      tnlocal(11,2) = 'CM'
      tnlocal(11,3) = 'CM'
      tnlocal(11,4) = 'OS'
c     CT-CM-CM-OS     vinyl ether [wlj2001]
      tnlocal(12,1) = 'CT'
      tnlocal(12,2) = 'CM'
      tnlocal(12,3) = 'CM'
      tnlocal(12,4) = 'OS'
c     HC-CM-CM-OH     vinyl alcohol [wlj2001]
      tnlocal(13,1) = 'HC~'
      tnlocal(13,2) = 'CM'
      tnlocal(13,3) = 'CM'
      tnlocal(13,4) = 'OHm'
c     CT-CM-CM-OH     vinyl alcohol [wlj2001]
      tnlocal(14,1) = 'CT'
      tnlocal(14,2) = 'CM'
      tnlocal(14,3) = 'CM'
      tnlocal(14,4) = 'OHm'
c     CT-C -NC-CT     imine [wlj2001]
      tnlocal(15,1) = 'CT'
      tnlocal(15,2) = 'C n'
      tnlocal(15,3) = 'NC'
      tnlocal(15,4) = 'CT'
c     HC-C -NC-HC     imine [wlj2001]
      tnlocal(16,1) = 'HC~'
      tnlocal(16,2) = 'C n'
      tnlocal(16,3) = 'NC'
      tnlocal(16,4) = 'HC~'
c     CT-C -NC-O?     oxime [wlj2001]
      tnlocal(17,1) = 'CT'
      tnlocal(17,2) = 'C n'
      tnlocal(17,3) = 'NC'
      tnlocal(17,4) = 'Owild'
c     HC-C -NC-O?     oxime [wlj2001]
      tnlocal(18,1) = 'HC~'
      tnlocal(18,2) = 'C n'
      tnlocal(18,3) = 'NC'
      tnlocal(18,4) = 'Owild'
c     HC-CM-C=-C= [wlj2004]
      tnlocal(19,1) = 'HC~'
      tnlocal(19,2) = 'CM'
      tnlocal(19,3) = 'C='
      tnlocal(19,4) = 'C='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (2.170d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 19
c     CA-C!-C!-CA     biphenyl [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'C!'
      tnlocal(1,3) = 'C!'
      tnlocal(1,4) = 'CA~'
c     CA-C!-C!-NC     biphenyl-like [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'C!'
      tnlocal(2,3) = 'C!'
      tnlocal(2,4) = 'NC'
c     CA-C!-C -CA     biphenyl [wlj2001]
      tnlocal(3,1) = 'CA~'
      tnlocal(3,2) = 'C!'
      tnlocal(3,3) = 'C k'
      tnlocal(3,4) = 'CA~'
c     CA-C!-NA-NB     biphenyl-like [wlj2001]
      tnlocal(4,1) = 'CA~'
      tnlocal(4,2) = 'C!'
      tnlocal(4,3) = 'NA'
      tnlocal(4,4) = 'NB'
c     CA-C!-NA-CW     biphenyl-like [wlj2001]
      tnlocal(5,1) = 'CA~'
      tnlocal(5,2) = 'C!'
      tnlocal(5,3) = 'NA'
      tnlocal(5,4) = 'CW'
c     CA-C!-CW-C?     biphenyl-like [wlj2001]
      tnlocal(6,1) = 'CA~'
      tnlocal(6,2) = 'C!'
      tnlocal(6,3) = 'CW'
      tnlocal(6,4) = 'Cwild'
c     CA-C!-CW-NA     biphenyl-like [wlj2001]
      tnlocal(7,1) = 'CA~'
      tnlocal(7,2) = 'C!'
      tnlocal(7,3) = 'CW'
      tnlocal(7,4) = 'NA'
c     CA-C!-CW-OS     biphenyl-like [wlj2001]
      tnlocal(8,1) = 'CA~'
      tnlocal(8,2) = 'C!'
      tnlocal(8,3) = 'CW'
      tnlocal(8,4) = 'OS'
c     CA-C!-CS-C?     biphenyl-like [wlj2001]
      tnlocal(9,1) = 'CA~'
      tnlocal(9,2) = 'C!'
      tnlocal(9,3) = 'CS'
      tnlocal(9,4) = 'Cwild'
c     CA-C!-CM-C?     biphenyl-like [wlj2001]
      tnlocal(10,1) = 'CA~'
      tnlocal(10,2) = 'C!'
      tnlocal(10,3) = 'CM'
      tnlocal(10,4) = 'Cwild'
c     CA-C!-CU-C?     biphenyl-like [wlj2001]
      tnlocal(11,1) = 'CA~'
      tnlocal(11,2) = 'C!'
      tnlocal(11,3) = 'CU'
      tnlocal(11,4) = 'Cwild'
c     CA-C!-CU-NB     biphenyl-like [wlj2001]
      tnlocal(12,1) = 'CA~'
      tnlocal(12,2) = 'C!'
      tnlocal(12,3) = 'CU'
      tnlocal(12,4) = 'NB'
c     CA-C!-CV-C?     biphenyl-like [wlj2001]
      tnlocal(13,1) = 'CA~'
      tnlocal(13,2) = 'C!'
      tnlocal(13,3) = 'CV'
      tnlocal(13,4) = 'Cwild'
c     CA-C!-CV-NB     biphenyl-like [wlj2001]
      tnlocal(14,1) = 'CA~'
      tnlocal(14,2) = 'C!'
      tnlocal(14,3) = 'CV'
      tnlocal(14,4) = 'NB'
c     CA-C!-CR-OS     biphenyl-like [wlj2001]
      tnlocal(15,1) = 'CA~'
      tnlocal(15,2) = 'C!'
      tnlocal(15,3) = 'CR'
      tnlocal(15,4) = 'OS'
c     CA-C!-CR-NA     biphenyl-like [wlj2001]
      tnlocal(16,1) = 'CA~'
      tnlocal(16,2) = 'C!'
      tnlocal(16,3) = 'CR'
      tnlocal(16,4) = 'NA'
c     CA-C!-CR-NB     biphenyl-like [wlj2001]
      tnlocal(17,1) = 'CA~'
      tnlocal(17,2) = 'C!'
      tnlocal(17,3) = 'CR'
      tnlocal(17,4) = 'NB'
c     CA-C!-CW-S      biphenyl-like [wlj2001]
      tnlocal(18,1) = 'CA~'
      tnlocal(18,2) = 'C!'
      tnlocal(18,3) = 'CW'
      tnlocal(18,4) = 'S'
c     CA-C!-CR-S      biphenyl-like [wlj2001]
      tnlocal(19,1) = 'CA~'
      tnlocal(19,2) = 'C!'
      tnlocal(19,3) = 'CR'
      tnlocal(19,4) = 'S'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.600d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (14.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     Cl-CM-CM-Cl     chloroalkene [wlj2001]
      tnlocal(1,1) = 'Cl'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CM'
      tnlocal(1,4) = 'Cl'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.300d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     HC-CT-P+-CT      phosphonium ion [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'P+'
      tnlocal(1,4) = 'CT'
c     HC-CT-CT-P+      phosphonium ion [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'P+'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.500d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-CT-P+-CT  phosphonium ion [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'P+'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (3.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     C -NC-OH-HO  oxime B3LYP/6-31G* [wlj2001]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'NC'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
c     C -NC-OS-CT     oxime   11/00 [wlj2001]
      tnlocal(2,1) = 'C n'
      tnlocal(2,2) = 'NC'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-3.50d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     CM-CM-OS-CT     vinyl ether [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
c     C=-CM-OS-CT     vinyl ether [wlj2001]
      tnlocal(2,1) = 'C='
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     CM-CM-CT-O?     allyl alcohols, ethers [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'CM'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'Owild'
c     C=-CM-CT-O?     allyl alcohols, ethers [wlj2001]
      tnlocal(2,1) = 'C='
      tnlocal(2,2) = 'CM'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'Owild'
c     CM-C=-CT-O?     allyl alcohols, ethers [wlj2001]
      tnlocal(3,1) = 'CM'
      tnlocal(3,2) = 'C='
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'Owild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.90d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     CM-CT-OH-HO     allyl alcohols [wlj2001]
      tnlocal(1,1) = 'CM'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OHm'
      tnlocal(1,4) = 'HO'
c     CA-CT-OH-HO     benzyl alcohols [wlj2001]
      tnlocal(2,1) = 'CA~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'OHm'
      tnlocal(2,4) = 'HO'
c     C=-CT-OH-HO     allyl alcohols [wlj2001]
      tnlocal(3,1) = 'C='
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'OHm'
      tnlocal(3,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.35d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     HC-CT-SY-OY     sulfone [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'SY'
      tnlocal(1,4) = 'OY'
c     HC-CT-SY-CT     sulfone [wlj2001]
      tnlocal(2,1) = 'HC~'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'SY'
      tnlocal(2,4) = 'CT'
c     HC-CT-SY-CA     sulfone [wlj2001]
      tnlocal(3,1) = 'HC~'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'SY'
      tnlocal(3,4) = 'CA~'
c     HC-CT-CT-SY     sulfone [wlj2001]
      tnlocal(4,1) = 'HC~'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'SY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CA-CA-SY-CT     sulfone 10/00 B3LYP PhSO2Me [wlj2001]
      tnlocal(1,1) = 'CA~'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'SY'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.50d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     CR-NA-CT-OS     imidazoles, indoles, purines [wlj2001]
      tnlocal(1,1) = 'CR'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
c     CR-N*-CT-OS     imidazoles, indoles, purines [wlj2001]
      tnlocal(2,1) = 'CR'
      tnlocal(2,2) = 'N*'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OS'
c     CK-NA-CT-OS'     imidazoles, indoles, purines [wlj2001]
      tnlocal(3,1) = 'CK'
      tnlocal(3,2) = 'NA'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'OS'
c     CK-N*-CT-OS'     imidazoles, indoles, purines [wlj2001]
      tnlocal(4,1) = 'CK'
      tnlocal(4,2) = 'N*'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.50d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     C?-NA-CT-OS     imidazoles, indoles, purines [wlj2001]
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'OS'
c     C?-N*-CT-OS     imidazoles, indoles, purines [wlj2001]
      tnlocal(2,1) = 'Cwild'
      tnlocal(2,2) = 'N*'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'OS'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.350d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     CR-NA-CT-CT     imidazoles, indoles, purines [wlj2001]
      tnlocal(1,1) = 'CR'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     CR-N*-CT-CT     imidazoles, indoles, purines [wlj2001]
      tnlocal(2,1) = 'CR'
      tnlocal(2,2) = 'N*'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'CT'
c     CK-NA-CT-CT     imidazoles, indoles, purines [wlj2001]
      tnlocal(3,1) = 'CK'
      tnlocal(3,2) = 'NA'
      tnlocal(3,3) = 'CT'
      tnlocal(3,4) = 'CT'
c     CK-N*-CT-CT     imidazoles, indoles, purines [wlj2001]
      tnlocal(4,1) = 'CK'
      tnlocal(4,2) = 'N*'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.350d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     C?-NA-CT-CT     imidazoles, indoles, purines [wlj2001]
      tnlocal(1,1) = 'Cwild'
      tnlocal(1,2) = 'NA'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     C?-N*-CT-CT'     imidazoles, indoles, purines [wlj2001]
      tnlocal(2,1) = 'Cwild'
      tnlocal(2,2) = 'N*'
      tnlocal(2,3) = 'CT'
      tnlocal(2,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.50d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     N?-CT-OS-CT     imidazoles, indoles, purines [wlj2001]
      tnlocal(1,1) = 'Nwild'
      tnlocal(1,2) = 'CT'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     O -C$-N$-H      small ring amides [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C$'
      tnlocal(1,3) = 'N$'
      tnlocal(1,4) = 'H'
c     CY-C$-N$-H      small ring amides [wlj2001]
      tnlocal(2,1) = 'CY'
      tnlocal(2,2) = 'C$'
      tnlocal(2,3) = 'N$'
      tnlocal(2,4) = 'H'
c     OS-C -N -H '     carbamates      [wlj2001]
      tnlocal(3,1) = 'OS'
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'H'
c     OH-C -N -H '     carbamates [wlj2001]
      tnlocal(4,1) = 'OHm'
      tnlocal(4,2) = 'C n'
      tnlocal(4,3) = 'N'
      tnlocal(4,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-1.396d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.427d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- C_n - N - CT - CT [jacs1996]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'CT'
c     C -N -CY-CY'     small ring       [wlj2001]
      tnlocal(2,1) = 'C n'
      tnlocal(2,2) = 'N'
      tnlocal(2,3) = 'CY'
      tnlocal(2,4) = 'CY'
c     C -N -CY-C$'     small ring       [wlj2001]
      tnlocal(3,1) = 'C n'
      tnlocal(3,2) = 'N'
      tnlocal(3,3) = 'CY'
      tnlocal(3,4) = 'C$'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (20.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O -C$-N$-CY'     small ring amides [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C$'
      tnlocal(1,3) = 'N$'
      tnlocal(1,4) = 'CY'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-2.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     N -C -OS-CT     carbamates [wlj2001]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'OS'
      tnlocal(1,4) = 'CT'
c     N -C -OH-HO     carbamates [wlj2001]
      tnlocal(2,1) = 'N'
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'OHm'
      tnlocal(2,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (5.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     S=-C -N -H      thioamides   wlj 01/00  [wlj2001]
      tnlocal(1,1) = 'S='
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     S=-C -N -CT     thioamides   fit to [wlj2001]
      tnlocal(1,1) = 'S='
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
c     S=-C -N -CA'     thioamides   MP3/6-31+G** [wlj2001]
      tnlocal(2,1) = 'S='
      tnlocal(2,2) = 'C n'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'CA~'
c     S=-C -N -CM'     thioamides   (Wiberg & Rush) [wlj2001]
      tnlocal(3,1) = 'S='
      tnlocal(3,2) = 'C n'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'CM'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.600d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-S -CA-CA     thioanisole  fit to MP4 [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'CA~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O -C -C -CT'     dicarbonyls Kahn & Bruice  [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.200d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O -C -C -HC'     dicarbonyls BMC 8,1881(2000) [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.800d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-0.760d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-C -C -CT'     dicarbonyls    [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.800d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     HC-C -C -HC'     dicarbonyls    [wlj2001]
      tnlocal(1,1) = 'HC~'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.700d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-1.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.0d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-C -C -CT'     dicarbonyls    [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.085d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -C -CT-HC'     dicarbonyls    [wlj2001]
      tnlocal(1,1) = 'C k'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'CT'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.300d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     N -C -C -HC'     dicarbonyls   [wlj2001]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'HC~'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-0.500d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.200d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     N -C -C -CT'     dicarbonyls  [wlj2001]
      tnlocal(1,1) = 'N'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.400d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (4.900d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -C -N -CT'     dicarbonyls   [wlj2001]
      tnlocal(1,1) = 'C k'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CT'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (1.600d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (3.200d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O -C -C -O '     dicarbonyls  [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C k'
      tnlocal(1,3) = 'C k'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (4.542d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.603d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (1.045d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     CT-C -N -OH     hydroxamic acids   [wlj2001]
      tnlocal(1,1) = 'CT'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'OHa'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (6.603d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     O -C -N -OH     hydroxamic acids [wlj2001]
      tnlocal(1,1) = 'O'
      tnlocal(1,2) = 'C n'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'OHa'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (5.519d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-6.700d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.581d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C -N -OH-HO'     hydroxamic acids [wlj2001]
      tnlocal(1,1) = 'C n'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (2.722d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (-5.154d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     H -N -OH-HO     hydroxamic acids [wlj2001]
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'N'
      tnlocal(1,3) = 'OHa'
      tnlocal(1,4) = 'HO'
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
      ntslocal = 20
c     ??-CA-CA-?? in aromatic ring [wlj2001]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CA~'
      tnlocal(1,3) = 'CA~'
      tnlocal(1,4) = 'wild'
c     ??-CB-CB-?? [wlj2001]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CB'
      tnlocal(2,3) = 'CB'
      tnlocal(2,4) = 'wild'
c     ??-CB-CS-?? [wlj2001]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CB'
      tnlocal(3,3) = 'CS'
      tnlocal(3,4) = 'wild'
c     ??-CS-C -?? [wlj2001]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CS'
      tnlocal(4,3) = 'C k'
      tnlocal(4,4) = 'wild'
c     ??-CR-CS-?? [wlj2001]
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CR'
      tnlocal(5,3) = 'CS'
      tnlocal(5,4) = 'wild'
c     ??-CS-CW-?? [wlj2001]
      tnlocal(6,1) = 'wild'
      tnlocal(6,2) = 'CS'
      tnlocal(6,3) = 'CW'
      tnlocal(6,4) = 'wild'
c     HA-CW-OS-CW [wlj2001]
      tnlocal(7,1) = 'HA'
      tnlocal(7,2) = 'CW'
      tnlocal(7,3) = 'OS'
      tnlocal(7,4) = 'CW'
c     ??-NC-CA-?? [wlj2001]
      tnlocal(8,1) = 'wild'
      tnlocal(8,2) = 'NC'
      tnlocal(8,3) = 'CA~'
      tnlocal(8,4) = 'wild'
c     CS-CW-OS-CW [wlj2001]
      tnlocal(9,1) = 'CS'
      tnlocal(9,2) = 'CW'
      tnlocal(9,3) = 'OS'
      tnlocal(9,4) = 'CW'
c     CA-NC-NC-CT [wlj2001]
      tnlocal(10,1) = 'CA~'
      tnlocal(10,2) = 'NC'
      tnlocal(10,3) = 'NC'
      tnlocal(10,4) = 'CT'
c     CA-NC-NC-CA [wlj2001]
      tnlocal(11,1) = 'CA~'
      tnlocal(11,2) = 'NC'
      tnlocal(11,3) = 'NC'
      tnlocal(11,4) = 'CA~'
c     CT-NC-NC-CT [wlj2001]
      tnlocal(12,1) = 'CT'
      tnlocal(12,2) = 'NC'
      tnlocal(12,3) = 'NC'
      tnlocal(12,4) = 'CT'
c     ??-NC-CB-?? [wlj2001]
      tnlocal(13,1) = 'wild'
      tnlocal(13,2) = 'NC'
      tnlocal(13,3) = 'CB'
      tnlocal(13,4) = 'wild'
c     ??-NA-CB-?? [wlj2001]
      tnlocal(14,1) = 'wild'
      tnlocal(14,2) = 'NA'
      tnlocal(14,3) = 'CB'
      tnlocal(14,4) = 'wild'
c     ??-NB-CB-?? [wlj2001]
      tnlocal(15,1) = 'wild'
      tnlocal(15,2) = 'NB'
      tnlocal(15,3) = 'CB'
      tnlocal(15,4) = 'wild'
c     C?-NB-CB-C? [wlj2001]
      tnlocal(16,1) = 'Cwild'
      tnlocal(16,2) = 'NB'
      tnlocal(16,3) = 'CB'
      tnlocal(16,4) = 'Cwild'
c     ??-N -CB-?? [wlj2001]
      tnlocal(17,1) = 'wild'
      tnlocal(17,2) = 'N'
      tnlocal(17,3) = 'CB'
      tnlocal(17,4) = 'wild'
c     ??-NC-CR-?? [wlj2001]
      tnlocal(18,1) = 'wild'
      tnlocal(18,2) = 'NC'
      tnlocal(18,3) = 'CR'
      tnlocal(18,4) = 'wild'
c     ??-NC-CQ-?? [wlj2001]
      tnlocal(19,1) = 'wild'
      tnlocal(19,2) = 'NC'
      tnlocal(19,3) = 'CQ'
      tnlocal(19,4) = 'wild'
c     ??-N -CQ-?? [wlj2001]
      tnlocal(20,1) = 'wild'
      tnlocal(20,2) = 'N'
      tnlocal(20,3) = 'CQ'
      tnlocal(20,4) = 'wild'
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
      nislocal = 3
c     --- Z - CA - X - Y [wlj2001]
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'CA~'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'wild'
c     --- I have learned that OPLS sometimes applies this improper to 
c     --- other sp2 centers, not just ones with a CA in the middle
c     --- Therefore a handful of other types have been added here
c     --- Their use appears to be optional
c     --- M.G. Martin 07-13-2004 from personal communication with 
c     --- W.L. Jorgensen on 05-20-2004.
c     --- Z - CM - X - Y [mgm2004]
      ninlocal(2,1) = 'wild'
      ninlocal(2,2) = 'CM'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'wild'
c     --- Z - C= - X - Y [mgm2004]
      ninlocal(3,1) = 'wild'
      ninlocal(3,2) = 'C='
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'wild'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,2)
      iclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(2) =  (21.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 3
c     --- Z - C - X - Y [wlj2001]
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
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,1)
      call twh_impstyle(GLB_SET,type,2)
      iclocal(1) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(2) =  (2.000d0)*(CON_KCALTOK)*(0.5d0)
      iclocal(3) =  (0.000d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- Z - N - X - Y [wlj2001]
      ninlocal(1,1) = 'wild'
      ninlocal(1,2) = 'N'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'wild'
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end
