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
      subroutine twh_ffcompassv1(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes,naatypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the compass      *
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
     &        ,nimprotypes,naatypes,nbitypes
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
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)
      integer aanslocal
      character*(FFNAMELEN) aanlocal
      dimension aanlocal(MAXSAMEAA,4)
      double precision aaclocal
      dimension aaclocal(MINAACOEFF:MAXAACOEFF)
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)
      double precision onepi

      ffname = 'COMPASSv1'
      call twh_classical_potential(GLB_SET,'9-6')
      call twh_classical_mixrule(GLB_SET,'Sixth Power')
      type = 0

c     one-four nonbond has no scale factor so just uses same params
c     References:
c     --- H. Sun, J. Phys. Chem. B (1998), 102, 7338-7364.
c     --- D. Rigby, H. Sun, B. E. Eichinger, Polymer International
c     44 (1997) 311-330..
c     --- S.W. Bunte, H. Sun, J. Phys. Chem. B (2000) 104 2477-2489.
c     --- H. Sun, P. Ren, J.R. Fried, Computational and Theoretical
c     Polymer Science, 8, 229-246 (1998).
c     --- J. Yang, Y. Ren, A. Tian, H. Sun, J. Phys. Chem. B 
c     104  4951-4957 (2000).
c     [jcc 2003]
c     M.J. McQuaid; J. Sun; D. Rigby; "Development and Validation of 
c     COMPASS Force Field Parameters for Molecules with Aliphatic
c     Azide Chains", J. Comp. Chem. 25, 61-71 (2003)
      
      
c     --- Parameters 1 to 5 from J. Phys. Chem. B (1998) 102 
c     7338-7364.
      
c     --- (c3a) aromatic carbon
      type = 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.91500d0
      nbclocal(2) = (0.0680d0)*(CON_KCALTOK)
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal) 
      nbnlocal ='c3a'
      cbnlocal ='c3a'
      canlocal ='c3a'
      ctnlocal ='c3a'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      
c     --- (c4) generic sp3 carbon bonded to less than 3 heavy atoms
      type = 2
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.8540d0
      nbclocal(2) = (0.0620d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      nbnlocal ='c4'
      cbnlocal ='c4'
      canlocal ='c4'
      ctnlocal ='c4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (c43) sp3 carbon bonded to exactly 3 heavy atoms (methine)
      type = 3
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.8540d0
      nbclocal(2) = (0.0400d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      nbnlocal ='c43'
      cbnlocal ='c4'
      canlocal ='c4'
      ctnlocal ='c4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (c44) sp3 carbon bonded to exactly 4 heavy atoms 
c     --- (quaternary)
      type = 4
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.8540d0
      nbclocal(2) = (0.0200d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      nbnlocal ='c44'
      cbnlocal ='c4'
      canlocal ='c4'
      ctnlocal ='c4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (h1) hydrogen, nonpolar
      type = 5
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.8780d0
      nbclocal(2) = (0.0230d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 1.0079d0
      nbnlocal ='h1'
      cbnlocal ='h1'
      canlocal ='h1'
      ctnlocal ='h1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- Parameters 6 to 8 from J. Phys. Chem. B (2000) 104 
c     2477-2489.
      
c     --- (n3o) nitrogen in nitro group bonded to 3 atoms
      type = 6
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.7600d0
      nbclocal(2) = (0.04800d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.00674d0
      nbnlocal ='n3o'
      cbnlocal ='n3o'
      canlocal ='n3o'
      ctnlocal ='n3o'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o12) oxygen in nitro group (-NO2) bonded to 1 atom
      type = 7
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.4000d0
      nbclocal(2) = (0.04800d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.99940d0
      nbnlocal ='o12'
      cbnlocal ='o1='
      canlocal ='o1='
      ctnlocal ='o1='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o2n) oxygen in nitrates bonded to 2 atoms
      type = 8
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.6500d0
      nbclocal(2) = (0.20000d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.99940d0
      nbnlocal ='o2n'
      cbnlocal ='o2n'
      canlocal ='o2n'
      ctnlocal ='o2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- Parameters 9 to 10 from Comp. and Theor. Polymer Sci. 
c     (1998)
      
c     --- (n2=)  nitrogen bonded to 2 atoms (from phosphazene)
      type = 9
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.830d0
      nbclocal(2) = (0.096d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.00674d0
      nbnlocal ='n2='
      cbnlocal ='n2='
      canlocal ='n2='
      ctnlocal ='n2='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (p4=)  phosphorous bonded to 4 atoms (from phosphazenes)
      type = 10
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.295d0
      nbclocal(2) = (0.065d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 30.974d0
      nbnlocal ='p4='
      cbnlocal ='p4='
      canlocal ='p4='
      ctnlocal ='p4='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' P'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- Parameters 11 to 28 from J. Phys. Chem. B (2000) 104 
c     4951-4957.
      
c     --- (he) helium
      type = 11
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.9000d0
      nbclocal(2) = (0.00500d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 4.003d0
      nbnlocal ='he'
      cbnlocal ='he'
      canlocal ='he'
      ctnlocal ='he'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'He'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (ne) neon
      type = 12
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.2000d0
      nbclocal(2) = (0.05500d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 20.180d0
      nbnlocal ='ne'
      cbnlocal ='ne'
      canlocal ='ne'
      ctnlocal ='ne'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'Ne'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (ar) argon
      type = 13
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.8800d0
      nbclocal(2) = (0.20000d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 39.948d0
      nbnlocal ='ar'
      cbnlocal ='ar'
      canlocal ='ar'
      ctnlocal ='ar'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'Ar'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (kr) krypton
      type = 14
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.3000d0
      nbclocal(2) = (0.28000d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 83.80d0
      nbnlocal ='kr'
      cbnlocal ='kr'
      canlocal ='kr'
      ctnlocal ='kr'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'Kr'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (xe) xenon
      type = 15
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.2600d0
      nbclocal(2) = (0.39000d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 131.29d0
      nbnlocal ='xe'
      cbnlocal ='xe'
      canlocal ='xe'
      ctnlocal ='xe'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'Xe'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (h1h) hydrogen in H2
      type = 16
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.4210d0
      nbclocal(2) = (0.02160d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 1.0079d0
      nbnlocal ='h1h'
      cbnlocal ='h1h'
      canlocal ='h1h'
      ctnlocal ='h1h'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (n1n) nitrogen in N2
      type = 17
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.8008d0
      nbclocal(2) = (0.05980d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.00674d0
      nbnlocal ='n1n'
      cbnlocal ='n1n'
      canlocal ='n1n'
      ctnlocal ='n1n'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (c1o) carbon in CO
      type = 18 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.0120d0
      nbclocal(2) = (0.05300d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      nbnlocal ='c1o'
      cbnlocal ='c1o'
      canlocal ='c1o'
      ctnlocal ='c1o'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o1o) oxygen in O2
      type = 19
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.4758d0
      nbclocal(2) = (0.07800d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o1o'
      cbnlocal ='o1o'
      canlocal ='o1o'
      ctnlocal ='o1o'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o1c) oxygen in CO
      type = 20
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.6020d0
      nbclocal(2) = (0.08500d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o1c'
      cbnlocal ='o1c'
      canlocal ='o1c'
      ctnlocal ='o1c'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (n1o) nitrogen in NO
      type = 21
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.4600d0
      nbclocal(2) = (0.12800d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.00674d0
      nbnlocal ='n1o'
      cbnlocal ='n1o'
      canlocal ='n1o'
      ctnlocal ='n1o'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o1n) oxygen in NO
      type = 22
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.3000d0
      nbclocal(2) = (0.15600d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o1n'
      cbnlocal ='o1n'
      canlocal ='o1n'
      ctnlocal ='o1n'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (c2=) carbon in CO2 and CS2
      type = 23
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.9150d0
      nbclocal(2) = (0.06800d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      nbnlocal ='c2='
      cbnlocal ='c2='
      canlocal ='c2='
      ctnlocal ='c2='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (s2=) sulfur in SO2
      type = 24
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.0470d0
      nbclocal(2) = (0.12500d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 32.066d0
      nbnlocal ='s2='
      cbnlocal ='s2='
      canlocal ='s2='
      ctnlocal ='s2='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' S'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (n2o) nitrogen in NO2
      type = 25
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.5290d0
      nbclocal(2) = (0.33300d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.00674d0
      nbnlocal ='n2o'
      cbnlocal ='n2o'
      canlocal ='2no'
      ctnlocal ='n2o'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o1=) oxygen in NO2 and SO2
      type = 26
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.4300d0
      nbclocal(2) = (0.19200d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o1='
      cbnlocal ='o1='
      canlocal ='o1='
      ctnlocal ='o1='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o1=*) oxygen in CO2
      type = 27
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.3600d0
      nbclocal(2) = (0.06700d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o1=*'
      cbnlocal ='o1=*'
      canlocal ='o1=*'
      ctnlocal ='o1=*'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (s1=) sulfur in CS2
      type = 28
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.0070d0
      nbclocal(2) = (0.31300d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 32.066d0
      nbnlocal ='s1='
      cbnlocal ='s1='
      canlocal ='s1='
      ctnlocal ='s1='
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' S'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp2'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- Parameters 29 to 33 from Polymer International 44 (1997) 
c     311-330
      
c     --- (c4o) alpha carbon in an alcohol
      type = 29 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.815d0
      nbclocal(2) = (0.068d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      nbnlocal ='c4o'
      cbnlocal ='c4'
      canlocal ='c4'
      ctnlocal ='c4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (h1o) Strongly polar hydrogen bonded to an oxygen
      type = 30
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.087d0
      nbclocal(2) = (0.008d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 1.0079d0
      nbnlocal ='h1o'
      cbnlocal ='h1'
      canlocal ='h1'
      ctnlocal ='h1'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' H'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o2) parameters listed, but atom type not defined
c     --- also listed in spectrochimica Acta where defined as
c     --- (o2) generic oxygen with two bonds attached
      type = 31
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.300d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o2'
      cbnlocal ='o2'
      canlocal ='o2'
      ctnlocal ='o2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 's'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o2e) Ether oxygen
      type = 32
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.300d0
      nbclocal(2) = (0.120d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o2e'
      cbnlocal ='o2e'
      canlocal ='o2'
      ctnlocal ='o2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (o2h) Hydroxyl oxygen
      type = 33         
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.580d0
      nbclocal(2) = (0.096d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 15.999d0
      nbnlocal ='o2h'
      cbnlocal ='o2h'
      canlocal ='o2'
      ctnlocal ='o2'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 

c     --- parameters from Sun Rigby Spectrochimica Acta (1997)

c     --- (si4) generic silicon with four bonds attached
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.405d0
      nbclocal(2) = (0.198d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 28.086d0
      nbnlocal ='si4'
      cbnlocal ='si4'
      canlocal ='si4'
      ctnlocal ='si4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'Si'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 

c     --- (si4c) a subset of Si4, non-hydrogen atom attached
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.290d0
      nbclocal(2) = (0.131d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 28.086d0
      nbnlocal ='si4c'
      cbnlocal ='si4'
      canlocal ='si4'
      ctnlocal ='si4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal = 'Si'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 

c     --- parameters from JCC 2003

c     --- (c4z) carbon, sp3, bonded to N3
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.6500d0
      nbclocal(2) = (0.0800d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 12.011d0
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbnlocal ='c4z'
      cbnlocal ='c4'
      canlocal ='c4'
      ctnlocal ='c4'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal =  ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)

c     --- (n2z) Nitrogen, first atom in N3
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.4000d0
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.007d0
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbnlocal ='n2z'
      cbnlocal ='n2z'
      canlocal ='n2z'
      ctnlocal ='n2z'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal =  ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)

c     --- (n2t) Nitrogen, central atom in N3
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.3000d0
      nbclocal(2) = (0.0500d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.007d0
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbnlocal ='n2t'
      cbnlocal ='n2t'
      canlocal ='n2t'
      ctnlocal ='n2t'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      elelocal =  ' N'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)


c     --- (n1z) Nitrogen, terminal atom in N3
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.5200d0
      nbclocal(2) = (0.0850d0)*(CON_KCALTOK)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      masslocal = 14.007d0
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbnlocal ='n1z'
      cbnlocal ='n1t'
      canlocal ='n1t'
      ctnlocal ='n1t'
      elelocal =  ' N'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Compass bond potentials                                        *
c     ******************************************************************
c     --- H. Sun J. Phys. Chem. B (1998) 102, 7338-7364.
c     --- see ffnonbond for complete list of references
c     --- harmonic, triple and quartic terms

c     --- c3a - c3a 
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4170d0
      vclocal(1) = (470.8361d0)*(CON_KCALTOK)
      vclocal(2) = (-627.6179d0)*(CON_KCALTOK)
      vclocal(3) = (1327.6345d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'c3a'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c3a - c4 
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.5010d0
      vclocal(1) = (321.9021d0)*(CON_KCALTOK)
      vclocal(2) = (-521.8208d0)*(CON_KCALTOK)
      vclocal(3) = (572.1628d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'c4'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c3a - h1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.0982d0
      vclocal(1) = (372.8251d0)*(CON_KCALTOK)
      vclocal(2) = (-803.4526d0)*(CON_KCALTOK)
      vclocal(3) = (894.3173d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'h1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - c4
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.5300d0
      vclocal(1) = (299.6700d0)*(CON_KCALTOK)
      vclocal(2) = (-501.7700d0)*(CON_KCALTOK)
      vclocal(3) = (679.8100d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'c4'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - h1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1010d0
      vclocal(1) = (345.0000d0)*(CON_KCALTOK)
      vclocal(2) = (-691.8900d0)*(CON_KCALTOK)
      vclocal(3) = (844.6000d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'h1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c3a - n3o
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4300d0
      vclocal(1) = (313.8329d0)*(CON_KCALTOK)
      vclocal(2) = (-568.6087d0)*(CON_KCALTOK)
      vclocal(3) = (600.9597d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'n3o'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - n3o
c     --- c4o - n3o
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4740d0
      vclocal(1) = (301.6051d0)*(CON_KCALTOK)
      vclocal(2) = (-535.7028d0)*(CON_KCALTOK)
      vclocal(3) = (555.0420d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'n3o'
      vnlocal(2,1) = 'c4o'
      vnlocal(2,2) = 'n3o'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - o2n
c     --- c4o - o2n
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4350d0
      vclocal(1) = (400.3954d0)*(CON_KCALTOK)
      vclocal(2) = (-835.1951d0)*(CON_KCALTOK)
      vclocal(3) = (1313.0142d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'o2n'
      vnlocal(2,1) = 'c4o'
      vnlocal(2,2) = 'o2n'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- h1 - n3o
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.0400d0
      vclocal(1) = (439.9346d0)*(CON_KCALTOK)
      vclocal(2) = (-943.7307d0)*(CON_KCALTOK)
      vclocal(3) = (1180.9318d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'h1'
      vnlocal(1,2) = 'n3o'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n3o - o1=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.2100d0
      vclocal(1) = (765.0664d0)*(CON_KCALTOK)
      vclocal(2) = (-2070.2830d0)*(CON_KCALTOK)
      vclocal(3) = (2793.3218d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n3o'
      vnlocal(1,2) = 'o1='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n3o - o2n
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4020d0
      vclocal(1) = (300.0000d0)*(CON_KCALTOK)
      vclocal(2) = (-1000.000d0)*(CON_KCALTOK)
      vclocal(3) = (2000.00d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n3o'
      vnlocal(1,2) = 'o2n'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c3a - n2=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.400d0
      vclocal(1) = (350.0000d0)*(CON_KCALTOK)
      vclocal(2) = (0.000d0)*(CON_KCALTOK)
      vclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'n2='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c3a - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.789d0
      vclocal(1) = (197.702d0)*(CON_KCALTOK)
      vclocal(2) = (-332.251d0)*(CON_KCALTOK)
      vclocal(3) = (325.716d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - n2=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.474d0
      vclocal(1) = (337.060d0)*(CON_KCALTOK)
      vclocal(2) = (-147.370d0)*(CON_KCALTOK)
      vclocal(3) = (213.633d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'n2='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.800d0
      vclocal(1) = (218.140d0)*(CON_KCALTOK)
      vclocal(2) = (-329.511d0)*(CON_KCALTOK)
      vclocal(3) = (290.349d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- cl1p - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 2.000d0
      vclocal(1) = (158.777d0)*(CON_KCALTOK)
      vclocal(2) = (-239.129d0)*(CON_KCALTOK)
      vclocal(3) = (210.084d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'cl1p'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- f1p - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.565d0
      vclocal(1) = (340.000d0)*(CON_KCALTOK)
      vclocal(2) = (-882.384d0)*(CON_KCALTOK)
      vclocal(3) = (1197.919d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'f1p'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- h1 - n2=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.031d0
      vclocal(1) = (540.112d0)*(CON_KCALTOK)
      vclocal(2) = (-1500.2952d0)*(CON_KCALTOK)
      vclocal(3) = (2431.008d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'h1'
      vnlocal(1,2) = 'n2='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- h1 - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.430d0
      vclocal(1) = (285.204d0)*(CON_KCALTOK)
      vclocal(2) = (-575.685d0)*(CON_KCALTOK)
      vclocal(3) = (677.846d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'h1'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n2= - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.598d0
      vclocal(1) = (393.006d0)*(CON_KCALTOK)
      vclocal(2) = (-751.405d0)*(CON_KCALTOK)
      vclocal(3) = (767.431d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n2='
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n3 - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.678d0
      vclocal(1) = (329.000d0)*(CON_KCALTOK)
      vclocal(2) = (-713.795d0)*(CON_KCALTOK)
      vclocal(3) = (902.919d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n3'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- o2 - p4=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.600d0
      vclocal(1) = (333.098d0)*(CON_KCALTOK)
      vclocal(2) = (-726.623d0)*(CON_KCALTOK)
      vclocal(3) = (924.620d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'o2'
      vnlocal(1,2) = 'p4='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- h1h - h1h
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 0.7412d0
      vclocal(1) = (414.2185d0)*(CON_KCALTOK)
      vclocal(2) = (-805.6549d0)*(CON_KCALTOK)
      vclocal(3) = (914.1296d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'h1h'
      vnlocal(1,2) = 'h1h'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n1n - n1n
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.0977d0
      vclocal(1) = (1651.3730d0)*(CON_KCALTOK)
      vclocal(2) = (-4069.3178d0)*(CON_KCALTOK)
      vclocal(3) = (5984.9629d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n1n'
      vnlocal(1,2) = 'n1n'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- o1o - o1o
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.2074d0
      vclocal(1) = (846.7150d0)*(CON_KCALTOK)
      vclocal(2) = (-2247.1760d0)*(CON_KCALTOK)
      vclocal(3) = (3478.9900d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'o1o'
      vnlocal(1,2) = 'o1o'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c1o - o1c
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1283d0
      vclocal(1) = (1368.7676d0)*(CON_KCALTOK)
      vclocal(2) = (-3157.0007d0)*(CON_KCALTOK)
      vclocal(3) = (4247.5298d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c1o'
      vnlocal(1,2) = 'o1c'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n1o - o1n
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1506d0
      vclocal(1) = (1147.8362d0)*(CON_KCALTOK)
      vclocal(2) = (-3167.7349d0)*(CON_KCALTOK)
      vclocal(3) = (5099.5811d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n1o'
      vnlocal(1,2) = 'o1n'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- o1= - s2=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4308d0
      vclocal(1) = (730.8387d0)*(CON_KCALTOK)
      vclocal(2) = (-1531.7910d0)*(CON_KCALTOK)
      vclocal(3) = (1859.7753d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'o1='
      vnlocal(1,2) = 's2='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c2= - o1=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1600d0
      vclocal(1) = (1161.3421d0)*(CON_KCALTOK)
      vclocal(2) = (-2564.5706d0)*(CON_KCALTOK)
      vclocal(3) = (3932.8735d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c2='
      vnlocal(1,2) = 'o1='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- n2o - o1=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1930d0
      vclocal(1) = (620.0000d0)*(CON_KCALTOK)
      vclocal(2) = (-1808.6018d0)*(CON_KCALTOK)
      vclocal(3) = (3077.5918d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'n2o'
      vnlocal(1,2) = 'o1='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c2= - s1=
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.5540d0
      vclocal(1) = (559.0065d0)*(CON_KCALTOK)
      vclocal(2) = (-1348.6633d0)*(CON_KCALTOK)
      vclocal(3) = (1248.8604d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c2='
      vnlocal(1,2) = 's1='
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- Parameters from Polymer International (1997) 44 311-330.

c     --- c3a - o2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.3768d0
      vclocal(1) = (428.8798d0)*(CON_KCALTOK)
      vclocal(2) = (-738.235d0)*(CON_KCALTOK)
      vclocal(3) = (1114.9655d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'o2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c3a - o2e
c     --- c3a - o2h
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.3768d0
      vclocal(1) = (428.8798d0)*(CON_KCALTOK)
      vclocal(2) = (-738.2351d0)*(CON_KCALTOK)
      vclocal(3) = (1114.9655d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'c3a'
      vnlocal(1,2) = 'o2e'
      vnlocal(2,1) = 'c3a'
      vnlocal(2,2) = 'o2h'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- c4 - o2e
c     --- c4 - o2h
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4200d0
      vclocal(1) = (400.3954d0)*(CON_KCALTOK)
      vclocal(2) = (-835.1951d0)*(CON_KCALTOK)
      vclocal(3) = (1313.0142d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'c4'
      vnlocal(1,2) = 'o2e'
      vnlocal(2,1) = 'c4'
      vnlocal(2,2) = 'o2h'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- h1 - o2h
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 0.9494d0
      vclocal(1) = (540.3633d0)*(CON_KCALTOK)
      vclocal(2) = (-1311.8663d0)*(CON_KCALTOK)
      vclocal(3) = (2132.4446d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'h1'
      vnlocal(1,2) = 'o2h'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- parameters from Sun; Rigby; Spectrochimica Acta Part A
c     --- 53 (1997) 1301-1323

c     o2 - si4
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.640d0
      vclocal(1) = (350.1232d0)*(CON_KCALTOK)
      vclocal(2) = (-517.3424d0)*(CON_KCALTOK)
      vclocal(3) = (673.7067d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'o2'
      vnlocal(1,2) = 'si4'
      call twh_set_vibnames(type,nvslocal,vnlocal)


c     --- parameters from [JCC 2003]

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1354d0
      vclocal(1) = (1337.7450d0)*(CON_KCALTOK)
      vclocal(2) = (-2675.4900d0)*(CON_KCALTOK)
      vclocal(3) = (3121.4049d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     n1t - n1t [jcc2003]
      vnlocal(1,1) = 'n1t'
      vnlocal(1,2) = 'n1t'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.2343d0
      vclocal(1) = (720.3345d0)*(CON_KCALTOK)
      vclocal(2) = (-1542.6689d0)*(CON_KCALTOK)
      vclocal(3) = (1799.7804d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     n2z - n2t [jcc2003]
      vnlocal(1,1) = 'n2z'
      vnlocal(1,2) = 'n2t'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.1354d0
      vclocal(1) = (1198.7450d0)*(CON_KCALTOK)
      vclocal(2) = (-2675.4900d0)*(CON_KCALTOK)
      vclocal(3) = (3121.4049d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     n2t - n1t [jcc2003]
      vnlocal(1,1) = 'n2t'
      vnlocal(1,2) = 'n1t'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.4814d0
      vclocal(1) = (324.4578d0)*(CON_KCALTOK)
      vclocal(2) = (-648.9156d0)*(CON_KCALTOK)
      vclocal(3) = (757.0681d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     n2z - c4 [jcc2003]
      vnlocal(1,1) = 'n2z'
      vnlocal(1,2) = 'c4'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_COMPASS)
      vclocal(0) = 1.0221d0
      vclocal(1) = (440.1623d0)*(CON_KCALTOK)
      vclocal(2) = (-960.3246d0)*(CON_KCALTOK)
      vclocal(3) = (1120.3787d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     n2z - h1 [jcc2003]
      vnlocal(1,1) = 'n2z'
      vnlocal(1,2) = 'h1'
      call twh_set_vibnames(type,nvslocal,vnlocal)


      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Compass bond angle potentials - quartic angle                  *
c     ******************************************************************
c     --- H. Sun J. Phys. Chem. B (1998) 102, 7338-7364.
c     --- see ffnonbond for full references
c     --- harmonic, triple and quartic terms
c     --- energy = br1*(theta-theta0)**2  + br2*(theta-theta0)**3 
c     --- + br3*(theta-theta0)**4 
c     --- + br4*(bond(ij)-bond0(ij))*(theta-theta0) 
c     --- + br5*(bond(jk)-bond0(jk))*(theta-theta0)
c     --- also need to set flag for bond-bond interactions
c     --- energybondbond = br6*(bond(ij)-bondo(ij))*(bond(jk)-bond0(jk))

c     --- c3a - c3a - c3a
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 118.9000d0
      bclocal(1) = (61.0226d0)*(CON_KCALTOK)
      bclocal(2) = (-34.9931d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (28.8708d0)*(CON_KCALTOK)
      bclocal(5) = (28.8708d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (68.2856d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'c3a'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c3a - c4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 120.0500d0
      bclocal(1) = (44.7148d0)*(CON_KCALTOK)
      bclocal(2) = (-22.7352d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (31.0771d0)*(CON_KCALTOK)
      bclocal(5) = (47.0579d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (12.0676d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c3a - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 117.9400d0
      bclocal(1) = (35.1558d0)*(CON_KCALTOK)
      bclocal(2) = (-12.4682d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (20.0033d0)*(CON_KCALTOK)
      bclocal(5) = (24.2183d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (1.0795d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c4 - c3a
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 111.0000d0
      bclocal(1) = (44.3234d0)*(CON_KCALTOK)
      bclocal(2) = (-9.4454d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'c3a'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c4 - c4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 108.4000d0
      bclocal(1) = (43.9594d0)*(CON_KCALTOK)
      bclocal(2) = (-8.3924d0)*(CON_KCALTOK)
      bclocal(3) = (-9.3379d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c4 - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 111.0000d0
      bclocal(1) = (44.3234d0)*(CON_KCALTOK)
      bclocal(2) = (-9.4454d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (26.4608d0)*(CON_KCALTOK)
      bclocal(5) = (11.7717d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (2.9168d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - c4 - c4 (generic sp3 C-C-C)
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 112.670d0
      bclocal(1) = (39.5160d0)*(CON_KCALTOK)
      bclocal(2) = (-7.4430d0)*(CON_KCALTOK)
      bclocal(3) = (-9.5583d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (8.0160d0)*(CON_KCALTOK)
      bclocal(5) = (8.0160d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - c4 - h1 (generic sp3 C-C-H)
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.7700d0
      bclocal(1) = (41.4530d0)*(CON_KCALTOK)
      bclocal(2) = (-10.6040d0)*(CON_KCALTOK)
      bclocal(3) = (5.1290d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (20.7540d0)*(CON_KCALTOK)
      bclocal(5) = (11.4210d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (3.3872d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c4 - h1 (generic sp3 H-C-H)
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 107.6600d0
      bclocal(1) = (39.6410d0)*(CON_KCALTOK)
      bclocal(2) = (-12.9210d0)*(CON_KCALTOK)
      bclocal(3) = (-2.4318d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (18.1030d0)*(CON_KCALTOK)
      bclocal(5) = (18.1030d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (5.3316d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c3a - n3o
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 118.8000d0
      bclocal(1) = (29.2436d0)*(CON_KCALTOK)
      bclocal(2) = (-8.8495d0)*(CON_KCALTOK)
      bclocal(3) = (-6.6020d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (30.5211d0)*(CON_KCALTOK)
      bclocal(5) = (59.8025d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (21.0495d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'n3o'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c4 - n3o
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 107.0000d0
      bclocal(1) = (54.9318d0)*(CON_KCALTOK)
      bclocal(2) = (-9.1333d0)*(CON_KCALTOK)
      bclocal(3) = (-11.5434d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (12.2491d0)*(CON_KCALTOK)
      bclocal(5) = (30.5314d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (3.37700d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'n3o'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c4 - o2n
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 108.7280d0
      bclocal(1) = (58.5446d0)*(CON_KCALTOK)
      bclocal(2) = (-10.8088d0)*(CON_KCALTOK)
      bclocal(3) = (-12.4006d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (4.6189d0)*(CON_KCALTOK)
      bclocal(5) = (55.3270d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (23.1979d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'o2n'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - n3o - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 117.7000d0
      bclocal(1) = (63.9404d0)*(CON_KCALTOK)
      bclocal(2) = (-18.4524d0)*(CON_KCALTOK)
      bclocal(3) = (-14.3129d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (40.3757d0)*(CON_KCALTOK)
      bclocal(5) = (92.1955d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (93.7948d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'n3o'
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - n3o - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 117.5000d0
      bclocal(1) = (64.5228d0)*(CON_KCALTOK)
      bclocal(2) = (-18.4582d0)*(CON_KCALTOK)
      bclocal(3) = (-14.4215d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (27.2141d0)*(CON_KCALTOK)
      bclocal(5) = (93.9927d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (48.1403d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'n3o'
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - n3o - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 115.7000d0
      bclocal(1) = (53.8034d0)*(CON_KCALTOK)
      bclocal(2) = (-14.1991d0)*(CON_KCALTOK)
      bclocal(3) = (-11.8708d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (-8.6275d0)*(CON_KCALTOK)
      bclocal(5) = (58.6036d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (14.8266d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'n3o'
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o1= - n3o - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 128.0000d0
      bclocal(1) = (95.1035d0)*(CON_KCALTOK)
      bclocal(2) = (-47.4240d0)*(CON_KCALTOK)
      bclocal(3) = (-27.9164d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (95.6936d0)*(CON_KCALTOK)
      bclocal(5) = (95.6936d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (265.7106d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o1='
      bnlocal(1,2) = 'n3o'
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - o2n - n3o
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 108.5000d0
      bclocal(1) = (55.7454d0)*(CON_KCALTOK)
      bclocal(2) = (-10.0067d0)*(CON_KCALTOK)
      bclocal(3) = (-6.2729d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'o2n'
      bnlocal(1,3) = 'n3o'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - c4 - o2n
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 105.0000d0
      bclocal(1) = (54.5381d0)*(CON_KCALTOK)
      bclocal(2) = (-8.3642d0)*(CON_KCALTOK)
      bclocal(3) = (-13.0838d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (2.6868d0)*(CON_KCALTOK)
      bclocal(5) = (20.4033d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (11.4318d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'o2n'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o2n - n3o - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 112.8000d0
      bclocal(1) = (85.5228d0)*(CON_KCALTOK)
      bclocal(2) = (-18.4582d0)*(CON_KCALTOK)
      bclocal(3) = (-14.4215d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (80.0000d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o2n'
      bnlocal(1,2) = 'n3o'
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c3a - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 120.000d0
      bclocal(1) = (60.000d0)*(CON_KCALTOK)
      bclocal(2) = (0.000d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - c3a - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 120.001d0
      bclocal(1) = (47.841d0)*(CON_KCALTOK)
      bclocal(2) = (-15.229d0)*(CON_KCALTOK)
      bclocal(3) = (-10.907d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - c4 - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 117.317d0
      bclocal(1) = (55.242d0)*(CON_KCALTOK)
      bclocal(2) = (0.000d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (19.244d0)*(CON_KCALTOK)
      bclocal(5) = (59.422d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (22.710d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c4 - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 107.499d0
      bclocal(1) = (62.731d0)*(CON_KCALTOK)
      bclocal(2) = (0.000d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (6.407d0)*(CON_KCALTOK)
      bclocal(5) = (46.373d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (5.664d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c4 - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.886d0
      bclocal(1) = (33.830d0)*(CON_KCALTOK)
      bclocal(2) = (-7.043d0)*(CON_KCALTOK)
      bclocal(3) = (-7.246d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (19.812d0)*(CON_KCALTOK)
      bclocal(5) = (16.940d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (1.050d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - n2= - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 117.200d0
      bclocal(1) = (37.262d0)*(CON_KCALTOK)
      bclocal(2) = (0.000d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (18.486d0)*(CON_KCALTOK)
      bclocal(5) = (7.837d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (12.563d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'n2='
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - n2= - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.910d0
      bclocal(1) = (31.091d0)*(CON_KCALTOK)
      bclocal(2) = (0.000d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (8.490d0)*(CON_KCALTOK)
      bclocal(5) = (8.490d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (1.457d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'n2='
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - n2= - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 120.000d0
      bclocal(1) = (26.068d0)*(CON_KCALTOK)
      bclocal(2) = (-8.298d0)*(CON_KCALTOK)
      bclocal(3) = (-5.943d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (40.063d0)*(CON_KCALTOK)
      bclocal(5) = (90.791d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (-18.287d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'n2='
      bnlocal(1,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- p4= - n2= - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 135.000d0
      bclocal(1) = (23.868d0)*(CON_KCALTOK)
      bclocal(2) = (-8.736d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (20.000d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'p4='
      bnlocal(1,2) = 'n2='
      bnlocal(1,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - n3 - p4=
c     --- h1 - n3 - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 120.083d0
      bclocal(1) = (25.001d0)*(CON_KCALTOK)
      bclocal(2) = (-6.117d0)*(CON_KCALTOK)
      bclocal(3) = (-5.457d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'n3'
      bnlocal(1,3) = 'p4='
      bnlocal(2,1) = 'h1'
      bnlocal(2,2) = 'n3'
      bnlocal(2,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - o2 - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 118.283d0
      bclocal(1) = (35.001d0)*(CON_KCALTOK)
      bclocal(2) = (-10.360d0)*(CON_KCALTOK)
      bclocal(3) = (-7.870d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - o2 - p4=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 117.000d0
      bclocal(1) = (26.031d0)*(CON_KCALTOK)
      bclocal(2) = (-5.828d0)*(CON_KCALTOK)
      bclocal(3) = (-5.620d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'p4='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - p4= - c3a
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.231d0
      bclocal(1) = (56.185d0)*(CON_KCALTOK)
      bclocal(2) = (-17.316d0)*(CON_KCALTOK)
      bclocal(3) = (-12.728d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'c3a'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - p4= - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 108.231d0
      bclocal(1) = (36.185d0)*(CON_KCALTOK)
      bclocal(2) = (-6.488d0)*(CON_KCALTOK)
      bclocal(3) = (-7.6460)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - p4= - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 109.600d0
      bclocal(1) = (63.062d0)*(CON_KCALTOK)
      bclocal(2) = (-19.740d0)*(CON_KCALTOK)
      bclocal(3) = (-14.329d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - p4= - n3
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 108.165d0
      bclocal(1) = (70.977d0)*(CON_KCALTOK)
      bclocal(2) = (-11.548d0)*(CON_KCALTOK)
      bclocal(3) = (-15.109d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - p4= - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 107.365d0
      bclocal(1) = (71.977d0)*(CON_KCALTOK)
      bclocal(2) = (-10.943d0)*(CON_KCALTOK)
      bclocal(3) = (-15.290d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - p4= - c4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 102.500d0
      bclocal(1) = (48.232d0)*(CON_KCALTOK)
      bclocal(2) = (-5.798d0)*(CON_KCALTOK)
      bclocal(3) = (-9.966d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (12.805d0)*(CON_KCALTOK)
      bclocal(5) = (12.805d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (6.246d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - p4= - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 102.900d0
      bclocal(1) = (52.071d0)*(CON_KCALTOK)
      bclocal(2) = (-6.468d0)*(CON_KCALTOK)
      bclocal(3) = (-10.773d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (11.126d0)*(CON_KCALTOK)
      bclocal(5) = (-19.470d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (3.882d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - p4= - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 119.300d0
      bclocal(1) = (47.366d0)*(CON_KCALTOK)
      bclocal(2) = (-14.641d0)*(CON_KCALTOK)
      bclocal(3) = (-10.736d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (-7.128d0)*(CON_KCALTOK)
      bclocal(5) = (26.353d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (1.072d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - p4= - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 101.408d0
      bclocal(1) = (39.695d0)*(CON_KCALTOK)
      bclocal(2) = (-5.134d0)*(CON_KCALTOK)
      bclocal(3) = (-8.227d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (20.000d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - p4= - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.033d0
      bclocal(1) = (45.978d0)*(CON_KCALTOK)
      bclocal(2) = (-14.052d0)*(CON_KCALTOK)
      bclocal(3) = (-10.399d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (-24.383d0)*(CON_KCALTOK)
      bclocal(5) = (72.925d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (12.570d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - p4= - n3
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 103.978d0
      bclocal(1) = (68.257d0)*(CON_KCALTOK)
      bclocal(2) = (-9.221d0)*(CON_KCALTOK)
      bclocal(3) = (-14.174d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - p4= - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  103.978d0
      bclocal(1) = (73.257d0)*(CON_KCALTOK)
      bclocal(2) = (-9.897d0)*(CON_KCALTOK)
      bclocal(3) = (-15.212d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- n2= - p4= - n2=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  125.000d0
      bclocal(1) = (90.523d0)*(CON_KCALTOK)
      bclocal(2) = (-20.801d0)*(CON_KCALTOK)
      bclocal(3) = (-19.602d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (20.000d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'n2='
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n2='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- n2= - p4= - n3
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  123.215d0
      bclocal(1) = (89.923d0)*(CON_KCALTOK)
      bclocal(2) = (-32.612d0)*(CON_KCALTOK)
      bclocal(3) = (-21.096d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'n2='
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- n2= - p4= - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  112.215d0
      bclocal(1) = (99.923d0)*(CON_KCALTOK)
      bclocal(2) = (-32.093d0)*(CON_KCALTOK)
      bclocal(3) = (-22.821d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'n2='
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- n3 - p4= - n3
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  107.100d0
      bclocal(1) = (85.769d0)*(CON_KCALTOK)
      bclocal(2) = (-5.779d0)*(CON_KCALTOK)
      bclocal(3) = (-17.489d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'n3'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'n3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- n3 - p4= - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  108.300d0
      bclocal(1) = (86.769d0)*(CON_KCALTOK)
      bclocal(2) = (-5.175d0)*(CON_KCALTOK)
      bclocal(3) = (-17.671d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'n3'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o2 - p4= - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  107.500d0
      bclocal(1) = (87.769d0)*(CON_KCALTOK)
      bclocal(2) = (-4.570d0)*(CON_KCALTOK)
      bclocal(3) = (-17.852d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o2'
      bnlocal(1,2) = 'p4='
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o1= - c2= - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  180.0000d0
      bclocal(1) = (57.1000d0)*(CON_KCALTOK)
      bclocal(2) = (0.0000d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (275.4350d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o1='
      bnlocal(1,2) = 'c2='
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- s1= - c2= - s1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  180.0000d0
      bclocal(1) = (48.0000d0)*(CON_KCALTOK)
      bclocal(2) = (0.0000d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (100.7369d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 's1='
      bnlocal(1,2) = 'c2='
      bnlocal(1,3) = 's1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o1= - n2o - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  134.1000d0
      bclocal(1) = (150.0000d0)*(CON_KCALTOK)
      bclocal(2) = (-82.1013d0)*(CON_KCALTOK)
      bclocal(3) = (-40.0005d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (-50.0000d0)*(CON_KCALTOK)
      bclocal(5) = (-50.0000d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (20.0000d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o1='
      bnlocal(1,2) = 'n2o'
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o1= - s2= - o1=
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  119.3000d0
      bclocal(1) = (115.2627d0)*(CON_KCALTOK)
      bclocal(2) = (-35.6278d0)*(CON_KCALTOK)
      bclocal(3) = (-26.1261d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (45.0585d0)*(CON_KCALTOK)
      bclocal(5) = (45.0585d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (20.0000d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o1='
      bnlocal(1,2) = 's2='
      bnlocal(1,3) = 'o1='
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- Parameters from Polymer International (1997) 44, 311-330.

c     --- c3a - c3a - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  123.4200d0
      bclocal(1) = (73.6781d0)*(CON_KCALTOK)
      bclocal(2) = (-21.6787d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (58.4790d0)*(CON_KCALTOK)
      bclocal(5) = (107.6806d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (48.4754d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - c4 - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  111.2700d0
      bclocal(1) = (54.5381d0)*(CON_KCALTOK)
      bclocal(2) = (-8.3642d0)*(CON_KCALTOK)
      bclocal(3) = (-13.0838d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (2.6868d0)*(CON_KCALTOK)
      bclocal(5) = (20.4033d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (11.4318d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c4 - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  108.7280d0
      bclocal(1) = (58.5446d0)*(CON_KCALTOK)
      bclocal(2) = (-10.8088d0)*(CON_KCALTOK)
      bclocal(3) = (-12.4006d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (4.6189d0)*(CON_KCALTOK)
      bclocal(5) = (55.3270d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (23.1979d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - o2 - c4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  102.9695d0
      bclocal(1) = (38.9739d0)*(CON_KCALTOK)
      bclocal(2) = (-6.2595d0)*(CON_KCALTOK)
      bclocal(3) = (-8.1710d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.false.)
      bclocal(6) = (0.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - o2 - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  108.1900d0
      bclocal(1) = (53.1250d0)*(CON_KCALTOK)
      bclocal(2) = (-8.5016d0)*(CON_KCALTOK)
      bclocal(3) = (0.0000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (53.8614d0)*(CON_KCALTOK)
      bclocal(5) = (23.9224d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (20.6577d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - o2 - c4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  104.5000d0
      bclocal(1) = (35.7454d0)*(CON_KCALTOK)
      bclocal(2) = (-10.0067d0)*(CON_KCALTOK)
      bclocal(3) = (-6.2729d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (-2.8112d0)*(CON_KCALTOK)
      bclocal(5) = (-2.8112d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (-7.1131d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - o2 - h1
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  105.8000d0
      bclocal(1) = (52.7061d0)*(CON_KCALTOK)
      bclocal(2) = (-12.1090d0)*(CON_KCALTOK)
      bclocal(3) = (-9.8681d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (28.5800d0)*(CON_KCALTOK)
      bclocal(5) = (18.9277d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (-9.6879d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - c3a - o2
c     --- only the bond-bond cross term was listed for this angle
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  0.0d0
      bclocal(1) = (0.0d0)*(CON_KCALTOK)
      bclocal(2) = (0.0d0)*(CON_KCALTOK)
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (4.5800d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'c3a'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o2 - c4 - o2
c     --- only the bond-bond cross term was listed for this angle
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  0.0d0
      bclocal(1) = (0.0d0)*(CON_KCALTOK)
      bclocal(2) = (0.0d0)*(CON_KCALTOK)
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      bclocal(4) = (0.0d0)*(CON_KCALTOK)
      bclocal(5) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (8.2983d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o2'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- parameters from Sun; Rigby; Spectrochimica Acta Part A
c     --- 53 (1997) 1301-1323

c     --- h1 - o2 - si4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  122.8d0
      bclocal(1) = (23.7764d0)*(CON_KCALTOK)
      bclocal(2) = (-19.8152d0)*(CON_KCALTOK)
      bclocal(3) = (9.6331d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (18.0902d0)*(CON_KCALTOK)
      bclocal(5) = (31.0726d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (6.3820d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'si4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- si4 - o2 - si4
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  159.0d0
      bclocal(1) = (8.5000d0)*(CON_KCALTOK)
      bclocal(2) = (-13.4188d0)*(CON_KCALTOK)
      bclocal(3) = (-4.1785d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (28.6686d0)*(CON_KCALTOK)
      bclocal(5) = (28.6686d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (41.1143d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'si4'
      bnlocal(1,2) = 'o2'
      bnlocal(1,3) = 'si4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c3a - si4 - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  114.9d0
      bclocal(1) = (23.0218d0)*(CON_KCALTOK)
      bclocal(2) = (-31.3993d0)*(CON_KCALTOK)
      bclocal(3) = (24.9814d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.false.)
      call twh_lbonbon(GLB_SET,type,.false.)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c3a'
      bnlocal(1,2) = 'si4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- c4 - si4 - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  114.9d0
      bclocal(1) = (23.0218d0)*(CON_KCALTOK)
      bclocal(2) = (-31.3993d0)*(CON_KCALTOK)
      bclocal(3) = (24.9814d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (6.4278d0)*(CON_KCALTOK)
      bclocal(5) = (20.5669d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (5.4896d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'c4'
      bnlocal(1,2) = 'si4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- h1 - si4 - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  107.4d0
      bclocal(1) = (57.6643d0)*(CON_KCALTOK)
      bclocal(2) = (-10.6506d0)*(CON_KCALTOK)
      bclocal(3) = (4.6274d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (6.4278d0)*(CON_KCALTOK)
      bclocal(5) = (20.5669d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (11.6183d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'h1'
      bnlocal(1,2) = 'si4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- o2 - si4 - o2
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) =  110.7d0
      bclocal(1) = (70.3069d0)*(CON_KCALTOK)
      bclocal(2) = (-6.9375d0)*(CON_KCALTOK)
      bclocal(3) = (0.000d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (23.4380d0)*(CON_KCALTOK)
      bclocal(5) = (23.4380d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (41.1143d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'o2'
      bnlocal(1,2) = 'si4'
      bnlocal(1,3) = 'o2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- angles from JCC 2003

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 171.6211d0
      bclocal(1) = (47.7899d0)*(CON_KCALTOK)
      bclocal(2) = (0.0d0)*(CON_KCALTOK)
      bclocal(3) = (0.0d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (25.5611d0)*(CON_KCALTOK)
      bclocal(5) = (1.2222d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (204.9909d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     n2z - n2t - n1t [jcc 2003]
      bnlocal(1,1) = 'n2z'
      bnlocal(1,2) = 'n2t'
      bnlocal(1,3) = 'n1t'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.0345d0
      bclocal(1) = (55.7635d0)*(CON_KCALTOK)
      bclocal(2) = (0.6618d0)*(CON_KCALTOK)
      bclocal(3) = (0.0022d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (37.4419d0)*(CON_KCALTOK)
      bclocal(5) = (141.1218d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (14.9026d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     n2t - n2z - h1 [jcc 2003]
      bnlocal(1,1) = 'n2t'
      bnlocal(1,2) = 'n2z'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 113.5017d0
      bclocal(1) = (82.6294d0)*(CON_KCALTOK)
      bclocal(2) = (0.9845d0)*(CON_KCALTOK)
      bclocal(3) = (0.0033d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (195.9722d0)*(CON_KCALTOK)
      bclocal(5) = (88.2679d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (84.2075d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     n2t - n2z - c4 [jcc 2003]
      bnlocal(1,1) = 'n2t'
      bnlocal(1,2) = 'n2z'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 107.9744d0
      bclocal(1) = (52.7803d0)*(CON_KCALTOK)
      bclocal(2) = (0.6615d0)*(CON_KCALTOK)
      bclocal(3) = (0.0023d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (61.9652d0)*(CON_KCALTOK)
      bclocal(5) = (3.3182d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (18.4621d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     n2z - c4 - h1 [jcc 2003]
      bnlocal(1,1) = 'n2z'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'h1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_COMPASS_AUTO)
      bclocal(0) = 110.9900d0
      bclocal(1) = (77.9387d0)*(CON_KCALTOK)
      bclocal(2) = (0.9499d0)*(CON_KCALTOK)
      bclocal(3) = (0.0033d0)*(CON_KCALTOK)
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(4) = (67.8888d0)*(CON_KCALTOK)
      bclocal(5) = (34.8803d0)*(CON_KCALTOK)
      call twh_lbonbon(GLB_SET,type,.true.)
      bclocal(6) = (36.9309d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     n2z - c4 - c4 [jcc 2003]
      bnlocal(1,1) = 'n2z'
      bnlocal(1,2) = 'c4'
      bnlocal(1,3) = 'c4'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return


c     ******************************************************************
c     * Compass torsions with cross terms                              *
c     ******************************************************************
c     --- H. Sun J. Phys. Chem. B (1998) 102, 7338-7364.
c     --- see ffnonbond for full references
c     --- cosine series plus cross terms with bond and angle
c     --- For bonding pattern I -- J -- K -- L
c     --- k0*(1 - Cos[phi]) + k1*(1 - Cos[2*phi]) + k2*(1 - Cos[3*phi])
c         + (b(ij) - b0)*(k3*Cos[phi] + k4*Cos[2*phi] + k5*Cos[3*phi])
c         + (b(jk) - b0)*(k6*Cos[phi] + k7*Cos[2*phi] + k8*Cos[3*phi])
c         + (b(kl) - b0)*(k9*Cos[phi] + k10*Cos[2*phi] + k11*Cos[3*phi])
c         + (theta(ijk) - theta0)*( k12*Cos[phi] + k13*Cos[2*phi] 
c                                 + k14*Cos[3*phi])
c         + (theta(jkl) - theta0)*( k15*Cos[phi] + k16*Cos[2*phi] 
c                                 + k17*Cos[3*phi])
c         + k18*(theta(ijk) - theta0)*(theta(jkl) - theta0)*Cos[phi]
c         + k19*(b(ij) - b0)*(b(kl) - b0)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'c3a')
c     --- regular torsion part
      tclocal(0) = (8.3667d0)*(CON_KCALTOK)
      tclocal(1) = (1.2000d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.1185d0)*(CON_KCALTOK)
      tclocal(4) = (6.3204d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (27.5989d0)*(CON_KCALTOK)
      tclocal(7) = (-2.3120d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (-0.1185d0)*(CON_KCALTOK)
      tclocal(10) = (6.3204d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (1.9767d0)*(CON_KCALTOK)
      tclocal(13) = (1.0239d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (1.9767d0)*(CON_KCALTOK)
      tclocal(16) = (1.0239d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (0.0d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (53.0000d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c3a - c3a
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'c3a'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (4.4072d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (-0.6918d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (9.1792d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (0.2421d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (3.8987d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (-4.4683d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-14.4097d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (2.5085d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c3a - c4
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (3.9661d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (-6.8958d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (-1.1521d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (-0.4669d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (2.5014d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (2.7147d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-4.8141d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (-6.2741d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c3a - h1
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (1.5590d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (-1.7970d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (3.9421d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (-0.4879d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (-0.1242d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (3.4601d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (4.4444d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.8743d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c3a - c3a - h1
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (2.3500d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (-0.6890d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (4.8228d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (-0.6890d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (2.4501d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (2.4501d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (0.3598d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (-1.7077d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c3a - c3a - h1
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (-0.2802d0)*(CON_KCALTOK)
      tclocal(1) = (-0.0678d0)*(CON_KCALTOK)
      tclocal(2) = (-0.0122d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- c3a - c3a - c4 - c3a
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'c3a'
c     --- c3a - c3a - c4 - c4
      tnlocal(2,1) = 'c3a'
      tnlocal(2,2) = 'c3a'
      tnlocal(2,3) = 'c4'
      tnlocal(2,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (-0.2802d0)*(CON_KCALTOK)
      tclocal(1) = (-0.0678d0)*(CON_KCALTOK)
      tclocal(2) = (-0.0122d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.5835d0)*(CON_KCALTOK)
      tclocal(4) = (1.1220d0)*(CON_KCALTOK)
      tclocal(5) = (0.3978d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-5.5679d0)*(CON_KCALTOK)
      tclocal(7) = (1.4083d0)*(CON_KCALTOK)
      tclocal(8) = (0.3010d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (1.3997d0)*(CON_KCALTOK)
      tclocal(10) = (0.7756d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.2251d0)*(CON_KCALTOK)
      tclocal(13) = (0.6548d0)*(CON_KCALTOK)
      tclocal(14) = (0.1237d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (4.6266d0)*(CON_KCALTOK)
      tclocal(16) = (0.1632d0)*(CON_KCALTOK)
      tclocal(17) = (0.0461d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-5.8888d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (-3.4826d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c4 - h1
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (-0.0228d0)*(CON_KCALTOK)
      tclocal(1) = (0.0280d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1863d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c4 - c4 - h1
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0514d0)*(CON_KCALTOK)
      tclocal(2) = (-0.143d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.0732d0)*(CON_KCALTOK)
      tclocal(4) = (0.0000d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-17.7870d0)*(CON_KCALTOK)
      tclocal(7) = (-7.1877d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (-0.0732d0)*(CON_KCALTOK)
      tclocal(10) = (0.0000d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.3886d0)*(CON_KCALTOK)
      tclocal(13) = (-0.3139d0)*(CON_KCALTOK)
      tclocal(14) = (0.1389d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.3886d0)*(CON_KCALTOK)
      tclocal(16) = (-0.3139d0)*(CON_KCALTOK)
      tclocal(17) = (0.1389d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-22.0450)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - c4 - c4  sp3 carbon backbone
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0316d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1681d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.2486d0)*(CON_KCALTOK)
      tclocal(4) = (0.2422d0)*(CON_KCALTOK)
      tclocal(5) = (-0.0925d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-14.8790d0)*(CON_KCALTOK)
      tclocal(7) = (-3.6581d0)*(CON_KCALTOK)
      tclocal(8) = (-0.3138d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0814d0)*(CON_KCALTOK)
      tclocal(10) = (0.0591d0)*(CON_KCALTOK)
      tclocal(11) = (0.2219d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-0.2454d0)*(CON_KCALTOK)
      tclocal(13) = (0.0000d0)*(CON_KCALTOK)
      tclocal(14) = (-0.1136d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.3113d0)*(CON_KCALTOK)
      tclocal(16) = (0.4516d0)*(CON_KCALTOK)
      tclocal(17) = (-0.1988d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-16.1640d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - c4 - h1  sp3 carbon backbone
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (-0.1432d0)*(CON_KCALTOK)
      tclocal(1) = (0.0617d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1530d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.2130d0)*(CON_KCALTOK)
      tclocal(4) = (0.3120d0)*(CON_KCALTOK)
      tclocal(5) = (0.0777d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-14.2610d0)*(CON_KCALTOK)
      tclocal(7) = (-0.5322d0)*(CON_KCALTOK)
      tclocal(8) = (-0.4864d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.2130d0)*(CON_KCALTOK)
      tclocal(10) = (0.3120d0)*(CON_KCALTOK)
      tclocal(11) = (0.0777d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-0.8085d0)*(CON_KCALTOK)
      tclocal(13) = (0.5569d0)*(CON_KCALTOK)
      tclocal(14) = (-0.2466d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (-0.8085d0)*(CON_KCALTOK)
      tclocal(16) = (0.5569d0)*(CON_KCALTOK)
      tclocal(17) = (-0.2466d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-12.5640d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - c4 - h1  sp3 carbon backbone
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (4.5000d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- * - c3a - c3a -*  (only when an explicit value of* does not 
c     --- exist)
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1530d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- * - c4 - c4 -*  (only when an explicit value of* does not 
c     --- exist)
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'n3o')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (7.2124d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (0.0000d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (0.0000d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (7.7594d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (0.0000d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-34.9681d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c3a - n3o
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'n3o'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'n3o')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (2.9126d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (0.0000d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (0.0000d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (-8.0369d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (0.0000d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (2.1508d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c3a - c3a - n3o
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'n3o'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'n3o')
      call twh_torstring(GLB_SET,type,4,'o1=')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (1.1600d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (0.0000d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (0.0000d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (0.0000d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (-3.4207d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-18.0436d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - n3o - o1=
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'n3o'
      tnlocal(1,4) = 'o1='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.3500d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - n3o - o1=
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'n3o'
      tnlocal(1,4) = 'o1='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'n3o')
      call twh_torstring(GLB_SET,type,4,'o1=')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.3500d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (0.0000d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (0.0000d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (-0.3086d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (1.0352d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-16.2615d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - n3o - o1=
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'n3o'
      tnlocal(1,4) = 'o1='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (-0.4000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - o2 - n3o
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'n3o'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'o1=')
      call twh_torstring(GLB_SET,type,2,'n3o')
      call twh_torstring(GLB_SET,type,3,'o2')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (0.0000d0)*(CON_KCALTOK)
      tclocal(1) = (2.0000d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (0.0000d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (0.0000d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-3.0000d0)*(CON_KCALTOK)
      tclocal(13) = (0.0000d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (0.0000d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (0.0d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- o1= - n3o - o2 - c4
      tnlocal(1,1) = 'o1='
      tnlocal(1,2) = 'n3o'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (5.477d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c3a - p4=
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'p4='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (2.270d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c3a - c3a - p4=
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'p4='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (-0.272d0)*(CON_KCALTOK)
      tclocal(1) = (1.190d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- c3a - c3a - p4= - h1
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'h1'
c     --- c3a - c3a - p4= - n2=
      tnlocal(2,1) = 'c3a'
      tnlocal(2,2) = 'c3a'
      tnlocal(2,3) = 'p4='
      tnlocal(2,4) = 'n2='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'n2=')
c     --- regular torsion part
      tclocal(0) = (0.097d0)*(CON_KCALTOK)
      tclocal(1) = (0.072d0)*(CON_KCALTOK)
      tclocal(2) = (-0.258d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-2.098d0)*(CON_KCALTOK)
      tclocal(13) = (1.861d0)*(CON_KCALTOK)
      tclocal(14) = (-1.689d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (-0.122d0)*(CON_KCALTOK)
      tclocal(16) = (1.893d0)*(CON_KCALTOK)
      tclocal(17) = (-0.567d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (0.0d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - c4 - n2=
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'n2='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'n2=')
c     --- regular torsion part
      tclocal(0) = (-0.151d0)*(CON_KCALTOK)
      tclocal(1) = (0.010d0)*(CON_KCALTOK)
      tclocal(2) = (-0.186d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-3.515d0)*(CON_KCALTOK)
      tclocal(7) = (-2.298d0)*(CON_KCALTOK)
      tclocal(8) = (-1.277d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-2.098d0)*(CON_KCALTOK)
      tclocal(13) = (1.861d0)*(CON_KCALTOK)
      tclocal(14) = (-1.689d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (-0.122d0)*(CON_KCALTOK)
      tclocal(16) = (1.893d0)*(CON_KCALTOK)
      tclocal(17) = (-0.567d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-27.506d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - c4 - n2=
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'n2='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'n2=')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (-5.072d0)*(CON_KCALTOK)
      tclocal(1) = (-0.498d0)*(CON_KCALTOK)
      tclocal(2) = (-0.438d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-2.380d0)*(CON_KCALTOK)
      tclocal(7) = (2.529d0)*(CON_KCALTOK)
      tclocal(8) = (-0.730d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-2.523d0)*(CON_KCALTOK)
      tclocal(13) = (2.848d0)*(CON_KCALTOK)
      tclocal(14) = (2.059d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (-3.692d0)*(CON_KCALTOK)
      tclocal(16) = (4.061d0)*(CON_KCALTOK)
      tclocal(17) = (-1.544d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-8.898d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - n2= - h1
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'n2='
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'n2=')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (1.266d0)*(CON_KCALTOK)
      tclocal(1) = (-0.774d0)*(CON_KCALTOK)
      tclocal(2) = (0.038d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-0.414d0)*(CON_KCALTOK)
      tclocal(7) = (-2.862d0)*(CON_KCALTOK)
      tclocal(8) = (0.007d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-1.895d0)*(CON_KCALTOK)
      tclocal(13) = (1.221d0)*(CON_KCALTOK)
      tclocal(14) = (-0.746d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.110d0)*(CON_KCALTOK)
      tclocal(16) = (0.065d0)*(CON_KCALTOK)
      tclocal(17) = (0.109d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-9.628d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - n2= - h1
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'n2='
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'p4=')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.068d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0d0)*(CON_KCALTOK)
      tclocal(7) = (0.0d0)*(CON_KCALTOK)
      tclocal(8) = (0.0d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0d0)*(CON_KCALTOK)
      tclocal(13) = (0.0d0)*(CON_KCALTOK)
      tclocal(14) = (0.0d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0d0)*(CON_KCALTOK)
      tclocal(16) = (0.0d0)*(CON_KCALTOK)
      tclocal(17) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-25.546d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - p4= - c4
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'p4=')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.068d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0d0)*(CON_KCALTOK)
      tclocal(7) = (0.0d0)*(CON_KCALTOK)
      tclocal(8) = (0.0d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0d0)*(CON_KCALTOK)
      tclocal(13) = (0.0d0)*(CON_KCALTOK)
      tclocal(14) = (0.0d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0d0)*(CON_KCALTOK)
      tclocal(16) = (0.0d0)*(CON_KCALTOK)
      tclocal(17) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-16.018d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - p4= - h1
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'p4=')
      call twh_torstring(GLB_SET,type,4,'n2=')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.068d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0d0)*(CON_KCALTOK)
      tclocal(7) = (0.0d0)*(CON_KCALTOK)
      tclocal(8) = (0.0d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0d0)*(CON_KCALTOK)
      tclocal(13) = (0.0d0)*(CON_KCALTOK)
      tclocal(14) = (0.0d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0d0)*(CON_KCALTOK)
      tclocal(16) = (0.0d0)*(CON_KCALTOK)
      tclocal(17) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-19.934d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - p4= - n2=
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'n2='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.350d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - n2= - p4= - c3a
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'c3a'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'n2=')
      call twh_torstring(GLB_SET,type,3,'p4=')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.369d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0d0)*(CON_KCALTOK)
      tclocal(7) = (0.0d0)*(CON_KCALTOK)
      tclocal(8) = (0.0d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0d0)*(CON_KCALTOK)
      tclocal(13) = (0.0d0)*(CON_KCALTOK)
      tclocal(14) = (0.0d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0d0)*(CON_KCALTOK)
      tclocal(16) = (0.0d0)*(CON_KCALTOK)
      tclocal(17) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-11.102d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - n2= - p4= - c4
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'n2=')
      call twh_torstring(GLB_SET,type,3,'p4=')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.350d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0d0)*(CON_KCALTOK)
      tclocal(4) = (0.0d0)*(CON_KCALTOK)
      tclocal(5) = (0.0d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0d0)*(CON_KCALTOK)
      tclocal(7) = (0.0d0)*(CON_KCALTOK)
      tclocal(8) = (0.0d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0d0)*(CON_KCALTOK)
      tclocal(10) = (0.0d0)*(CON_KCALTOK)
      tclocal(11) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0d0)*(CON_KCALTOK)
      tclocal(13) = (0.0d0)*(CON_KCALTOK)
      tclocal(14) = (0.0d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0d0)*(CON_KCALTOK)
      tclocal(16) = (0.0d0)*(CON_KCALTOK)
      tclocal(17) = (0.0d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-3.788d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - n2= - p4= - h1
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.350d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - n2= - p4= - o2
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
c     --- null torsion
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 5
c     --- p4= - n2= - p4= - h1
      tnlocal(1,1) = 'p4='
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'h1'
c     --- h1 - n3 - p4= - h1
      tnlocal(2,1) = 'h1'
      tnlocal(2,2) = 'n3'
      tnlocal(2,3) = 'p4='
      tnlocal(2,4) = 'h1'
c     --- h1 - n3 - p4= - n2=
      tnlocal(3,1) = 'h1'
      tnlocal(3,2) = 'n3'
      tnlocal(3,3) = 'p4='
      tnlocal(3,4) = 'n2='
c     --- h1 - n3 - p4= - o2
      tnlocal(4,1) = 'h1'
      tnlocal(4,2) = 'n3'
      tnlocal(4,3) = 'p4='
      tnlocal(4,4) = 'o2'
c     ---* - n3 - p4= -*
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'n3'
      tnlocal(5,3) = 'p4='
      tnlocal(5,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (1.800d0)*(CON_KCALTOK)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = (2.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- p4= - n2= - p4= - n2=
      tnlocal(1,1) = 'p4='
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'n2='
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (5.708d0)*(CON_KCALTOK)
      tclocal(1) = (2.118d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- h1 - o2 - p4= - h1
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'o2'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'h1'
c     --- h1 - o2 - p4= - n2=
      tnlocal(2,1) = 'h1'
      tnlocal(2,2) = 'o2'
      tnlocal(2,3) = 'p4='
      tnlocal(2,4) = 'n2='
c     --- h1 - o2 - p4= - o2
      tnlocal(3,1) = 'h1'
      tnlocal(3,2) = 'o2'
      tnlocal(3,3) = 'p4='
      tnlocal(3,4) = 'o2'
c     ---* - o2 - p4= -* (check to see if explicitly declared elsewhere)
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'o2'
      tnlocal(4,3) = 'p4='
      tnlocal(4,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (1.000d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ---* - c3a - n2= -*
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'n2='
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (-0.272d0)*(CON_KCALTOK)
      tclocal(1) = (1.190d0)*(CON_KCALTOK)
      tclocal(2) = (0.000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ---* - c3a - p4= -*
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.020d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ---* - c4 - n2= -*
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'n2='
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.068d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ---* - c4 - p4= -*
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.350d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     ---* - n2= - p4= -*
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'n2='
      tnlocal(1,3) = 'p4='
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'o2')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (4.8498d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (0.2655d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (4.8255d0)*(CON_KCALTOK)
      tclocal(8) = (0.0d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (4.8905d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (10.0155d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (1.7404d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-21.0247d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (-2.2436d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - c3a - o2
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'c3a')
      call twh_torstring(GLB_SET,type,4,'o2')
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (1.7234d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.0000d0)*(CON_KCALTOK)
      tclocal(4) = (-1.5867d0)*(CON_KCALTOK)
      tclocal(5) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.0000d0)*(CON_KCALTOK)
      tclocal(7) = (5.5432d0)*(CON_KCALTOK)
      tclocal(8) = (0.0000d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.0000d0)*(CON_KCALTOK)
      tclocal(10) = (4.2641d0)*(CON_KCALTOK)
      tclocal(11) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.0000d0)*(CON_KCALTOK)
      tclocal(13) = (1.8729d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.0000d0)*(CON_KCALTOK)
      tclocal(16) = (2.5706d0)*(CON_KCALTOK)
      tclocal(17) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (4.2296d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (2.0517d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c3a - c3a - o2
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'c3a'
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (1.5000d0)*(CON_KCALTOK)
      tclocal(2) = (0.0000d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - o2 - c4
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c3a')
      call twh_torstring(GLB_SET,type,2,'c3a')
      call twh_torstring(GLB_SET,type,3,'o2')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (-0.6900d0)*(CON_KCALTOK)
      tclocal(1) = (0.5097d0)*(CON_KCALTOK)
      tclocal(2) = (0.0095d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.9000d0)*(CON_KCALTOK)
      tclocal(4) = (-1.3456d0)*(CON_KCALTOK)
      tclocal(5) = (1.1900d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (1.1580d0)*(CON_KCALTOK)
      tclocal(7) = (3.2697d0)*(CON_KCALTOK)
      tclocal(8) = (3.5132d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (3.4132d0)*(CON_KCALTOK)
      tclocal(10) = (0.5873d0)*(CON_KCALTOK)
      tclocal(11) = (-0.1323d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-5.1360d0)*(CON_KCALTOK)
      tclocal(13) = (-1.0122d0)*(CON_KCALTOK)
      tclocal(14) = (0.0000d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (4.6852d0)*(CON_KCALTOK)
      tclocal(16) = (0.0230d0)*(CON_KCALTOK)
      tclocal(17) = (-0.5980d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-4.6072d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (1.1590d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c3a - c3a - o2 - h1
      tnlocal(1,1) = 'c3a'
      tnlocal(1,2) = 'c3a'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'o2')
c     --- regular torsion part
      tclocal(0) = (0.7137d0)*(CON_KCALTOK)
      tclocal(1) = (0.2660d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2545d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.3190d0)*(CON_KCALTOK)
      tclocal(4) = (0.4411d0)*(CON_KCALTOK)
      tclocal(5) = (-0.7174d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-21.8842d0)*(CON_KCALTOK)
      tclocal(7) = (-7.6764d0)*(CON_KCALTOK)
      tclocal(8) = (-0.6868d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (1.1538d0)*(CON_KCALTOK)
      tclocal(10) = (0.8409d0)*(CON_KCALTOK)
      tclocal(11) = (-0.9138d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.5623d0)*(CON_KCALTOK)
      tclocal(13) = (-0.3041d0)*(CON_KCALTOK)
      tclocal(14) = (-0.4015d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.9672d0)*(CON_KCALTOK)
      tclocal(16) = (-0.7566d0)*(CON_KCALTOK)
      tclocal(17) = (-1.2331d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-29.0420d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - c4 - o2
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
c     --- there is a typo in the original paper for the 
c     --- angle-angle-torsion interaction.  The negative sign was 
c     --- dropped by mistake.  It has been corrected here.
c     --- personal communication David Rigby 08-09-2001 M.G. Martin
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'o2')
c     --- regular torsion part
      tclocal(0) = (-0.1435d0)*(CON_KCALTOK)
      tclocal(1) = (0.2530d0)*(CON_KCALTOK)
      tclocal(2) = (-0.0905d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (0.9681d0)*(CON_KCALTOK)
      tclocal(4) = (0.9551d0)*(CON_KCALTOK)
      tclocal(5) = (0.0436d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-16.7975d0)*(CON_KCALTOK)
      tclocal(7) = (-1.2296d0)*(CON_KCALTOK)
      tclocal(8) = (-0.2750d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.5903d0)*(CON_KCALTOK)
      tclocal(10) = (0.6669d0)*(CON_KCALTOK)
      tclocal(11) = (0.8584d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (2.3668d0)*(CON_KCALTOK)
      tclocal(13) = (2.4920d0)*(CON_KCALTOK)
      tclocal(14) = (-1.0122d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (-0.1892d0)*(CON_KCALTOK)
      tclocal(16) = (0.4918d0)*(CON_KCALTOK)
      tclocal(17) = (0.7273d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-20.2006d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - c4 - o2
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'o2')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'c4')
      call twh_torstring(GLB_SET,type,4,'o2')
c     --- regular torsion part
      tclocal(0) = (1.1000d0)*(CON_KCALTOK)
      tclocal(1) = (-0.0500d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1441d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (1.0165d0)*(CON_KCALTOK)
      tclocal(4) = (0.7553d0)*(CON_KCALTOK)
      tclocal(5) = (-0.4609d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-17.2585d0)*(CON_KCALTOK)
      tclocal(7) = (-3.6157d0)*(CON_KCALTOK)
      tclocal(8) = (-0.8364d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (1.0165d0)*(CON_KCALTOK)
      tclocal(10) = (0.7553d0)*(CON_KCALTOK)
      tclocal(11) = (-0.4609d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (0.5511d0)*(CON_KCALTOK)
      tclocal(13) = (0.9737d0)*(CON_KCALTOK)
      tclocal(14) = (-0.6673d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.5511d0)*(CON_KCALTOK)
      tclocal(16) = (0.9737d0)*(CON_KCALTOK)
      tclocal(17) = (-0.6673d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-14.0484d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- o2 - c4 - c4 -o2
      tnlocal(1,1) = 'o2'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'o2')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (-0.4000d0)*(CON_KCALTOK)
      tclocal(1) = (-0.4028d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2450d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.2456d0)*(CON_KCALTOK)
      tclocal(4) = (1.0517d0)*(CON_KCALTOK)
      tclocal(5) = (-0.7795d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-5.9288d0)*(CON_KCALTOK)
      tclocal(7) = (-2.7007d0)*(CON_KCALTOK)
      tclocal(8) = (-0.3175d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.4741d0)*(CON_KCALTOK)
      tclocal(10) = (1.2635d0)*(CON_KCALTOK)
      tclocal(11) = (0.5576d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-2.7466d0)*(CON_KCALTOK)
      tclocal(13) = (1.4877d0)*(CON_KCALTOK)
      tclocal(14) = (-0.8955d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.5676d0)*(CON_KCALTOK)
      tclocal(16) = (0.9450d0)*(CON_KCALTOK)
      tclocal(17) = (0.0703d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-19.0059d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - o2 - c4
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'c4')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'o2')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (-0.6732d0)*(CON_KCALTOK)
      tclocal(1) = (-0.4778d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1670d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.5800d0)*(CON_KCALTOK)
      tclocal(4) = (0.9004d0)*(CON_KCALTOK)
      tclocal(5) = (0.000d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (1.2472d0)*(CON_KCALTOK)
      tclocal(7) = (0.0000d0)*(CON_KCALTOK)
      tclocal(8) = (0.7485d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.000d0)*(CON_KCALTOK)
      tclocal(10) = (0.5343d0)*(CON_KCALTOK)
      tclocal(11) = (0.9025d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-3.5903d0)*(CON_KCALTOK)
      tclocal(13) = (2.5225d0)*(CON_KCALTOK)
      tclocal(14) = (0.4888d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.8726d0)*(CON_KCALTOK)
      tclocal(16) = (-0.3577d0)*(CON_KCALTOK)
      tclocal(17) = (0.3888d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-12.1038d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- c4 - c4 - o2 - h1
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.9513d0)*(CON_KCALTOK)
      tclocal(1) = (0.1155d0)*(CON_KCALTOK)
      tclocal(2) = (0.0720d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - o2 - c3a
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'c3a'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'o2')
      call twh_torstring(GLB_SET,type,4,'c4')
c     --- regular torsion part
      tclocal(0) = (0.5302d0)*(CON_KCALTOK)
      tclocal(1) = (0.0000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.3966d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-0.6054d0)*(CON_KCALTOK)
      tclocal(4) = (1.3339d0)*(CON_KCALTOK)
      tclocal(5) = (0.9648d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (-6.8007d0)*(CON_KCALTOK)
      tclocal(7) = (-4.6546d0)*(CON_KCALTOK)
      tclocal(8) = (-1.4101d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (-0.1620d0)*(CON_KCALTOK)
      tclocal(10) = (0.1564d0)*(CON_KCALTOK)
      tclocal(11) = (-1.1408d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-1.8234d0)*(CON_KCALTOK)
      tclocal(13) = (1.6393d0)*(CON_KCALTOK)
      tclocal(14) = (0.5144d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (-0.7777d0)*(CON_KCALTOK)
      tclocal(16) = (0.4340d0)*(CON_KCALTOK)
      tclocal(17) = (-0.6653d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-16.4438d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - o2 - c4
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,5)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_torstring(GLB_SET,type,1,'h1')
      call twh_torstring(GLB_SET,type,2,'c4')
      call twh_torstring(GLB_SET,type,3,'o2')
      call twh_torstring(GLB_SET,type,4,'h1')
c     --- regular torsion part
      tclocal(0) = (0.1863d0)*(CON_KCALTOK)
      tclocal(1) = (-0.4338d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2121d0)*(CON_KCALTOK)
c     --- bond(ij) - torsion cross
      tclocal(3) = (-1.7554d0)*(CON_KCALTOK)
      tclocal(4) = (1.3145d0)*(CON_KCALTOK)
      tclocal(5) = (0.2263d0)*(CON_KCALTOK)
c     --- bond(jk) - torsion cross
      tclocal(6) = (0.000d0)*(CON_KCALTOK)
      tclocal(7) = (0.9241d0)*(CON_KCALTOK)
      tclocal(8) = (-0.5889d0)*(CON_KCALTOK)
c     --- bond(kl) - torsion cross
      tclocal(9) = (0.2493d0)*(CON_KCALTOK)
      tclocal(10) = (0.6803d0)*(CON_KCALTOK)
      tclocal(11) = (0.000d0)*(CON_KCALTOK)
c     --- angle(ijk) - torsion cross
      tclocal(12) = (-3.4060d0)*(CON_KCALTOK)
      tclocal(13) = (1.6396d0)*(CON_KCALTOK)
      tclocal(14) = (0.0737d0)*(CON_KCALTOK)
c     --- angle(jkl) - torsion cross
      tclocal(15) = (0.000d0)*(CON_KCALTOK)
      tclocal(16) = (-0.2810d0)*(CON_KCALTOK)
      tclocal(17) = (-0.5944d0)*(CON_KCALTOK)
c     --- angle(ijk) - angle(jkl) - torsion cross
      tclocal(18) = (-10.5093d0)*(CON_KCALTOK)
c     --- bond(ij) - bond(kl) cross
      tclocal(19) = (0.0d0)*CON_KCALTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- h1 - c4 - o2 - h1
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'o2'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.050d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- h1 - o2 - si4 - c4
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'o2'
      tnlocal(1,3) = 'si4'
      tnlocal(1,4) = 'c4'
c     --- h1 - o2 - si4 - h1
      tnlocal(2,1) = 'h1'
      tnlocal(2,2) = 'o2'
      tnlocal(2,3) = 'si4'
      tnlocal(2,4) = 'h1'
c     --- h1 - o2 - si4 - o2
      tnlocal(3,1) = 'h1'
      tnlocal(3,2) = 'o2'
      tnlocal(3,3) = 'si4'
      tnlocal(3,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.010d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- si4 - o2 - si4 - c4
      tnlocal(1,1) = 'si4'
      tnlocal(1,2) = 'o2'
      tnlocal(1,3) = 'si4'
      tnlocal(1,4) = 'c4'
c     --- si4 - o2 - si4 - h1
      tnlocal(2,1) = 'si4'
      tnlocal(2,2) = 'o2'
      tnlocal(2,3) = 'si4'
      tnlocal(2,4) = 'h1'
c     --- * - o2 - si4 - *
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'o2'
      tnlocal(3,3) = 'si4'
      tnlocal(3,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (-0.225d0)*(CON_KCALTOK)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = (-0.010d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- si4 - o2 - si4 - o2
      tnlocal(1,1) = 'si4'
      tnlocal(1,2) = 'o2'
      tnlocal(1,3) = 'si4'
      tnlocal(1,4) = 'o2'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- torsions from [JCC 2003]

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2637d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     h1 - n2z - n2t - n1t  [JCC2003]
      tnlocal(1,1) = 'h1'
      tnlocal(1,2) = 'n2z'
      tnlocal(1,3) = 'n2t'
      tnlocal(1,4) = 'n1t'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(2) = (-0.1823d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     c4 - n2z - n2t - n1t  [JCC2003]
      tnlocal(1,1) = 'c4'
      tnlocal(1,2) = 'n2z'
      tnlocal(1,3) = 'n2t'
      tnlocal(1,4) = 'n1t'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2181d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     n2t - n2z - c4 - h1  [JCC2003]
      tnlocal(1,1) = 'n2t'
      tnlocal(1,2) = 'n2z'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2021d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     n2t - n2z - c4 - c4  [JCC2003]
      tnlocal(1,1) = 'n2t'
      tnlocal(1,2) = 'n2z'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'c4'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,6)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
c     --- regular torsion part
      tclocal(0) = (0.0d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(CON_KCALTOK)
      tclocal(2) = (-0.2259d0)*(CON_KCALTOK)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     n2z - c4 - c4 - h1  [JCC2003]
      tnlocal(1,1) = 'n2z'
      tnlocal(1,2) = 'c4'
      tnlocal(1,3) = 'c4'
      tnlocal(1,4) = 'h1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Compass improper harmonic torsion                              *
c     * used exclusively to control out of plane                       *
c     ******************************************************************
c     --- listed as IJKL where J is the central atom in a planar center
c     --- torsion goes J-I-K-L 

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (7.1794d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - c3a
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'c3a'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (7.8153d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - c4
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'c4'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (4.8912d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - h1
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'h1'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.9194d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - n3o
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'n3o'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (36.2612d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - n3o - o1= - o1=
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'n3o'
      ninlocal(1,3) = 'o1='
      ninlocal(1,4) = 'o1='
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (44.3062d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c4 - n3o - o1= - o1=
      ninlocal(1,1) = 'c4'
      ninlocal(1,2) = 'n3o'
      ninlocal(1,3) = 'o1='
      ninlocal(1,4) = 'o1='
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (38.5581d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- h1 - n3o - o1= - o1=
      ninlocal(1,1) = 'h1'
      ninlocal(1,2) = 'n3o'
      ninlocal(1,3) = 'o1='
      ninlocal(1,4) = 'o1='
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (45.0000d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- o1= - n3o - o1= - o2
      ninlocal(1,1) = 'o1='
      ninlocal(1,2) = 'n3o'
      ninlocal(1,3) = 'o1='
      ninlocal(1,4) = 'o2'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (8.000d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - n2=
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'n2='
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (6.709d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - p4=
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'p4='
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (13.0421d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- c3a - c3a - c3a - o2
      ninlocal(1,1) = 'c3a'
      ninlocal(1,2) = 'c3a'
      ninlocal(1,3) = 'c3a'
      ninlocal(1,4) = 'o2'
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

c     --- H. Sun J. Phys. Chem. B (1998) 102, 7338-7364.
c     --- aacoeff*(theta(KJI) - theta0(KJI))*(theta(KJL) - theta0(KJL))
c                 L
c           K -- J
c                 I
c     --- defined from atom J to K to I to L
c     --- recorded below as I-J-K-L as in the Sun paper
c     --- in towhee_input  this is assigned from atom J with 
c     --- K,I,L in list

      type = 1
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (2.0403d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - c3a - h1
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c3a'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 2
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (3.0118d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - c3a - h1
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c3a'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 3
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-1.8202d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c3a - c4 - c4 - h1
      aanslocal = 1
      aanlocal(1,1) = 'c3a'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 4
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-0.1729d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - c4 - c4 ([C-C]-C & [C-C]-C sp3 carbon chain)
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'c4'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 5
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-1.3199d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - c4 - h1 ([C-C]-C & [C-C]-H sp3 carbon chain)
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 6
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-0.4825d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - c4 - h1 ([C-C]-H & [C-C]-H sp3 carbon chain)
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 7
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (1.0827d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c3a - c4 - h1 - c4
      aanslocal = 1
      aanlocal(1,1) = 'c3a'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'c4'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 8
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (2.3794d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c3a - c4 - h1 - h1
      aanslocal = 1
      aanlocal(1,1) = 'c3a'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 9
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (0.1184d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - h1 - c4 ([H-C]-C & [H-C]-C sp3 carbon chain)
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'c4'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 10 
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (0.2738d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - h1 - h1 ([H-C]-C & [H-C]-H sp3 carbon chain)
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 11
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-0.3157d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - h1 - h1 ([H-C]-H & [H-C]-H sp3 carbon chain)
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 12
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (1.091d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - c4 - n2=
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'n2='
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 13
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (2.753d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - h1 - n2=
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'n2='
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 14
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (1.768d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - h1 - n2=
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'n2='
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 15
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-1.306d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - n2= - h1
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'n2='
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 16
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-2.947d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - n2= - h1
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'n2='
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

c     --- parameters from Polymer International 44,(1997)

      type = 17
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-0.8330d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - c4 - o2
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'o2'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 18
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (2.5926d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - c4 - o2 [C-C]-H and [C-C]-O
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'c4'
      aanlocal(1,4) = 'o2'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 19
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (3.9177d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - h1 - o2 [H-C]-C and [H-C]-O
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'o2'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 20 
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (2.4259d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - h1 - o2 [H-C]-H and [H-C]-O
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'h1'
      aanlocal(1,4) = 'o2'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 21
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (-3.5744d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - o2 - c4
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'o2'
      aanlocal(1,4) = 'c4'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 22
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (0.1689d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- c4 - c4 - o2 - h1 [O-C]-C and [O-C]-H
      aanslocal = 1
      aanlocal(1,1) = 'c4'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'o2'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

      type = 23
      call twh_aastyle(GLB_SET,type,1)
      aaclocal(0) = (2.1283d0)*(CON_KCALTOK)
      call twh_set_aacoeff(type,aaclocal)
      call twh_aaff(GLB_SET,type,ffname)
c     --- h1 - c4 - o2 - h1 [O-C]-H and [O-C]-H
      aanslocal = 1
      aanlocal(1,1) = 'h1'
      aanlocal(1,2) = 'c4'
      aanlocal(1,3) = 'o2'
      aanlocal(1,4) = 'h1'
      call twh_set_aanames(type,aanslocal,aanlocal)

c     --- Null angle-angle groups - there is absolutely no reason to 
c     --- actually use these in a simulation, but they are listed here 
c     --- for completeness
c     --- c3a - c3a - c3a - o2
c     --- c3a - c3a - o2 - c3a
      
      naatypes = type
      call twh_ffcheck(lfailure,5,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [jpcb 1998]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'c3a'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [jpcb 1998]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'c4'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.1268d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [jpcb 1998]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'h1'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [jpcb 1998]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'c4'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.0530d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [jpcb 1998]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'h1'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.2225d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [spec Acta 1997]
      qbinlocal(1) = 'o2'
      qbinlocal(2) = 'si4'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.3110d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [JCC 2003]
      qbinlocal(1) = 'n2z'
      qbinlocal(2) = 'c4'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.3350d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [JCC 2003]
      qbinlocal(1) = 'n2z'
      qbinlocal(2) = 'h1'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.3860d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [JCC 2003]
      qbinlocal(1) = 'n2t'
      qbinlocal(2) = 'n1t'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.2470d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [JCC 2003]
      qbinlocal(1) = 'n2t'
      qbinlocal(2) = 'n2z'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0420d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Polymer International 1997]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'o2e'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0420d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Polymer International 1997]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'o2h'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.1600d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Polymer International 1997]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'o2e'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.1600d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Polymer International 1997]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'o2h'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.4200d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Polymer International 1997]
      qbinlocal(1) = 'h1'
      qbinlocal(2) = 'o2'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.4200d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Polymer International 1997]
      qbinlocal(1) = 'h1'
      qbinlocal(2) = 'o2h'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.199d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'n2='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.060d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.345d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'n2='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.050d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.120d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'cl1p'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.180d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'f1p'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.328d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'h1'
      qbinlocal(2) = 'n2='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.050d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'h1'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'n2='
      qbinlocal(2) = 'n2='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.025d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'n2='
      qbinlocal(2) = 'n3'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.043d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'n2='
      qbinlocal(2) = 'o2'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.350d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'n2='
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.120d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'n3'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.140d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [CTPS 1998]
      qbinlocal(1) = 'o2'
      qbinlocal(2) = 'p4='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.2390d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'c3a'
      qbinlocal(2) = 'n3o'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.2100d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'n3o'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.3170d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'c4'
      qbinlocal(2) = 'o2n'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.2100d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'c4o'
      qbinlocal(2) = 'n3o'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.3170d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'c4o'
      qbinlocal(2) = 'o2n'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.1880d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'h1'
      qbinlocal(2) = 'n3o'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.4280d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'n3o'
      qbinlocal(2) = 'o1='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0010d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Bunte JPCB 2000]
      qbinlocal(1) = 'n3o'
      qbinlocal(2) = 'o2n'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.0203d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'c1o'
      qbinlocal(2) = 'o1c'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.4000d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'c2='
      qbinlocal(2) = 'o1='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0258d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'c2='
      qbinlocal(2) = 's1='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0730d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'n2o'
      qbinlocal(2) = 'o1='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'h1h'
      qbinlocal(2) = 'h1h'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'n1n'
      qbinlocal(2) = 'n1n'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0288d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'n1o'
      qbinlocal(2) = 'o1n'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.2351d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'o1='
      qbinlocal(2) = 's2='
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     [Yang JPCB 2000]
      qbinlocal(1) = 'o1o'
      qbinlocal(2) = 'o1o'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return
      
      return
      end
