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
      subroutine twh_ffcharmm19(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the charmm19     *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_NBCOEFF
#define FUNCTION_GET_NBNAME
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      integer ii,type,isame
      integer itrial,i,j,jlo,jhi,jj,ij
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
      double precision onepi

      ffname = 'Charmm19'
      type = 0
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     **************************************************************
c     * Charmm19 force field                                       *
c     * contains all of the parameters listed in the paper of      *
c     * E. Neria; S. Fischer; M. Karplus "Simulation of activation *
c     * free energies in molecular systems" J. Chem. Phys. 105 (5) *
c     * 1902-1921 (1996)                                           *
c     * I am using the terminology found in that paper for the     *
c     * naming conventions                                         *
c     **************************************************************

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C: Carbonyl carbon [neria1996]
      nbnlocal ='C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.365d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0486d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH: Extended aliphatic carbon with one hydrogen -CH- group
c     [neria1996]
      nbnlocal ='CH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.235d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1142d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C(2): Extended aliphatic carbon with two hydrogens -CH2- group
c     [neria1996]
      nbnlocal ='C(2)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.165d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1811d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     C(3): Extended methyl terminal -CH3 group [neria1996]
      nbnlocal ='C(3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.100d0)*(CON_KCALTOK)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CR: Extended aromatic carbon with one hydrogen [neria1996]
      nbnlocal ='CR'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0498d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     H: Polar hydrogen [neria1996]
      nbnlocal ='H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0498d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     HC: Polar hydrogen in Arg., Lys., and N term [neria1996]
      nbnlocal ='HC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.8000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0498d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     HT: Water hydrogen, modified TIP3P model [neria1996]
      nbnlocal ='HT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     N: Nitrogen with no hydrogens [neria1996]
      nbnlocal ='N'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     N1: Nitrogen bound to one hydrogen [neria1996]
      nbnlocal ='N1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     N2: Nitrogen bound to two hydrogens [neria1996]
      nbnlocal ='N2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     N3: Nitrogen bound to three hydrogens [neria1996]
      nbnlocal ='N3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     NC2: Nitrogen bound to two hydrogens in ARG [neria1996]
      nbnlocal ='NC2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     NR: Nitrogen in an aromatic ring with no hydrogens [neria1996]
      nbnlocal ='NR'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1591d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     O: carbonyl oxygen [neria1996] 
      nbnlocal ='O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.6469d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     OC: carboxyl oxygen [neria1996]
      nbnlocal ='OC'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1591d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     OH: Hydroxyl oxygen [neria1996] 
      nbnlocal ='OH'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.6000d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1591d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (OT) water oxygen, modified TIP3P [neria1996]
      nbnlocal ='OT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.890d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0430d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     (S) sulfur
      nbnlocal ='S'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     --- setup all of the cross terms explicitly 
      do itrial = 1,2
         do i = 1, ntypes
            if ( itrial .eq. 1 ) then
c               --- we are setting up all of the dissimilar cross terms
               jlo = i+1
               jhi = ntypes
            else
c              --- we are setting up all of the self terms
               jlo = i
               jhi = i
            endif
            do j = jlo,jhi
               ii = twh_getnbtype(i,i)
               jj = twh_getnbtype(j,j)
               ij = twh_getnbtype(i,j)
c              --- only valid for the charmm19 force field (LJ)
c              --- Lorentz-Berthelot rules 
c              --- arimethic mean of sigma
               nbclocal(1)=(0.5d0*(twh_get_nbcoeff(1,ii)
     &              + twh_get_nbcoeff(1,jj)))
c              --- geometric mean of epsilon
               nbclocal(2) = dsqrt(twh_get_nbcoeff(2,ii)
     &              *twh_get_nbcoeff(2,jj) )
c              --- 1-4 nonbonded interactions
               nbclocal(3) = (0.5d0*(twh_get_nbcoeff(3,ii)
     &              + twh_get_nbcoeff(3,jj)))
               nbclocal(4) = dsqrt(twh_get_nbcoeff(4,ii)
     &              *twh_get_nbcoeff(4,jj) )
c              --- Charmm19 hacks for water
               if ( twh_get_nbname(i) .eq. 'OT'
     &              .and. twh_get_nbname(j) .eq. 'OT' ) then
c                 --- OT OT water hack
                  nbclocal(1) = 3.5365d0
                  nbclocal(2) = 0.152073d0*CON_KCALTOK
c                 --- don't worry about 1-4 for water
                  nbclocal(3) = 0.0d0
                  nbclocal(4) = 0.0d0
               elseif ( twh_get_nbname(i) .eq. 'HT'
     &                 .and. twh_get_nbname(j) .eq. 'HT' ) then
c                 --- HT HT water hack
                  nbclocal(1) = 0.4490
                  nbclocal(2) = 0.04598d0*CON_KCALTOK
c                 --- don't worry about 1-4 for water
                  nbclocal(3) = 0.0d0
                  nbclocal(4) = 0.0d0
               elseif ( (twh_get_nbname(i) .eq. 'HT' 
     &                 .and. twh_get_nbname(j) .eq. 'OT')
     &                 .or. ( twh_get_nbname(i) .eq. 'OT'
     &                 .and. twh_get_nbname(j) .eq. 'HT') ) then
c                --- OT HT water hack
                  nbclocal(1) = 1.9927d0
                  nbclocal(2) = 0.08363d0*CON_KCALTOK
c                 --- don't worry about 1-4 for water
                  nbclocal(3) = 0.0d0
                  nbclocal(4) = 0.0d0
               endif
               call twh_set_nbcoeff(ij,nbclocal)
            enddo
         enddo
      enddo

c     ******************************************************************
c     * Charmm19                                                       *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - C [neria1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.52d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- C - CH [neria1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CH'
c     --- C - C(2) [neria1996]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'C(2)'
c     --- C - C(3) [neria1996]
      vnlocal(3,1) = 'C'
      vnlocal(3,2) = 'C(3)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- C - CR [neria1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CR'
c     --- C - OH [neria1996]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'OH'
c     --- CR - CR [neria1996]
      vnlocal(3,1) = 'CR'
      vnlocal(3,2) = 'CR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.33d0
      vclocal(1) = (471.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 4
c     --- C - N [neria1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'N'
c     --- C - N1 [neria1996]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'N1'
c     --- C - N2 [neria1996]
      vnlocal(3,1) = 'C'
      vnlocal(3,2) = 'N2'
c     --- C - NR [neria1996]
      vnlocal(4,1) = 'C'
      vnlocal(4,2) = 'NR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.33d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C - NC2 [neria1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'NC2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.23d0
      vclocal(1) = (580.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C - O [neria1996]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O'
c     --- C - OC [neria1996]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'OC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.53d0
      vclocal(1) = (225.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CH - CH [neria1996]
      vnlocal(1,1) = 'CH'
      vnlocal(1,2) = 'CH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.52d0
      vclocal(1) = (225.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CH - C(2) [neria1996]
      vnlocal(1,1) = 'CH'
      vnlocal(1,2) = 'C(2)'
c     --- CH - C(3) [neria1996]
      vnlocal(2,1) = 'CH'
      vnlocal(2,2) = 'C(3)'
c     --- C(2) - C(2) [neria1996]
      vnlocal(3,1) = 'C(2)'
      vnlocal(3,2) = 'C(2)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.45d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 8
c     --- CH - N [neria1996]
      vnlocal(1,1) = 'CH'
      vnlocal(1,2) = 'N'
c     --- CH - N1 [neria1996]
      vnlocal(2,1) = 'CH'
      vnlocal(2,2) = 'N1'
c     --- CH - N2 [neria1996]
      vnlocal(3,1) = 'CH'
      vnlocal(3,2) = 'N2'
c     --- CH - N3 [neria1996]
      vnlocal(4,1) = 'CH'
      vnlocal(4,2) = 'N3'
c     --- C(2) - N [neria1996]
      vnlocal(5,1) = 'C(2)'
      vnlocal(5,2) = 'N'
c     --- C(2) - N1 [neria1996]
      vnlocal(6,1) = 'C(2)'
      vnlocal(6,2) = 'N1'
c     --- C(2) - N2 [neria1996]
      vnlocal(7,1) = 'C(2)'
      vnlocal(7,2) = 'N2'
c     --- C(2) - N3 [neria1996]
      vnlocal(8,1) = 'C(2)'
      vnlocal(8,2) = 'N3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.42d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CH - OH [neria1996]
      vnlocal(1,1) = 'CH'
      vnlocal(1,2) = 'OH'
c     --- C(2) - OH [neria1996]
      vnlocal(2,1) = 'C(2)'
      vnlocal(2,2) = 'OH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (225.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C(2) - C(3) [neria1996]
      vnlocal(1,1) = 'C(2)'
      vnlocal(1,2) = 'C(3)'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.45d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C(2) - CR [neria1996]
      vnlocal(1,1) = 'C(2)'
      vnlocal(1,2) = 'CR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C(2) - S [neria1996]
      vnlocal(1,1) = 'C(2)'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.49d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C(3) - N1 [neria1996]
      vnlocal(1,1) = 'C(3)'
      vnlocal(1,2) = 'N1'
c     --- C(3) - NR [neria1996]
      vnlocal(2,1) = 'C(3)'
      vnlocal(2,2) = 'NR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.77d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- C(3) - S [neria1996]
      vnlocal(1,1) = 'C(3)'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.305d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CR - N1 [neria1996]
      vnlocal(1,1) = 'CR'
      vnlocal(1,2) = 'N1'
c     --- CR - NR [neria1996]
      vnlocal(2,1) = 'CR'
      vnlocal(2,2) = 'NR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.98d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- H - N1 [neria1996]
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'N1'
c     --- H - N2 [neria1996]
      vnlocal(2,1) = 'H'
      vnlocal(2,2) = 'N2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.96d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- H - OH [neria1996]
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'OH'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.98d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HC - N1 [neria1996]
      vnlocal(1,1) = 'HC'
      vnlocal(1,2) = 'N1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.04d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HC - N3 [neria1996]
      vnlocal(1,1) = 'HC'
      vnlocal(1,2) = 'N3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.00d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HC - NC2 [neria1996]
      vnlocal(1,1) = 'HC'
      vnlocal(1,2) = 'NC2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9572d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HT - OT [neria1996]
      vnlocal(1,1) = 'HT'
      vnlocal(1,2) = 'OT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.43d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OC - S [neria1996]
      vnlocal(1,1) = 'OC'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.02d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- S - S [neria1996]
      vnlocal(1,1) = 'S'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)
      
      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * charmm19 angles simple harmonic                                *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - C - C [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'C'
c     --- C - CH - C(3) [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CH'
      bnlocal(2,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - C - C(2) [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'C(2)'
c     --- C - C - C(3) [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - C - CR [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - C - N1 [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - C - NR [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - C - OH [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OH'
c     --- C - C - O [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O'
c     --- CR - C - OH [neria1996]
      bnlocal(3,1) = 'CR'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CH - CH [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'CH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- C - CH - C(2) [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'C(2)'
c     --- C - NR - CR [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'NR'
      bnlocal(2,3) = 'CR'
c     --- C(3) - NR - C [neria1996]
      bnlocal(3,1) = 'C(3)'
      bnlocal(3,2) = 'NR'
      bnlocal(3,3) = 'C'
c     --- C(3) - NR - CR [neria1996]
      bnlocal(4,1) = 'C(3)'
      bnlocal(4,2) = 'NR'
      bnlocal(4,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- C - CH - N [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'N'
c     --- C - CH - N1 [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CH'
      bnlocal(2,3) = 'N1'
c     --- C - CH - N2 [neria1996]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'CH'
      bnlocal(3,3) = 'N2'
c     --- C - CH - N3 [neria1996]
      bnlocal(4,1) = 'C'
      bnlocal(4,2) = 'CH'
      bnlocal(4,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - C(2) - CH [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'CH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - C(2) - C(2) [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'C(2)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - C(2) - N1 [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'N1'
c     --- C - C(2) - N2 [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C(2)'
      bnlocal(2,3) = 'N2'
c     --- C - C(2) - N3 [neria1996]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'C(2)'
      bnlocal(3,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR - C [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.0d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR - C(2) [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'C(2)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.0d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR - CR [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR - N1 [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'N1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR - NR [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - N - CH [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CH'
c     --- C - N - C(2) [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'N'
      bnlocal(2,3) = 'C(2)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (77.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - N1 - CH [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N1'
      bnlocal(1,3) = 'CH'
c     --- C - N1 - C(2) [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'N1'
      bnlocal(2,3) = 'C(2)'
c     --- C - N1 - C(3) [neria1996]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'N1'
      bnlocal(3,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- C - N1 - CR [neria1996]
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'N1'
      bnlocal(isame,3) = 'CR'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - N1 - H [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N1'
      bnlocal(1,3) = 'H'
c     --- C - N2 - H [neria1996]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'N2'
      bnlocal(2,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 7
c     --- C - NC2 - HC [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'HC'
c     --- CH - N1 - H [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'N1'
      bnlocal(2,3) = 'H'
c     --- CH - N2 - H [neria1996]
      bnlocal(3,1) = 'CH'
      bnlocal(3,2) = 'N2'
      bnlocal(3,3) = 'H'
c     --- C(2) - N1 - H [neria1996]
      bnlocal(4,1) = 'C(2)'
      bnlocal(4,2) = 'N1'
      bnlocal(4,3) = 'H'
c     --- C(2) - N2 - H [neria1996]
      bnlocal(5,1) = 'C(2)'
      bnlocal(5,2) = 'N2'
      bnlocal(5,3) = 'H'
c     --- C(3) - N1 - H [neria1996]
      bnlocal(6,1) = 'C(3)'
      bnlocal(6,2) = 'N1'
      bnlocal(6,3) = 'H'
c     --- CR - N1 - H [neria1996]
      bnlocal(7,1) = 'CR'
      bnlocal(7,2) = 'N1'
      bnlocal(7,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 102.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NR - C [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NR'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - OH - H [neria1996]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'OH'
      bnlocal(1,3) = 'H'
c     --- CH - CH - N2 [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'CH'
      bnlocal(2,3) = 'N2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.5d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 8
c     --- CH - C - N [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
c     --- CH - C - N1 [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'N1'
c     --- C(2) - C - N [neria1996]
      bnlocal(3,1) = 'C(2)'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'N'
c     --- C(2) - C - N1 [neria1996]
      bnlocal(4,1) = 'C(2)'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'N1'
c     --- C(2) - C - N2 [neria1996]
      bnlocal(5,1) = 'C(2)'
      bnlocal(5,2) = 'C'
      bnlocal(5,3) = 'N2'
c     --- C(2) - C - NC2 [neria1996]
      bnlocal(6,1) = 'C(2)'
      bnlocal(6,2) = 'C'
      bnlocal(6,3) = 'NC2'
c     --- C(3) - C - N [neria1996]
      bnlocal(7,1) = 'C(3)'
      bnlocal(7,2) = 'C'
      bnlocal(7,3) = 'N'
c     --- C(3) - C - N1 [neria1996]
      bnlocal(8,1) = 'C(3)'
      bnlocal(8,2) = 'C'
      bnlocal(8,3) = 'N1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH - C - O [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
c     --- C(3) - C - O [neria1996]
      bnlocal(2,1) = 'C(3)'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH - C - OC [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CH - C - OH [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OH'
c     --- C(2) - C - OH [neria1996]
      bnlocal(2,1) = 'C(2)'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'OH'
c     --- O - C - OH [neria1996]
      bnlocal(3,1) = 'O'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH - CH - C(2) [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'C(2)'
c     --- CH - C(2) - C(2) [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'C(2)'
      bnlocal(2,3) = 'C(2)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- CH - CH - C(3) [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'C(3)'
c     --- CH - C(2) - OH [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'C(2)'
      bnlocal(2,3) = 'OH'
c     --- C(2) - C(2) - C(3) [neria1996]
      bnlocal(3,1) = 'C(2)'
      bnlocal(3,2) = 'C(2)'
      bnlocal(3,3) = 'C(3)'
c     --- C(3) - C(2) - OH [neria1996]
      bnlocal(4,1) = 'C(3)'
      bnlocal(4,2) = 'C(2)'
      bnlocal(4,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH - CH - N1 [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'N1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH - CH - N3 [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH - CH - OH [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH - C(2) - CH [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'CH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH - C(2) - C(3) [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH - C(2) - S [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'S'
c     --- C(2) - C(2) - S [neria1996]
      bnlocal(2,1) = 'C(2)'
      bnlocal(2,2) = 'C(2)'
      bnlocal(2,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH - N - C(2) [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'C(2)'
c     --- CH - N - C(3) [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'N'
      bnlocal(2,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CH - N1 - C(3) [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'N1'
      bnlocal(1,3) = 'C(3)'
c     --- CH - N2 - C(2) [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'N2'
      bnlocal(2,3) = 'C(2)'
c     --- C(2) - N1 - C(3) [neria1996]
      bnlocal(3,1) = 'C(2)'
      bnlocal(3,2) = 'N1'
      bnlocal(3,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
c     --- CH - N3 - HC [neria1996]
      bnlocal(1,1) = 'CH'
      bnlocal(1,2) = 'N3'
      bnlocal(1,3) = 'HC'
c     --- CH - N3 - C(2) [neria1996]
      bnlocal(2,1) = 'CH'
      bnlocal(2,2) = 'N3'
      bnlocal(2,3) = 'C(2)'
c     --- CH - OH - H [neria1996]
      bnlocal(3,1) = 'CH'
      bnlocal(3,2) = 'OH'
      bnlocal(3,3) = 'H'
c     --- C(2) - OH - H [neria1996]
      bnlocal(4,1) = 'C(2)'
      bnlocal(4,2) = 'OH'
      bnlocal(4,3) = 'H'
c     --- C(2) - N3 - HC [neria1996]
      bnlocal(5,1) = 'C(2)'
      bnlocal(5,2) = 'N3'
      bnlocal(5,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C - CR [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C - NR [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.6d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C - O [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C - 0C [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - CH - CH(3) [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - CH - N [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
c     --- C(2) - CH - N1 [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'N1'
c     --- C(2) - CH - N2 [neria1996]
      bnlocal(2,1) = 'C(2)'
      bnlocal(2,2) = 'CH'
      bnlocal(2,3) = 'N2'
c     --- C(2) - CH - N3 [neria1996]
      bnlocal(3,1) = 'C(2)'
      bnlocal(3,2) = 'CH'
      bnlocal(3,3) = 'N3'
c     --- CR - N1 - CR [neria1996]
      bnlocal(4,1) = 'CR'
      bnlocal(4,2) = 'N1'
      bnlocal(4,3) = 'CR'
c     --- CR - NR - CR [neria1996]
      bnlocal(5,1) = 'CR'
      bnlocal(5,2) = 'NR'
      bnlocal(5,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C(2) - C(2) [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'C(2)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C(2) - N [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C(2) - N1 [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'N1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C(2) - C(2) - N2 [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'N2'
c     --- C(3) - CH - N2 [neria1996]
      bnlocal(2,1) = 'C(3)'
      bnlocal(2,2) = 'CH'
      bnlocal(2,3) = 'N2'
c     --- C(3) - CH - N3 [neria1996]
      bnlocal(3,1) = 'C(3)'
      bnlocal(3,2) = 'CH'
      bnlocal(3,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- C(2) - C(2) - N3 [neria1996]
      isame = isame + 1
      bnlocal(isame,1) = 'C(2)'
      bnlocal(isame,2) = 'C(2)'
      bnlocal(isame,3) = 'N3'
c     --- CR - C - N1 [neria1996]
      isame = isame + 1
      bnlocal(isame,1) = 'CR'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'N1'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - C(2) - N3 [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'C(2)'
      bnlocal(1,3) = 'N3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 99.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - S - C(3) [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.2d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(2) - S - S [neria1996]
      bnlocal(1,1) = 'C(2)'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(3) - CH - C(3) [neria1996]
      bnlocal(1,1) = 'C(3)'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'C(3)'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(3) - CH - N1 [neria1996]
      bnlocal(1,1) = 'C(3)'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'N1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.5d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C(3) - CH - OH [neria1996]
      bnlocal(1,1) = 'C(3)'
      bnlocal(1,2) = 'CH'
      bnlocal(1,3) = 'OH'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CR - C - CR [neria1996]
      bnlocal(1,1) = 'CR'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CR - C - NR [neria1996]
      bnlocal(1,1) = 'CR'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CR - CR - CR [neria1996]
      bnlocal(1,1) = 'CR'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.0d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- H - N2 - H [neria1996]
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'N2'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- HC - N2 - HC [neria1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'N2'
      bnlocal(1,3) = 'HC'
c     --- HC - NC2 - HC [neria1996]
      bnlocal(2,1) = 'HC'
      bnlocal(2,2) = 'NC2'
      bnlocal(2,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HC - N3 - HC [neria1996]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'N3'
      bnlocal(1,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.52d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HT - OT - HT [neria1996]
      bnlocal(1,1) = 'HT'
      bnlocal(1,2) = 'OT'
      bnlocal(1,3) = 'HT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.0d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- N - C - O [neria1996]
      bnlocal(1,1) = 'N'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- N1 - C - NR [neria1996]
      bnlocal(1,1) = 'N1'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
c     --- NC2 - C - N1 [neria1996]
      bnlocal(2,1) = 'NC2'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'N1'
c     --- NC2 - C - NC2 [neria1996]
      bnlocal(3,1) = 'NC2'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'NC2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- N1 - C - O [neria1996]
      bnlocal(1,1) = 'N1'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
c     --- N2 - C - O [neria1996]
      bnlocal(2,1) = 'N2'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- N1 - CR - N1 [neria1996]
      bnlocal(1,1) = 'N1'
      bnlocal(1,2) = 'CR'
      bnlocal(1,3) = 'N1'
c     --- N1 - CR - NR [neria1996]
      bnlocal(2,1) = 'N1'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- OC - C - OC [neria1996]
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- OC - S - OC [neria1996]
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)
      
      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm19 - the original version from the Neria paper           *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (10.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- CH - C - N - CH [neria1996]
      tnlocal(1,1) = 'CH'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CH'
c     --- C(2) - C - N - CH [neria1996]
      tnlocal(2,1) = 'C(2)'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'CH'
c     --- X - C - CR - X [neria1996]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'CR'
      tnlocal(3,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (2.500d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- CR - C - C - C [neria1996]
      tnlocal(1,1) = 'CR'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'C'
c     --- CR - C - C - N1 [neria1996]
      tnlocal(2,1) = 'CR'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'C'
      tnlocal(2,4) = 'N1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (5.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- CR - C - C - CR [neria1996]
      tnlocal(1,1) = 'CR'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'CR'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
c     --- Null torsions
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- X - C - CH - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CH'
      tnlocal(1,4) = 'wild'
c     --- X - C - C(2) - X [neria1996]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'C(2)'
      tnlocal(2,4) = 'wild'
c     --- null - null - null - null [MGM hack]
c     --- this is used for the torsions that are not specified by
c     --- charmm19, but need to be listed for the logic of Towhee
c     --- to function properly
      tnlocal(3,1) = 'null'
      tnlocal(3,2) = 'null'
      tnlocal(3,3) = 'null'
      tnlocal(3,4) = 'null'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (8.200d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- X - C - N - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'wild'
c     --- X - C - N1 - X [neria1996]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'N1'
      tnlocal(2,4) = 'wild'
c     --- X - C - N2 - X [neria1996]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'N2'
      tnlocal(3,4) = 'wild'
c     --- X - C - NC2 - X [neria1996]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'NC2'
      tnlocal(4,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.800d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C - OH - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'OH'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.600d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- X - CH - CH - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH'
      tnlocal(1,3) = 'CH'
      tnlocal(1,4) = 'wild'
c     --- X - CH - C(2) - X [neria1996]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CH'
      tnlocal(2,3) = 'C(2)'
      tnlocal(2,4) = 'wild'
c     --- X - C(2) - C(2) - X [neria1996]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C(2)'
      tnlocal(3,3) = 'C(2)'
      tnlocal(3,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.30d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- X - CH - N - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'wild'
c     --- X - CH - N1 - X [neria1996]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CH'
      tnlocal(2,3) = 'N1'
      tnlocal(2,4) = 'wild'
c     --- X - C(2) - N - X [neria1996]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C(2)'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'wild'
c     --- X - C(2) - N1 - X [neria1996]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'C(2)'
      tnlocal(4,3) = 'N1'
      tnlocal(4,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.800d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - CH - N2 - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH'
      tnlocal(1,3) = 'N2'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.600d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 3
c     --- X - CH - N3 - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH'
      tnlocal(1,3) = 'N3'
      tnlocal(1,4) = 'wild'
c     --- X - C(2) - N2 - X [neria1996]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C(2)'
      tnlocal(2,3) = 'N2'
      tnlocal(2,4) = 'wild'
c     --- X - C(2) - N3 - X [neria1996]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C(2)'
      tnlocal(3,3) = 'N3'
      tnlocal(3,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.500d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - CH - OH - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH'
      tnlocal(1,3) = 'OH'
      tnlocal(1,4) = 'wild'
c     --- X - C(2) - OH - X [neria1996]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C(2)'
      tnlocal(2,3) = 'OH'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (1.200d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - C(2) - S - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C(2)'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (4.000d0)*(CON_KCALTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - S - S - X [neria1996]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- below is charmm19 torsions that are listed in charmm as 
c     --- impropers.  I have modified these from the Neria paper so 
c     --- that they explicitly list all the atoms to avoid any confusion
c     --- with existing regular torsions

      type = type + 1
      call twh_torstyle(GLB_SET,type,1)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(0) = (100.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 4
c     --- C - N1 - CR - N1 (from C - X - X - N1) [neria1996]
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'N1'
      tnlocal(1,3) = 'CR'
      tnlocal(1,4) = 'N1'
c     --- C - NR - CR - N1 (from C - X - X - N1) [neria1996]
      tnlocal(2,1) = 'C'
      tnlocal(2,2) = 'NR'
      tnlocal(2,3) = 'CR'
      tnlocal(2,4) = 'N1'
c     --- C - C - C - N1 (from C - X - X - N1) [neria1996]
      tnlocal(3,1) = 'C'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'C'
      tnlocal(3,4) = 'N1'
c     --- C - C - CR - N1 (from C - X - X - N1) [neria1996]
      tnlocal(4,1) = 'C'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'CR'
      tnlocal(4,4) = 'N1'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,1)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(0) = (25.000d0)*(CON_KCALTOK)
      tclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 21
c     --- N1 - CR - N1 - CR (from CR - X - X - N1) [neria1996]
      tnlocal(1,1) = 'N1'
      tnlocal(1,2) = 'CR'
      tnlocal(1,3) = 'N1'
      tnlocal(1,4) = 'CR'
c     --- CR - N1 - CR - C (from C - X - X - CR) [neria1996]
      tnlocal(2,1) = 'CR'
      tnlocal(2,2) = 'N1'
      tnlocal(2,3) = 'CR'
      tnlocal(2,4) = 'C'
c     --- N1 - CR - C - N1 (from N1 - X - X - N1) [neria1996]
      tnlocal(3,1) = 'N1'
      tnlocal(3,2) = 'CR'
      tnlocal(3,3) = 'C'
      tnlocal(3,4) = 'N1'
c     --- CR - C - N1 - CR (from CR - X - X - CR) [neria1996]
      tnlocal(4,1) = 'CR'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'N1'
      tnlocal(4,4) = 'CR'
c     --- C - N1 - CR - NR (from NR - X - X - C) [neria1996]
      tnlocal(5,1) = 'C'
      tnlocal(5,2) = 'N1'
      tnlocal(5,3) = 'CR'
      tnlocal(5,4) = 'NR'
c     --- N1 - CR - NR - CR (from CR - X - X - N1) [neria1996]
      tnlocal(6,1) = 'N1'
      tnlocal(6,2) = 'CR'
      tnlocal(6,3) = 'NR'
      tnlocal(6,4) = 'CR'
c     --- CR - NR - CR - C (from C - X - X - CR) [neria1996]
      tnlocal(7,1) = 'CR'
      tnlocal(7,2) = 'NR'
      tnlocal(7,3) = 'CR'
      tnlocal(7,4) = 'C'
c     --- NR - CR - C - N1 (from N1 - X - X - NR) [neria1996]
      tnlocal(8,1) = 'NR'
      tnlocal(8,2) = 'CR'
      tnlocal(8,3) = 'C'
      tnlocal(8,4) = 'N1'
c     --- NR - CR - N1 - CR (from NR - X - X - CR) [neria1996]
      tnlocal(9,1) = 'NR'
      tnlocal(9,2) = 'CR'
      tnlocal(9,3) = 'N1'
      tnlocal(9,4) = 'CR'
c     --- N1 - CR - C - NR (from N1 - X - X - NR) [neria1996]
      tnlocal(10,1) = 'N1'
      tnlocal(10,2) = 'CR'
      tnlocal(10,3) = 'C'
      tnlocal(10,4) = 'NR'
c     --- CR - C - NR - CR (from CR - X - X - CR) [neria1996]
      tnlocal(11,1) = 'CR'
      tnlocal(11,2) = 'C'
      tnlocal(11,3) = 'NR'
      tnlocal(11,4) = 'CR'
c     --- C - CR - CR - CR (from C - X - X - CR) [neria1996]
      tnlocal(12,1) = 'C'
      tnlocal(12,2) = 'CR'
      tnlocal(12,3) = 'CR'
      tnlocal(12,4) = 'CR'
c     --- CR - CR - CR - CR (from CR - X - X - CR) [neria1996]
      tnlocal(13,1) = 'CR'
      tnlocal(13,2) = 'CR'
      tnlocal(13,3) = 'CR'
      tnlocal(13,4) = 'CR'
c     --- CR - C - CR - CR (from CR - X - X - CR) [neria1996]
      tnlocal(14,1) = 'CR'
      tnlocal(14,2) = 'C'
      tnlocal(14,3) = 'CR'
      tnlocal(14,4) = 'CR'
c     --- C - CR - N1 - C (from C - X - X - C) [neria1996]
      tnlocal(15,1) = 'C'
      tnlocal(15,2) = 'CR'
      tnlocal(15,3) = 'N1'
      tnlocal(15,4) = 'C'
c     --- C - C - C - CR (from C - X - X - CR) [neria1996]
      tnlocal(16,1) = 'C'
      tnlocal(16,2) = 'C'
      tnlocal(16,3) = 'C'
      tnlocal(16,4) = 'CR'
c     --- N1 - C - CR - CR (from N1 - C - CR - CR) [neria1996]
      tnlocal(17,1) = 'N1'
      tnlocal(17,2) = 'C'
      tnlocal(17,3) = 'CR'
      tnlocal(17,4) = 'CR'
c     --- CR - CR - C - C (from CR - X - X - C) [neria1996]
      tnlocal(18,1) = 'CR'
      tnlocal(18,2) = 'CR'
      tnlocal(18,3) = 'C'
      tnlocal(18,4) = 'C'
c     --- CR - C - C - N1 (from CR - X - X - N1) [neria1996]
      tnlocal(19,1) = 'CR'
      tnlocal(19,2) = 'C'
      tnlocal(19,3) = 'C'
      tnlocal(19,4) = 'N1'
c     --- CR - C - C - CR (from CR - X - X - CR) [neria1996]
      tnlocal(20,1) = 'CR'
      tnlocal(20,2) = 'C'
      tnlocal(20,3) = 'C'
      tnlocal(20,4) = 'CR'
c     --- C - CR - CR - C (from C - X - X - C) [neria1996]
      tnlocal(21,1) = 'C'
      tnlocal(21,2) = 'CR'
      tnlocal(21,3) = 'CR'
      tnlocal(21,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)
      
      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm19 improper harmonic torsion                             *
c     ******************************************************************

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (90.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 7
c     --- C - C - CR - C(2) [neria1996]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'C'
      ninlocal(1,3) = 'CR'
      ninlocal(1,4) = 'C(2)'
c     --- C - CR - C - C(2) [neria1996]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'CR'
      ninlocal(2,3) = 'C'
      ninlocal(2,4) = 'C(2)'
c     --- C - CR - CR - C(2) [neria1996]
      ninlocal(3,1) = 'C'
      ninlocal(3,2) = 'CR'
      ninlocal(3,3) = 'CR'
      ninlocal(3,4) = 'C(2)'
c     --- C - CR - N1 - C(2) [neria1996]
      ninlocal(4,1) = 'C'
      ninlocal(4,2) = 'CR'
      ninlocal(4,3) = 'N1'
      ninlocal(4,4) = 'C(2)'
c     --- C - N1 - CR - C(2) [neria1996]
      ninlocal(5,1) = 'C'
      ninlocal(5,2) = 'N1'
      ninlocal(5,3) = 'CR'
      ninlocal(5,4) = 'C(2)'
c     --- C - X - X - C(2) [neria1996]
      ninlocal(6,1) = 'C'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'C(2)'
c     --- C - X - X - C(3) [neria1996]
      ninlocal(7,1) = 'C'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'wild'
      ninlocal(7,4) = 'C(3)'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (150.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- C - CR - CR - OH [neria1996]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'CR'
      ninlocal(1,3) = 'CR'
      ninlocal(1,4) = 'OH'
c     --- C - X - X - OH [neria1996]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'OH'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (45.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 10
c     --- C - H - H - N2 [neria1996]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'H'
      ninlocal(1,3) = 'H'
      ninlocal(1,4) = 'N2'
c     --- H - X - X - O [neria1996]
      ninlocal(2,1) = 'H'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'O'
c     --- N - CH - C(2) - C [neria1996]
      ninlocal(3,1) = 'N'
      ninlocal(3,2) = 'CH'
      ninlocal(3,3) = 'C(2)'
      ninlocal(3,4) = 'C'
c     --- N - X - X - C(2) [neria1996]
      ninlocal(4,1) = 'N'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'C(2)'
c     --- N1 - X - X - CH [neria1996]
      ninlocal(5,1) = 'N1'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'CH'
c     --- N1 - X - X - C(2) [neria1996]
      ninlocal(6,1) = 'N1'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'C(2)'
c     --- N1 - X - X - C(3) [neria1996]
      ninlocal(7,1) = 'N1'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'wild'
      ninlocal(7,4) = 'C(3)'
c     --- N1 - X - X - H [neria1996]
      ninlocal(8,1) = 'N1'
      ninlocal(8,2) = 'wild'
      ninlocal(8,3) = 'wild'
      ninlocal(8,4) = 'H'
c     --- N2 - X - X - H [neria1996]
      ninlocal(9,1) = 'N2'
      ninlocal(9,2) = 'wild'
      ninlocal(9,3) = 'wild'
      ninlocal(9,4) = 'H'
c     --- NC2 - X - X - HC [neria1996]
      ninlocal(10,1) = 'NC2'
      ninlocal(10,2) = 'wild'
      ninlocal(10,3) = 'wild'
      ninlocal(10,4) = 'HC'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (100.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 5
c     --- C - OC - OC - CH [neria1996]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'OC'
      ninlocal(1,3) = 'OC'
      ninlocal(1,4) = 'CH'
c     --- C - OC - OC - C(2) [neria1996]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'OC'
      ninlocal(2,3) = 'OC'
      ninlocal(2,4) = 'C(2)'
c     --- C - X - X - N1 [neria1996]
      ninlocal(3,1) = 'C'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'N1'
c     --- C - X - X - O [neria1996]
      ninlocal(4,1) = 'C'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'O'
c     --- C - X - X - OC [neria1996]
      ninlocal(5,1) = 'C'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'OC'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (25.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 9
c     --- C - X - X - C [neria1996]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'C'
c     --- C - X - X - CR [neria1996]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'CR'
c     --- CR - X - X - CR [neria1996]
      ninlocal(3,1) = 'CR'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'CR'
c     --- CR - X - X - N1 [neria1996]
      ninlocal(4,1) = 'CR'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'N1'
c     --- N1 - X - X - N1 [neria1996]
      ninlocal(5,1) = 'N1'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'N1'
c     --- N1 - X - X - NR [neria1996]
      ninlocal(6,1) = 'N1'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'NR'
c     --- NR - X - X - C [neria1996]
      ninlocal(7,1) = 'NR'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'wild'
      ninlocal(7,4) = 'C'
c     --- NR - X - X - CR [neria1996]
      ninlocal(8,1) = 'NR'
      ninlocal(8,2) = 'wild'
      ninlocal(8,3) = 'wild'
      ninlocal(8,4) = 'CR'
c     --- NR - X - X - C(3) [neria1996]
      ninlocal(9,1) = 'NR'
      ninlocal(9,2) = 'wild'
      ninlocal(9,3) = 'wild'
      ninlocal(9,4) = 'C(3)'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (75.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- C - X - X - H [neria1996]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'H'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (55.00d0)*(CON_KCALTOK)
      iclocal(1) = (35.26439d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 3
c     --- CH - X - X - CH [neria1996]
      ninlocal(1,1) = 'CH'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'CH'
c     --- CH - X - X - C(2) [neria1996]
      ninlocal(2,1) = 'CH'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'C(2)'
c     --- CH - X - X - C(3) [neria1996]
      ninlocal(3,1) = 'CH'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'C(3)'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end
