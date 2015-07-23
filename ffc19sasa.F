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
      subroutine twh_ffc19sasa(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the c19sasa      *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-21-2014 by M.G. martin                        *
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
      integer ii,type
      integer itrial,i,j,jlo,jhi,jj,ij,isame
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal,ntslocal,icoeff
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

      ffname = 'C19sasa'
      type = 0
      call twh_classical_potential(GLB_SET,'12-6 plus solvation')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     **************************************************************
c     * C19SASA force field                                        *
c     * contains all of the parameters listed in the parameter     *
c     * that is issued with Charmm28 "param19_eef1.inp"            *
c     * and all of the solvation parameters are from               *
c     * P. Ferrara, J. Apostolakis, and A. Caflisch,               *
c     * Proteins 46, 24-33 (2002).                                 *
c     **************************************************************

c     --- we do not use coefficients 8-11 even though they are output
c     --- for this nonbonded potential.  Zero them out
      nbclocal(8) = 0.0d0
      nbclocal(9) = 0.0d0
      nbclocal(10) = 0.0d0
      nbclocal(11) = 0.0d0

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.100d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1200d0)*(CON_KCALTOK)
      nbclocal(3) = (1.9d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(4) = (0.100d0)*(CON_KCALTOK)
c     Solvation parameters
      nbclocal(5) = 1.720d0
      nbclocal(6) = 1.554d0
      nbclocal(7) = 0.012d0
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
c     Solvation parameters
      nbclocal(5) = 1.800d0
      nbclocal(6) = 1.276d0
      nbclocal(7) = 0.012d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH1E: Extended aliphatic carbon with one hydrogen -CH- group
c     [param19_eef1.inp, (same as CH - neria1996)]
      nbnlocal ='CH1E'
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
c     Solvation parameters
      nbclocal(5) = 1.900d0
      nbclocal(6) = 1.045d0
      nbclocal(7) = 0.012d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH2E: Extended aliphatic carbon with two hydrogens -CH2- group
c     [param19_eef1.inp, (same as C(2) - neria1996)]
      nbnlocal ='CH2E'
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
c     Solvation parameters
      nbclocal(5) = 2.000d0
      nbclocal(6) = 0.880d0
      nbclocal(7) = 0.012d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CH3E: Extended methyl terminal -CH3 group 
c     [param19_eef1.inp, (same as C(3) - neria1996)]
      nbnlocal ='CH3E'
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
c     Solvation parameters
      nbclocal(5) = 1.800d0
      nbclocal(6) = 1.073d0
      nbclocal(7) = 0.012d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     CR1E: Extended aromatic carbon with one hydrogen 
c     [param19_eef1.inp, (same as CR - neria1996)]
      nbnlocal ='CR1E'
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
c     Solvation parameters
      nbclocal(5) = 1.550d0
      nbclocal(6) = 1.028d0
      nbclocal(7) = -0.060d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     NH1: Nitrogen bound to one hydrogen 
c     [param19_eef1.inp (same as N1 - neria1996)]
      nbnlocal ='NH1'
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
c     Solvation parameters
      nbclocal(5) = 1.550d0
      nbclocal(6) = 1.028d0
      nbclocal(7) = -0.060d0
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
      nbclocal(2) = (0.2384d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
c     Solvation parameters
      nbclocal(5) = 1.600d0
      nbclocal(6) = 1.215d0
      nbclocal(7) = -0.060d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     NH2: Nitrogen bound to two hydrogens 
c     [param19_eef1.inp (same as N2 - neria1996)]
      nbnlocal ='NH2'
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
c     Solvation parameters
      nbclocal(5) = 1.600d0
      nbclocal(6) = 1.215d0
      nbclocal(7) = -0.060d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     NH3: Nitrogen bound to three hydrogens 
c     [param19_eef1.inp (same as N3 - neria1996)]
      nbnlocal ='NH3'
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
c     Solvation parameters
      nbclocal(5) = 1.550d0
      nbclocal(6) = 1.028d0
      nbclocal(7) = -0.060d0
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
c     Solvation parameters
      nbclocal(5) = 1.550d0
      nbclocal(6) = 1.028d0
      nbclocal(7) = -0.060d0
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
      nbclocal(2) = (0.1591d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
c     Solvation parameters
      nbclocal(5) = 1.520d0
      nbclocal(6) = 1.080d0
      nbclocal(7) = -0.060d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     OH1: Hydroxyl oxygen 
c     [param19_eef1.inp (same as OH - neria1996)] 
      nbnlocal ='OH1'
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
c     Solvation parameters
      nbclocal(5) = 1.500d0
      nbclocal(6) = 0.926d0
      nbclocal(7) = -0.060d0
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
c     Solvation parameters
      nbclocal(5) = 1.700d0
      nbclocal(6) = 0.922d0
      nbclocal(7) = -0.060d0
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
      nbclocal(1) = (1.890d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0430d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
c     Solvation parameters
      nbclocal(5) = 1.800d0
      nbclocal(6) = 1.121d0
      nbclocal(7) = 0.012d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     S: sulfur [param19_eef1.inp]
      nbnlocal ='S'
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
c     Solvation parameters
      nbclocal(5) = 1.800d0
      nbclocal(6) = 1.121d0
      nbclocal(7) = 0.012d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     SH1E: extended sulfur with one hydrogen
c     [param19_eef1.inp]
      nbnlocal ='SH1E'
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
c     Solvation parameters
      nbclocal(5) = 1.100d0
      nbclocal(6) = 1.128d0
      nbclocal(7) = 0.000d0
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
c     Solvation parameters
      nbclocal(5) = 1.100d0
      nbclocal(6) = 1.128d0
      nbclocal(7) = 0.000d0
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
c     CR: aromatic carbon with no hydrogens (found in the amino 
c     acids which have a ring (TYR, PHE, HIS, TRP, ARG))      
c     [param19_eef1.inp, introduced in lazaridis1999]
      nbnlocal ='CR'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.490d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0262d0)*(CON_KCALTOK)
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
c     CT: aliphatic carbon [param19_eef1.inp]
      nbnlocal ='CT'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.490d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0262d0)*(CON_KCALTOK)
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
c     CM: Carbon monoxide carbon [param19_eef1.inp]
      nbnlocal ='CM'
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
      nbclocal(1) = (1.4680d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0450d0)*(CON_KCALTOK)
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
c     HA: aliphatic hydrogen [param19_eef1.inp]
      nbnlocal ='HA'
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
c     NP: pyrole nitrogen [param19_eef1.inp]
      nbnlocal ='NP'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (1.7398d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0758d0)*(CON_KCALTOK)
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
c     OH2: ST2 water oxygen [param19_eef1.inp]
      nbnlocal ='OH2'
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
c     OT: water oxygen, modified TIP3P [neria1996]
      nbnlocal ='OT'
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
c     OS: esther oxygen [param19_eef1.inp]
      nbnlocal ='OS'
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
c     OM: Oxygen in carbon monoxide [param19_eef1.inp]
      nbnlocal ='OM'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.650d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 32.066d0
      elelocal = 'Fe'
      bplocal = 'sp'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     FE: iron
c     [param19_eef1.inp]
      nbnlocal ='FE'
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
c              --- reset these parameters manually
               do icoeff = 1,4
                  call twh_nbcoeff(GLB_SET,icoeff,ij,nbclocal(icoeff))
               enddo
            enddo
         enddo
      enddo

c     ******************************************************************
c     * Charmm19                                                       *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.52d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 6
c     --- C - CH1E [param19_eef1.inp]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'CH1E'
c     --- C - CH2E [param19_eef1.inp]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'CH2E'
c     --- C - CH3E [param19_eef1.inp]
      vnlocal(3,1) = 'C'
      vnlocal(3,2) = 'CH3E'
c     --- CR - CH1E [param19_eef1.inp]
      vnlocal(4,1) = 'CR'
      vnlocal(4,2) = 'CH1E'
c     --- CR - CH2E [param19_eef1.inp]
      vnlocal(5,1) = 'CR'
      vnlocal(5,2) = 'CH2E'
c     --- CR - CH3E [param19_eef1.inp]
      vnlocal(6,1) = 'CR'
      vnlocal(6,2) = 'CH3E'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 8
c     --- C - C [param19_eef1.inp]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'C'
c     --- C - CR1E [param19_eef1.inp]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'CR1E'
c     --- C - OH1 [param19_eef1.inp]
      vnlocal(3,1) = 'C'
      vnlocal(3,2) = 'OH1'
c     --- CR1E - CR1E [param19_eef1.inp]
      vnlocal(4,1) = 'CR1E'
      vnlocal(4,2) = 'CR1E'
c     --- CR - CR [param19_eef1.inp]
      vnlocal(5,1) = 'CR'
      vnlocal(5,2) = 'CR'
c     --- CR - C [param19_eef1.inp]
      vnlocal(6,1) = 'CR'
      vnlocal(6,2) = 'C'
c     --- CR - CR1E [param19_eef1.inp]
      vnlocal(7,1) = 'CR'
      vnlocal(7,2) = 'CR1E'
c     --- CR - OH1 [param19_eef1.inp]
      vnlocal(8,1) = 'CR'
      vnlocal(8,2) = 'OH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.53d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C - CT [param19_eef1.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'CT'
c     --- C - CT [param19_eef1.inp]
      isame = isame + 1
      vnlocal(isame,1) = 'CR'
      vnlocal(isame,2) = 'CT'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.33d0
      vclocal(1) = (471.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 10
c     --- C - N [param19_eef1.inp]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'N'
c     --- C - NH1 [param19_eef1.inp]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'NH1'
c     --- C - NH2 [param19_eef1.inp]
      vnlocal(3,1) = 'C'
      vnlocal(3,2) = 'NH2'
c     --- C - NR [param19_eef1.inp]
      vnlocal(4,1) = 'C'
      vnlocal(4,2) = 'NR'
c     --- C - NP [param19_eef1.inp]
      vnlocal(5,1) = 'C'
      vnlocal(5,2) = 'NP'
c     --- CR - N [param19_eef1.inp]
      vnlocal(6,1) = 'CR'
      vnlocal(6,2) = 'N'
c     --- CR - NH1 [param19_eef1.inp]
      vnlocal(7,1) = 'CR'
      vnlocal(7,2) = 'NH1'
c     --- CR - NH2 [param19_eef1.inp]
      vnlocal(8,1) = 'CR'
      vnlocal(8,2) = 'NH2'
c     --- CR - NR [param19_eef1.inp]
      vnlocal(9,1) = 'CR'
      vnlocal(9,2) = 'NR'
c     --- CR - NP [param19_eef1.inp]
      vnlocal(10,1) = 'CR'
      vnlocal(10,2) = 'NP'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.33d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C - NC2 [param19_eef1.inp]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'NC2'
c     --- CR - NC2 [param19_eef1.inp]
      vnlocal(2,1) = 'CR'
      vnlocal(2,2) = 'NC2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.23d0
      vclocal(1) = (580.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 4
c     --- C - O [param19_eef1.inp]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'O'
c     --- C - OC [param19_eef1.inp]
      vnlocal(2,1) = 'C'
      vnlocal(2,2) = 'OC'
c     --- CR - O [param19_eef1.inp]
      vnlocal(3,1) = 'CR'
      vnlocal(3,2) = 'O'
c     --- CR - OC [param19_eef1.inp]
      vnlocal(4,1) = 'CR'
      vnlocal(4,2) = 'OC'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.43d0
      vclocal(1) = (292.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- C - OS [param19_eef1.inp]
      vnlocal(1,1) = 'C'
      vnlocal(1,2) = 'OS'
c     --- CR - OS [param19_eef1.inp]
      vnlocal(2,1) = 'CR'
      vnlocal(2,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.53d0
      vclocal(1) = (225.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CH1E - CH1E [param19_eef1.inp]
      vnlocal(1,1) = 'CH1E'
      vnlocal(1,2) = 'CH1E'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.52d0
      vclocal(1) = (225.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CH1E - CH2E [param19_eef1.inp]
      vnlocal(1,1) = 'CH1E'
      vnlocal(1,2) = 'CH2E'
c     --- CH1E - CH3E [param19_eef1.inp]
      vnlocal(2,1) = 'CH1E'
      vnlocal(2,2) = 'CH3E'
c     --- CH2E - CH2E [param19_eef1.inp]
      vnlocal(3,1) = 'CH2E'
      vnlocal(3,2) = 'CH2E'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.45d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 13
c     --- CH1E - N [param19_eef1.inp]
      vnlocal(1,1) = 'CH1E'
      vnlocal(1,2) = 'N'
c     --- CH1E - NH1 [param19_eef1.inp]
      vnlocal(2,1) = 'CH1E'
      vnlocal(2,2) = 'NH1'
c     --- CH1E - NH2 [param19_eef1.inp]
      vnlocal(3,1) = 'CH1E'
      vnlocal(3,2) = 'NH2'
c     --- CH1E - NH3 [param19_eef1.inp]
      vnlocal(4,1) = 'CH1E'
      vnlocal(4,2) = 'NH3'
c     --- CH2E - N [param19_eef1.inp]
      vnlocal(5,1) = 'CH2E'
      vnlocal(5,2) = 'N'
c     --- CH2E - NH1 [param19_eef1.inp]
      vnlocal(6,1) = 'CH2E'
      vnlocal(6,2) = 'NH1'
c     --- CH2E - NH2 [param19_eef1.inp]
      vnlocal(7,1) = 'CH2E'
      vnlocal(7,2) = 'NH2'
c     --- CH2E - NH3 [param19_eef1.inp]
      vnlocal(8,1) = 'CH2E'
      vnlocal(8,2) = 'NH3'
c     --- CT - N [param19_eef1.inp]
      vnlocal(9,1) = 'CT'
      vnlocal(9,2) = 'N'
c     --- CT - NH1 [param19_eef1.inp]
      vnlocal(10,1) = 'CT'
      vnlocal(10,2) = 'NH1'
c     --- CT - NH2 [param19_eef1.inp]
      vnlocal(11,1) = 'CT'
      vnlocal(11,2) = 'NH2'
c     --- CT - NH3 [param19_eef1.inp]
      vnlocal(12,1) = 'CT'
      vnlocal(12,2) = 'NH3'
c     --- CT - NC2 [param19_eef1.inp]
      vnlocal(13,1) = 'CT'
      vnlocal(13,2) = 'NC2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.42d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CH1E - OH1 [param19_eef1.inp]
      vnlocal(1,1) = 'CH1E'
      vnlocal(1,2) = 'OH1'
c     --- CH2E - OH1 [param19_eef1.inp]
      vnlocal(2,1) = 'CH2E'
      vnlocal(2,2) = 'OH1'
c     --- CT - OH1 [param19_eef1.inp]
      vnlocal(3,1) = 'CT'
      vnlocal(3,2) = 'OH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.54d0
      vclocal(1) = (225.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CH2E - CH3E [param19_eef1.inp]
      vnlocal(1,1) = 'CH2E'
      vnlocal(1,2) = 'CH3E'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.45d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CH2E - CR1E [param19_eef1.inp]
      vnlocal(1,1) = 'CH2E'
      vnlocal(1,2) = 'CR1E'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- CH2E - SH1E [param19_eef1.inp]
      vnlocal(1,1) = 'CH2E'
      vnlocal(1,2) = 'SH1E'
c     --- CH2E - S [param19_eef1.inp]
      vnlocal(2,1) = 'CH2E'
      vnlocal(2,2) = 'S'
c     --- CT - S [param19_eef1.inp]
      vnlocal(3,1) = 'CT'
      vnlocal(3,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.49d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CH3E - NH1 [param19_eef1.inp]
      vnlocal(1,1) = 'CH3E'
      vnlocal(1,2) = 'NH1'
c     --- CH3E - NR [param19_eef1.inp]
      vnlocal(2,1) = 'CH3E'
      vnlocal(2,2) = 'NR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.77d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CH3E - S [param19_eef1.inp]
      vnlocal(1,1) = 'CH3E'
      vnlocal(1,2) = 'S'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (292.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CH3E - OS [param19_eef1.inp]
      vnlocal(1,1) = 'CH3E'
      vnlocal(1,2) = 'OS'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.128d0
      vclocal(1) = (1115.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CM - OM [param19_eef1.inp]
      vnlocal(1,1) = 'CM'
      vnlocal(1,2) = 'OM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.305d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
c     --- CR1E - NH1 [param19_eef1.inp]
      vnlocal(1,1) = 'CR1E'
      vnlocal(1,2) = 'NH1'
c     --- CR1E - NR [param19_eef1.inp]
      vnlocal(2,1) = 'CR1E'
      vnlocal(2,2) = 'NR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.53d0
      vclocal(1) = (200.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- CT - CT [param19_eef1.inp]
      vnlocal(1,1) = 'CT'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.79d0
      vclocal(1) = (258.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- FE - CM [param19_eef1.inp]
      vnlocal(1,1) = 'FE'
      vnlocal(1,2) = 'CM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.09d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- FE - NP [param19_eef1.inp]
      vnlocal(1,1) = 'FE'
      vnlocal(1,2) = 'NP'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.98d0
      vclocal(1) = (65.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- FE - NR [param19_eef1.inp]
      vnlocal(1,1) = 'FE'
      vnlocal(1,2) = 'NR'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.80d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- FE - OM [param19_eef1.inp]
      vnlocal(1,1) = 'FE'
      vnlocal(1,2) = 'OM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.96d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- H - OH1 [param19_eef1.inp]
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'OH1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.00d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- H - OH2 [param19_eef1.inp]
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'OH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.6329931d0
      vclocal(1) = (0.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- H - H [param19_eef1.inp]
      vnlocal(1,1) = 'H'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.98d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
c     --- HC - NH1 [param19_eef1.inp]
      vnlocal(1,1) = 'HC'
      vnlocal(1,2) = 'NH1'
c     --- H - NH1 [param19_eef1.inp]
      vnlocal(2,1) = 'H'
      vnlocal(2,2) = 'NH1'
c     --- H - NH2 [param19_eef1.inp]
      vnlocal(3,1) = 'H'
      vnlocal(3,2) = 'NH2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.04d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HC - NH3 [param19_eef1.inp]
      vnlocal(1,1) = 'HC'
      vnlocal(1,2) = 'NH3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.00d0
      vclocal(1) = (405.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HC - NC2 [param19_eef1.inp]
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
c     --- HT - OT [param19_eef1.inp]
      vnlocal(1,1) = 'HT'
      vnlocal(1,2) = 'OT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5139d0
      vclocal(1) = (0.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HT - HT [param19_eef1.inp]
      vnlocal(1,1) = 'HT'
      vnlocal(1,2) = 'HT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (350.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HA - C [param19_eef1.inp]
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'C'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0800d0
      vclocal(1) = (300.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- HA - CT [param19_eef1.inp]
      vnlocal(1,1) = 'HA'
      vnlocal(1,2) = 'CT'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.23d0
      vclocal(1) = (600.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OM - OM [param19_eef1.inp]
      vnlocal(1,1) = 'OM'
      vnlocal(1,2) = 'OM'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.43d0
      vclocal(1) = (400.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- OC - S [param19_eef1.inp]
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
c     --- S - S [param19_eef1.inp]
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
      nbslocal = 3
c     --- C - C - C [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'C'
c     --- C - CH1E - CH3E [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CH1E'
      bnlocal(2,3) = 'CH3E'
c     --- CR - CR - CR [param19_eef1.inp]
      bnlocal(3,1) = 'CR'
      bnlocal(3,2) = 'CR'
      bnlocal(3,3) = 'CR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - C - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CH2E'
c     --- CR - CR - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'CH2E'
c     --- C - C - CH3E [param19_eef1.inp]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - C - CR1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CR1E'
c     --- CR - CR - CR1E [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - C - CT [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - C - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NH1'
c     --- CR - CR - NH1 [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - C - NR [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
c     --- C - C - NP [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NP'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - C - OH1 [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OH1'
c     --- C - C - O [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O'
c     --- CR1E - C - OH1 [param19_eef1.inp]
      bnlocal(3,1) = 'CR1E'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'OH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CH1E - CH1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'CH1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 7
c     --- C - CH1E - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'CH2E'
c     --- C - NR - CR1E [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'NR'
      bnlocal(2,3) = 'CR1E'
c     --- CH3E - NR - C [param19_eef1.inp]
      bnlocal(3,1) = 'CH3E'
      bnlocal(3,2) = 'NR'
      bnlocal(3,3) = 'C'
c     --- CH3E - NR - CR1E [param19_eef1.inp]
      bnlocal(4,1) = 'CH3E'
      bnlocal(4,2) = 'NR'
      bnlocal(4,3) = 'CR1E'
c     --- C - CT - CT [param19_eef1.inp]
      bnlocal(5,1) = 'C'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'CT'
c     --- C - CT - HA [param19_eef1.inp]
      bnlocal(6,1) = 'C'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'HA'
c     --- CR - NR - CR1E [param19_eef1.inp]
      bnlocal(7,1) = 'CR'
      bnlocal(7,2) = 'NR'
      bnlocal(7,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- C - CH1E - N [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'N'
c     --- C - CH1E - NH1 [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CH1E'
      bnlocal(2,3) = 'NH1'
c     --- C - CH1E - NH2 [param19_eef1.inp]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'CH1E'
      bnlocal(3,3) = 'NH2'
c     --- C - CH1E - NH3 [param19_eef1.inp]
      bnlocal(4,1) = 'C'
      bnlocal(4,2) = 'CH1E'
      bnlocal(4,3) = 'NH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - CH2E - CH1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'CH1E'
c     --- CR - CH2E - CH1E [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CH2E'
      bnlocal(2,3) = 'CH1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.0d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CH2E - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'CH2E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 7
c     --- C - CH2E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'NH1'
c     --- C - CH2E - NH2 [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'CH2E'
      bnlocal(2,3) = 'NH2'
c     --- C - CH2E - NH3 [param19_eef1.inp]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'CH2E'
      bnlocal(3,3) = 'NH3'
c     --- C - CT - N [param19_eef1.inp]
      bnlocal(4,1) = 'C'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'N'
c     --- C - CT - NH1 [param19_eef1.inp]
      bnlocal(5,1) = 'C'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'NH1'
c     --- C - CT - NH2 [param19_eef1.inp]
      bnlocal(6,1) = 'C'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'NH2'
c     --- C - CT - NH3 [param19_eef1.inp]
      bnlocal(7,1) = 'C'
      bnlocal(7,2) = 'CT'
      bnlocal(7,3) = 'NH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 126.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR1E - C [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.0d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - CR1E - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'CH2E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 119.0d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - CR1E - CR1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'CR1E'
c     --- CR - CR1E - CR1E [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CR1E'
      bnlocal(2,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - CR1E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'NH1'
c     --- CR - CR1E - NH1 [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CR1E'
      bnlocal(2,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - CR1E - NR [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'NR'
c     --- CR - CR1E - NR [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'CR1E'
      bnlocal(2,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (80.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
c     --- C - N - CH1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CH1E'
c     --- C - N - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'N'
      bnlocal(2,3) = 'CH2E'
c     --- C - N - CT [param19_eef1.inp]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'N'
      bnlocal(3,3) = 'CT'
c     --- C - NC2 - CT [param19_eef1.inp]
      bnlocal(4,1) = 'C'
      bnlocal(4,2) = 'NC2'
      bnlocal(4,3) = 'CT'
c     --- C - NH1 - CT [param19_eef1.inp]
      bnlocal(5,1) = 'C'
      bnlocal(5,2) = 'NH1'
      bnlocal(5,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (77.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- C - NH1 - CH1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'CH1E'
c     --- C - NH1 - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'NH1'
      bnlocal(2,3) = 'CH2E'
c     --- C - NH1 - CH3E [param19_eef1.inp]
      bnlocal(3,1) = 'C'
      bnlocal(3,2) = 'NH1'
      bnlocal(3,3) = 'CH3E'
c     --- CR - NH1 - CH2E [param19_eef1.inp]
      bnlocal(4,1) = 'CR'
      bnlocal(4,2) = 'NH1'
      bnlocal(4,3) = 'CH2E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - NH1 - CR1E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'CR1E'
c     --- CR - NH1 - CR1E [param19_eef1.inp]
      bnlocal(2,1) = 'CR'
      bnlocal(2,2) = 'NH1'
      bnlocal(2,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- C - NH1 - H [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'H'
c     --- C - NH2 - H [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'NH2'
      bnlocal(2,3) = 'H'
c     --- CR - NH1 - H [param19_eef1.inp]
      bnlocal(3,1) = 'CR'
      bnlocal(3,2) = 'NH1'
      bnlocal(3,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 10
c     --- C - NC2 - HC [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'HC'
c     --- CH1E - NH1 - H [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'NH1'
      bnlocal(2,3) = 'H'
c     --- CH1E - NH2 - H [param19_eef1.inp]
      bnlocal(3,1) = 'CH1E'
      bnlocal(3,2) = 'NH2'
      bnlocal(3,3) = 'H'
c     --- CH2E - NH1 - H [param19_eef1.inp]
      bnlocal(4,1) = 'CH2E'
      bnlocal(4,2) = 'NH1'
      bnlocal(4,3) = 'H'
c     --- CH2E - NH2 - H [param19_eef1.inp]
      bnlocal(5,1) = 'CH2E'
      bnlocal(5,2) = 'NH2'
      bnlocal(5,3) = 'H'
c     --- CH3E - NH1 - H [param19_eef1.inp]
      bnlocal(6,1) = 'CH3E'
      bnlocal(6,2) = 'NH1'
      bnlocal(6,3) = 'H'
c     --- CR1E - NH1 - H [param19_eef1.inp]
      bnlocal(7,1) = 'CR1E'
      bnlocal(7,2) = 'NH1'
      bnlocal(7,3) = 'H'
c     --- CR - NC2 - HC [param19_eef1.inp]
      bnlocal(8,1) = 'CR'
      bnlocal(8,2) = 'NC2'
      bnlocal(8,3) = 'HC'
c     --- CT - NC2 - HC [param19_eef1.inp]
      bnlocal(9,1) = 'CT'
      bnlocal(9,2) = 'NC2'
      bnlocal(9,3) = 'HC'
c     --- CT - NH1 - H [param19_eef1.inp]
      bnlocal(10,1) = 'CT'
      bnlocal(10,2) = 'NH1'
      bnlocal(10,3) = 'H'
c     --- CT - NH2 - H [param19_eef1.inp]
      bnlocal(11,1) = 'CT'
      bnlocal(11,2) = 'NH2'
      bnlocal(11,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 102.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- C - NR - C [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NR'
      bnlocal(1,3) = 'C'
c     --- C - NP - C [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'NP'
      bnlocal(2,3) = 'C'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 8
c     --- C - OH1 - H [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'OH1'
      bnlocal(1,3) = 'H'
c     --- CH1E - CH1E - NH2 [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'CH1E'
      bnlocal(2,3) = 'NH2'
c     --- CT - CT - OH1 [param19_eef1.inp]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'OH1'
c     --- HA - CT - N [param19_eef1.inp]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'N'
c     --- HA - CT - NC2 [param19_eef1.inp]
      bnlocal(5,1) = 'HA'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'NC2'
c     --- HA - CT - NH1 [param19_eef1.inp]
      bnlocal(6,1) = 'HA'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'NH1'
c     --- HA - CT - NH3 [param19_eef1.inp]
      bnlocal(7,1) = 'HA'
      bnlocal(7,2) = 'CT'
      bnlocal(7,3) = 'NH3'
c     --- HA - CT - OH1 [param19_eef1.inp]
      bnlocal(8,1) = 'HA'
      bnlocal(8,2) = 'CT'
      bnlocal(8,3) = 'OH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.5d0
      bclocal(1) = (20.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 12
c     --- CH1E - C - N [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'N'
c     --- CH1E - C - NH1 [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NH1'
c     --- CH2E - C - N [param19_eef1.inp]
      bnlocal(3,1) = 'CH2E'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'N'
c     --- CH2E - C - NH1 [param19_eef1.inp]
      bnlocal(4,1) = 'CH2E'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'NH1'
c     --- CH2E - C - NH2 [param19_eef1.inp]
      bnlocal(5,1) = 'CH2E'
      bnlocal(5,2) = 'C'
      bnlocal(5,3) = 'NH2'
c     --- CH2E - C - NC2 [param19_eef1.inp]
      bnlocal(6,1) = 'CH2E'
      bnlocal(6,2) = 'C'
      bnlocal(6,3) = 'NC2'
c     --- CH3E - C - N [param19_eef1.inp]
      bnlocal(7,1) = 'CH3E'
      bnlocal(7,2) = 'C'
      bnlocal(7,3) = 'N'
c     --- CH3E - C - NH1 [param19_eef1.inp]
      bnlocal(8,1) = 'CH3E'
      bnlocal(8,2) = 'C'
      bnlocal(8,3) = 'NH1'
c     --- CT - C - N [param19_eef1.inp]
      bnlocal(9,1) = 'CT'
      bnlocal(9,2) = 'C'
      bnlocal(9,3) = 'N'
c     --- CT - C - NH1 [param19_eef1.inp]
      bnlocal(10,1) = 'CT'
      bnlocal(10,2) = 'C'
      bnlocal(10,3) = 'NH1'
c     --- CT - C - NH2 [param19_eef1.inp]
      bnlocal(11,1) = 'CT'
      bnlocal(11,2) = 'C'
      bnlocal(11,3) = 'NH2'
c     --- CH2E - CR - NH1 [param19_eef1.inp]
      bnlocal(12,1) = 'CH2E'
      bnlocal(12,2) = 'CR'
      bnlocal(12,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- CH1E - C - O [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
c     --- CH3E - C - O [param19_eef1.inp]
      bnlocal(2,1) = 'CH3E'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'O'
c     --- CT - C - O [param19_eef.inp]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'O'
c     --- HA - C - O [param19_eef.inp]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH1E - C - OC [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- CH1E - C - OH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OH1'
c     --- CH2E - C - OH1 [param19_eef1.inp]
      bnlocal(2,1) = 'CH2E'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'OH1'
c     --- O - C - OH1 [param19_eef1.inp]
      bnlocal(3,1) = 'O'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'OH1'
c     --- CT - C - OH1 [param19_eef1.inp]
      bnlocal(4,1) = 'CT'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'OH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH1E - CH1E - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'CH2E'
c     --- CH1E - CH2E - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'CH2E'
      bnlocal(2,3) = 'CH2E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 5
c     --- CH1E - CH1E - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'CH3E'
c     --- CH1E - CH2E - OH1 [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'CH2E'
      bnlocal(2,3) = 'OH1'
c     --- CH2E - CH2E - CH3E [param19_eef1.inp]
      bnlocal(3,1) = 'CH2E'
      bnlocal(3,2) = 'CH2E'
      bnlocal(3,3) = 'CH3E'
c     --- CH3E - CH2E - OH1 [param19_eef1.inp]
      bnlocal(4,1) = 'CH3E'
      bnlocal(4,2) = 'CH2E'
      bnlocal(4,3) = 'OH1'
c     --- CT - CT - CT [param19_eef1.inp]
      bnlocal(5,1) = 'CT'
      bnlocal(5,2) = 'CT'
      bnlocal(5,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH1E - CH1E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH1E - CH1E - NH3 [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'NH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH1E - CH1E - OH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'OH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 117.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH1E - CH2E - CH1E [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'CH1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH1E - CH2E - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- CH1E - CH2E - S [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'S'
c     --- CH1E - CH2E - SH1E [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'CH2E'
      bnlocal(2,3) = 'SH1E'
c     --- CH2E - CH2E - S [param19_eef1.inp]
      bnlocal(3,1) = 'CH2E'
      bnlocal(3,2) = 'CH2E'
      bnlocal(3,3) = 'S'
c     --- CT - CT - S [param19_eef1.inp]
      bnlocal(4,1) = 'CT'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CH1E - N - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'N'
      bnlocal(1,3) = 'CH2E'
c     --- CH1E - N - CH3E [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'N'
      bnlocal(2,3) = 'CH3E'
c     --- CT - N - CT [param19_eef1.inp]
      bnlocal(3,1) = 'CT'
      bnlocal(3,2) = 'N'
      bnlocal(3,3) = 'CT'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CH1E - NH1 - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'NH1'
      bnlocal(1,3) = 'CH3E'
c     --- CH1E - NH2 - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'NH2'
      bnlocal(2,3) = 'CH2E'
c     --- CH2E - NH1 - CH3E [param19_eef1.inp]
      bnlocal(3,1) = 'CH2E'
      bnlocal(3,2) = 'NH1'
      bnlocal(3,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (35.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 7
c     --- CH1E - NH3 - HC [param19_eef1.inp]
      bnlocal(1,1) = 'CH1E'
      bnlocal(1,2) = 'NH3'
      bnlocal(1,3) = 'HC'
c     --- CH1E - NH3 - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'CH1E'
      bnlocal(2,2) = 'NH3'
      bnlocal(2,3) = 'CH2E'
c     --- CH1E - OH1 - H [param19_eef1.inp]
      bnlocal(3,1) = 'CH1E'
      bnlocal(3,2) = 'OH1'
      bnlocal(3,3) = 'H'
c     --- CH2E - OH1 - H [param19_eef1.inp]
      bnlocal(4,1) = 'CH2E'
      bnlocal(4,2) = 'OH1'
      bnlocal(4,3) = 'H'
c     --- CH2E - NH3 - HC [param19_eef1.inp]
      bnlocal(5,1) = 'CH2E'
      bnlocal(5,2) = 'NH3'
      bnlocal(5,3) = 'HC'
c     --- CT - NH3 - HC [param19_eef1.inp]
      bnlocal(6,1) = 'CT'
      bnlocal(6,2) = 'NH3'
      bnlocal(6,3) = 'HC'
c     --- CT - OH1 - H [param19_eef1.inp]
      bnlocal(7,1) = 'CT'
      bnlocal(7,2) = 'OH1'
      bnlocal(7,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.5d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH2E - C - CR1E [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CR1E'
c     --- CH2E - CR - CR1E [param19_eef1.inp]
      bnlocal(2,1) = 'CH2E'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 116.0d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH2E - C - NR [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
c     --- CH2E - CR - NR [param19_eef1.inp]
      bnlocal(2,1) = 'CH2E'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.6d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - C - O [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 118.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH2E - C - 0C [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'OC'
c     --- CT - C - 0C [param19_eef1.inp]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - CH1E - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - CH1E - N [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 9
c     --- CH2E - CH1E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'NH1'
c     --- CH2E - CH1E - NH2 [param19_eef1.inp]
      bnlocal(2,1) = 'CH2E'
      bnlocal(2,2) = 'CH1E'
      bnlocal(2,3) = 'NH2'
c     --- CH2E - CH1E - NH3 [param19_eef1.inp]
      bnlocal(3,1) = 'CH2E'
      bnlocal(3,2) = 'CH1E'
      bnlocal(3,3) = 'NH3'
c     --- CR1E - NH1 - CR1E [param19_eef1.inp]
      bnlocal(4,1) = 'CR1E'
      bnlocal(4,2) = 'NH1'
      bnlocal(4,3) = 'CR1E'
c     --- CR1E - NR - CR1E [param19_eef1.inp]
      bnlocal(5,1) = 'CR1E'
      bnlocal(5,2) = 'NR'
      bnlocal(5,3) = 'CR1E'
c     --- CT - CT - NC2 [param19_eef1.inp]
      bnlocal(6,1) = 'CT'
      bnlocal(6,2) = 'CT'
      bnlocal(6,3) = 'NC2'
c     --- CT - CT - NH1 [param19_eef1.inp]
      bnlocal(7,1) = 'CT'
      bnlocal(7,2) = 'CT'
      bnlocal(7,3) = 'NH1'
c     --- CT - CT - NH2 [param19_eef1.inp]
      bnlocal(8,1) = 'CT'
      bnlocal(8,2) = 'CT'
      bnlocal(8,3) = 'NH2'
c     --- CT - CT - NH3 [param19_eef1.inp]
      bnlocal(9,1) = 'CT'
      bnlocal(9,2) = 'CT'
      bnlocal(9,3) = 'NH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (45.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - CH2E - CH2E [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'CH2E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 105.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CH2E - CH2E - N [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'N'
c     --- CT - CT - N [param19_eef1.inp]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'N'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - CH2E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CH2E - CH2E - NH2 [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'NH2'
c     --- CH3E - CH1E - NH2 [param19_eef1.inp]
      bnlocal(2,1) = 'CH3E'
      bnlocal(2,2) = 'CH1E'
      bnlocal(2,3) = 'NH2'
c     --- CH3E - CH1E - NH3 [param19_eef1.inp]
      bnlocal(3,1) = 'CH3E'
      bnlocal(3,2) = 'CH1E'
      bnlocal(3,3) = 'NH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CH2E - CH2E - NH3 [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'CH2E'
      bnlocal(1,3) = 'NH3'
c     --- CR1E - C - NH1 [param19_eef1.inp]
      bnlocal(2,1) = 'CR1E'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NH1'
c     --- CR1E - CR - NH1 [param19_eef1.inp]
      bnlocal(3,1) = 'CR1E'
      bnlocal(3,2) = 'CR'
      bnlocal(3,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- CH2E - N - CH3E [param19_eef1.inp]
      isame = isame + 1
      bnlocal(isame,1) = 'CH2E'
      bnlocal(isame,2) = 'N'
      bnlocal(isame,3) = 'CH3E'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 99.5d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - S - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.2d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH2E - S - S [param19_eef1.inp]
      bnlocal(1,1) = 'CH2E'
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
c     --- CH3E - CH1E - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'CH3E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH3E - CH1E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH3E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'NH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.5d0
      bclocal(1) = (60.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CH3E - CH1E - OH1 [param19_eef1.inp]
      bnlocal(1,1) = 'CH3E'
      bnlocal(1,2) = 'CH1E'
      bnlocal(1,3) = 'OH1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- CR1E - C - CR1E [param19_eef1.inp]
      bnlocal(1,1) = 'CR1E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CR1E'
c     --- CR1E - CR - CR1E [param19_eef1.inp]
      bnlocal(2,1) = 'CR1E'
      bnlocal(2,2) = 'CR'
      bnlocal(2,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
c     --- CR1E - C - NR [param19_eef1.inp]
      bnlocal(1,1) = 'CR1E'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
c     --- CR1E - C - NP [param19_eef1.inp]
      bnlocal(2,1) = 'CR1E'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NP'
c     --- CR1E - CR - NR [param19_eef1.inp]
      bnlocal(3,1) = 'CR1E'
      bnlocal(3,2) = 'CR'
      bnlocal(3,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (90.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CR1E - CR1E - CR1E [param19_eef1.inp]
      bnlocal(1,1) = 'CR1E'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'CR1E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.0d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- H - NH2 - H [param19_eef1.inp]
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'NH2'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.0d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 6
c     --- HC - NC2 - HC [param19_eef1.inp]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NC2'
      bnlocal(1,3) = 'HC'
c     --- C - C - HA [param19_eef1.inp]
      bnlocal(2,1) = 'C'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'HA'
c     --- HA - C - NR [param19_eef1.inp]
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'NR'
c     --- HA - C - NH2 [param19_eef1.inp]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'C'
      bnlocal(4,3) = 'NH2'
c     --- HA - C - NH1 [param19_eef1.inp]
      bnlocal(5,1) = 'HA'
      bnlocal(5,2) = 'C'
      bnlocal(5,3) = 'NH1'
c     --- HC - NC2 - HC [param19_eef1.inp]
      bnlocal(6,1) = 'HC'
      bnlocal(6,2) = 'NC2'
      bnlocal(6,3) = 'HC'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (40.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 4
c     --- HC - NH3 - HC [param19_eef1.inp]
      bnlocal(1,1) = 'HC'
      bnlocal(1,2) = 'NH3'
      bnlocal(1,3) = 'HC'
c     --- CT - CT - HA [param19_eef1.inp]
      bnlocal(2,1) = 'CT'
      bnlocal(2,2) = 'CT'
      bnlocal(2,3) = 'HA'
c     --- HA - CT - HA [param19_eef1.inp]
      bnlocal(3,1) = 'HA'
      bnlocal(3,2) = 'CT'
      bnlocal(3,3) = 'HA'
c     --- HA - CT - S [param19_eef1.inp]
      bnlocal(4,1) = 'HA'
      bnlocal(4,2) = 'CT'
      bnlocal(4,3) = 'S'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 104.52d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- HT - OT - HT [param19_eef1.inp]
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
c     --- N - C - O [param19_eef1.inp]
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
      nbslocal = 6
c     --- NH1 - C - NR [param19_eef1.inp]
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'NR'
c     --- NC2 - C - NH1 [param19_eef1.inp]
      bnlocal(2,1) = 'NC2'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'NH1'
c     --- NC2 - C - NC2 [param19_eef1.inp]
      bnlocal(3,1) = 'NC2'
      bnlocal(3,2) = 'C'
      bnlocal(3,3) = 'NC2'
c     --- NC2 - CR - NH1 [param19_eef1.inp]
      bnlocal(4,1) = 'NC2'
      bnlocal(4,2) = 'CR'
      bnlocal(4,3) = 'NH1'
c     --- NC2 - CR - NC2 [param19_eef1.inp]
      bnlocal(5,1) = 'NC2'
      bnlocal(5,2) = 'CR'
      bnlocal(5,3) = 'NC2'
c     --- OS - C - O [param19_eef1.inp]
      bnlocal(6,1) = 'OS'
      bnlocal(6,2) = 'C'
      bnlocal(6,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.3d0
      bclocal(1) = (70.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- OS - C - CH1E [param19_eef1.inp]
      bnlocal(1,1) = 'OS'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'CH1E'
c     --- OS - C - CH2E [param19_eef1.inp]
      bnlocal(2,1) = 'OS'
      bnlocal(2,2) = 'C'
      bnlocal(2,3) = 'CH2E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 121.0d0
      bclocal(1) = (65.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- NH1 - C - O [param19_eef1.inp]
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'C'
      bnlocal(1,3) = 'O'
c     --- NH2 - C - O [param19_eef1.inp]
      bnlocal(2,1) = 'NH2'
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
c     --- NH1 - CR1E - NH1 [param19_eef1.inp]
      bnlocal(1,1) = 'NH1'
      bnlocal(1,2) = 'CR1E'
      bnlocal(1,3) = 'NH1'
c     --- NH1 - CR1E - NR [param19_eef1.inp]
      bnlocal(2,1) = 'NH1'
      bnlocal(2,2) = 'CR1E'
      bnlocal(2,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.5d0
      bclocal(1) = (85.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- OC - C - OC [param19_eef1.inp]
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
c     --- OC - S - OC [param19_eef1.inp]
      bnlocal(1,1) = 'OC'
      bnlocal(1,2) = 'S'
      bnlocal(1,3) = 'OC'
      call twh_set_bendnames(type,nbslocal,bnlocal)
      
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- FE - CM - OM [param19_eef1.inp]
      bnlocal(1,1) = 'FE'
      bnlocal(1,2) = 'CM'
      bnlocal(1,3) = 'OM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (0.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- FE - OM - OM [param19_eef1.inp]
      bnlocal(1,1) = 'FE'
      bnlocal(1,2) = 'OM'
      bnlocal(1,3) = 'OM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 128.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - NP - FE [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'NP'
      bnlocal(1,3) = 'FE'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 124.8d0
      bclocal(1) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- CR1E - NR - FE [param19_eef1.inp]
      bnlocal(1,1) = 'CR1E'
      bnlocal(1,2) = 'NR'
      bnlocal(1,3) = 'FE'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
c     --- NP - FE - CM [param19_eef1.inp]
      bnlocal(1,1) = 'NP'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'CM'
c     --- NP - FE - OM [param19_eef1.inp]
      bnlocal(2,1) = 'NP'
      bnlocal(2,2) = 'FE'
      bnlocal(2,3) = 'OM'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NP - FE - NP [param19_eef1.inp]
      bnlocal(1,1) = 'NP'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'NP'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 115.0d0
      bclocal(1) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- NP - FE - NR [param19_eef1.inp]
      bnlocal(1,1) = 'NP'
      bnlocal(1,2) = 'FE'
      bnlocal(1,3) = 'NR'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 120.5d0
      bclocal(1) = (46.50d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- C - OS - CH3E [param19_eef1.inp]
      bnlocal(1,1) = 'C'
      bnlocal(1,2) = 'OS'
      bnlocal(1,3) = 'CH3E'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47122d0
      bclocal(1) = (55.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- H - OH2 - H [param19_eef1.inp]
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'OH2'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * C19eef1 - the original version from the Neria paper           *
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
c     --- CH1E - C - N - CH1E [param19_eef1.inp]
      tnlocal(1,1) = 'CH1E'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'CH1E'
c     --- CH2E - C - N - CH1E [param19_eef1.inp]
      tnlocal(2,1) = 'CH2E'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'N'
      tnlocal(2,4) = 'CH1E'
c     --- X - C - CR1E - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'CR1E'
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
      ntslocal = 4
c     --- CR1E - C - C - C [param19_eef1.inp]
      tnlocal(1,1) = 'CR1E'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'C'
c     --- CR1E - C - C - NH1 [param19_eef1.inp]
      tnlocal(2,1) = 'CR1E'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'C'
      tnlocal(2,4) = 'NH1'
c     --- CR1E - CR - CR - CR [param19_eef1.inp]
      tnlocal(3,1) = 'CR1E'
      tnlocal(3,2) = 'CR'
      tnlocal(3,3) = 'CR'
      tnlocal(3,4) = 'CR'
c     --- CR1E - CR - CR - NH1 [param19_eef1.inp]
      tnlocal(4,1) = 'CR1E'
      tnlocal(4,2) = 'CR'
      tnlocal(4,3) = 'CR'
      tnlocal(4,4) = 'NH1'
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
      ntslocal = 2
c     --- CR1E - C - C - CR1E [param19_eef1.inp]
      tnlocal(1,1) = 'CR1E'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'CR1E'
c     --- CR1E - CR - CR - CR1E [param19_eef1.inp]
      tnlocal(2,1) = 'CR1E'
      tnlocal(2,2) = 'CR'
      tnlocal(2,3) = 'CR'
      tnlocal(2,4) = 'CR1E'
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
      ntslocal = 5
c     --- X - C - CH1E - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'CH1E'
      tnlocal(1,4) = 'wild'
c     --- X - C - CH2E - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'CH2E'
      tnlocal(2,4) = 'wild'
c     --- X - CR - CH2E - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CR'
      tnlocal(3,3) = 'CH2E'
      tnlocal(3,4) = 'wild'
c     --- X - C - CT - X [param19_eef1.inp]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'wild'
c     --- null - null - null - null [MGM hack]
c     --- this is used for the torsions that are not specified by
c     --- charmm19, but need to be listed for the logic of Towhee
c     --- to function properly
      tnlocal(5,1) = 'null'
      tnlocal(5,2) = 'null'
      tnlocal(5,3) = 'null'
      tnlocal(5,4) = 'null'
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
      ntslocal = 6
c     --- X - C - N - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'wild'
c     --- X - C - NH1 - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'NH1'
      tnlocal(2,4) = 'wild'
c     --- X - C - NH2 - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'C'
      tnlocal(3,3) = 'NH2'
      tnlocal(3,4) = 'wild'
c     --- X - C - NC2 - X [param19_eef1.inp]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'C'
      tnlocal(4,3) = 'NC2'
      tnlocal(4,4) = 'wild'
c     --- X - CR - NC2 - X [param19_eef1.inp]
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CR'
      tnlocal(5,3) = 'NC2'
      tnlocal(5,4) = 'wild'
c     --- X - CR - NH1 - X [param19_eef1.inp]
      tnlocal(6,1) = 'wild'
      tnlocal(6,2) = 'CR'
      tnlocal(6,3) = 'NH1'
      tnlocal(6,4) = 'wild'
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
      ntslocal = 2
c     --- X - C - OH1 - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'OH1'
      tnlocal(1,4) = 'wild'
c     --- X - C - OS - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'C'
      tnlocal(2,3) = 'OS'
      tnlocal(2,4) = 'wild'
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
      ntslocal = 4
c     --- X - CH1E - CH1E - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH1E'
      tnlocal(1,3) = 'CH1E'
      tnlocal(1,4) = 'wild'
c     --- X - CH1E - CH2E - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CH1E'
      tnlocal(2,3) = 'CH2E'
      tnlocal(2,4) = 'wild'
c     --- X - CH2E - CH2E - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CH2E'
      tnlocal(3,3) = 'CH2E'
      tnlocal(3,4) = 'wild'
c     --- X - CT - CT - X [param19_eef1.inp]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'CT'
      tnlocal(4,4) = 'wild'
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
      ntslocal = 7
c     --- X - CH1E - N - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH1E'
      tnlocal(1,3) = 'N'
      tnlocal(1,4) = 'wild'
c     --- X - CH1E - NH1 - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CH1E'
      tnlocal(2,3) = 'NH1'
      tnlocal(2,4) = 'wild'
c     --- X - CH2E - N - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CH2E'
      tnlocal(3,3) = 'N'
      tnlocal(3,4) = 'wild'
c     --- X - CH2E - NH1 - X [param19_eef1.inp]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CH2E'
      tnlocal(4,3) = 'NH1'
      tnlocal(4,4) = 'wild'
c     --- X - CT - N - X [param19_eef1.inp]
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'N'
      tnlocal(5,4) = 'wild'
c     --- X - CT - NC2 - X [param19_eef1.inp]
      tnlocal(6,1) = 'wild'
      tnlocal(6,2) = 'CT'
      tnlocal(6,3) = 'NC2'
      tnlocal(6,4) = 'wild'
c     --- X - CT - NH1 - X [param19_eef1.inp]
      tnlocal(7,1) = 'wild'
      tnlocal(7,2) = 'CT'
      tnlocal(7,3) = 'NH1'
      tnlocal(7,4) = 'wild'
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
c     --- X - CH1E - NH2 - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH1E'
      tnlocal(1,3) = 'NH2'
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
      ntslocal = 5
c     --- X - CH1E - NH3 - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH1E'
      tnlocal(1,3) = 'NH3'
      tnlocal(1,4) = 'wild'
c     --- X - CH2E - NH2 - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CH2E'
      tnlocal(2,3) = 'NH2'
      tnlocal(2,4) = 'wild'
c     --- X - CH2E - NH3 - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CH2E'
      tnlocal(3,3) = 'NH3'
      tnlocal(3,4) = 'wild'
c     --- X - CT - NH2 - X [param19_eef1.inp]
      tnlocal(4,1) = 'wild'
      tnlocal(4,2) = 'CT'
      tnlocal(4,3) = 'NH2'
      tnlocal(4,4) = 'wild'
c     --- X - CT - NH3 - X [param19_eef1.inp]
      tnlocal(5,1) = 'wild'
      tnlocal(5,2) = 'CT'
      tnlocal(5,3) = 'NH3'
      tnlocal(5,4) = 'wild'
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
      ntslocal = 3
c     --- X - CH1E - OH1 - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH1E'
      tnlocal(1,3) = 'OH1'
      tnlocal(1,4) = 'wild'
c     --- X - CH2E - OH1 - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CH2E'
      tnlocal(2,3) = 'OH1'
      tnlocal(2,4) = 'wild'
c     --- X - CT - OH1 - X [param19_eef1.inp]
      tnlocal(3,1) = 'wild'
      tnlocal(3,2) = 'CT'
      tnlocal(3,3) = 'OH1'
      tnlocal(3,4) = 'wild'
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
      ntslocal = 2
c     --- X - CH2E - S - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CH2E'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'wild'
c     --- X - CT - S - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CT'
      tnlocal(2,3) = 'S'
      tnlocal(2,4) = 'wild'
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
c     --- X - S - S - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'S'
      tnlocal(1,3) = 'S'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.050d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- X - FE - NR - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'FE'
      tnlocal(1,3) = 'NR'
      tnlocal(1,4) = 'wild'
c     --- X - FE - CM - X [param19_eef1.inp]
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'FE'
      tnlocal(2,3) = 'CM'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.4d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.000d0)*(CON_KCALTOK)
      tclocal(2) = 4.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- X - FE - NR - X [param19_eef1.inp]
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'FE'
      tnlocal(1,3) = 'OM'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- below is c19eef1 torsions that are listed in charmm as 
c     --- impropers.  I have modified these from param19_eef1.inp so 
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
c     --- CR - NH1 - CR1E - NH1 
c     --- (from CR - X - X - NH1) [param19_eef1.inp]
      tnlocal(1,1) = 'CR'
      tnlocal(1,2) = 'NH1'
      tnlocal(1,3) = 'CR1E'
      tnlocal(1,4) = 'NH1'
c     --- CR - NR - CR1E - NH1 
c     --- (from CR - X - X - NH1) [param19_eef1.inp]
      tnlocal(2,1) = 'CR'
      tnlocal(2,2) = 'NR'
      tnlocal(2,3) = 'CR1E'
      tnlocal(2,4) = 'NH1'
c     --- CR - CR - CR - NH1 
c     --- (from CR - X - X - NH1) [param19_eef1.inp]
      tnlocal(3,1) = 'CR'
      tnlocal(3,2) = 'CR'
      tnlocal(3,3) = 'CR'
      tnlocal(3,4) = 'NH1'
c     --- CR - CR - CR1E - NH1 
c     --- (from CR - X - X - NH1) [param19_eef1.inp]
      tnlocal(4,1) = 'CR'
      tnlocal(4,2) = 'CR'
      tnlocal(4,3) = 'CR1E'
      tnlocal(4,4) = 'NH1'
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
c     --- NH1 - CR1E - NH1 - CR1E 
c     --- (from CR1E - X - X - NH1) [param19_eef1.inp]
      tnlocal(1,1) = 'NH1'
      tnlocal(1,2) = 'CR1E'
      tnlocal(1,3) = 'NH1'
      tnlocal(1,4) = 'CR1E'
c     --- CR1E - NH1 - CR1E - CR 
c     --- (from CR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(2,1) = 'CR1E'
      tnlocal(2,2) = 'NH1'
      tnlocal(2,3) = 'CR1E'
      tnlocal(2,4) = 'CR'
c     --- NH1 - CR1E - CR - NH1 
c     --- (from NH1 - X - X - NH1) [param19_eef1.inp]
      tnlocal(3,1) = 'NH1'
      tnlocal(3,2) = 'CR1E'
      tnlocal(3,3) = 'CR'
      tnlocal(3,4) = 'NH1'
c     --- CR1E - CR - NH1 - CR1E 
c     --- (from CR1E - X - X - CR1E) [param19_eef1.inp]
      tnlocal(4,1) = 'CR1E'
      tnlocal(4,2) = 'CR'
      tnlocal(4,3) = 'NH1'
      tnlocal(4,4) = 'CR1E'
c     --- CR - NH1 - CR1E - NR 
c     --- (from NR - X - X - CR) [param19_eef1.inp]
      tnlocal(5,1) = 'CR'
      tnlocal(5,2) = 'NH1'
      tnlocal(5,3) = 'CR1E'
      tnlocal(5,4) = 'NR'
c     --- NH1 - CR1E - NR - CR1E 
c     --- (from CR1E - X - X - NH1) [param19_eef1.inp]
      tnlocal(6,1) = 'NH1'
      tnlocal(6,2) = 'CR1E'
      tnlocal(6,3) = 'NR'
      tnlocal(6,4) = 'CR1E'
c     --- CR1E - NR - CR1E - CR 
c     --- (from CR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(7,1) = 'CR1E'
      tnlocal(7,2) = 'NR'
      tnlocal(7,3) = 'CR1E'
      tnlocal(7,4) = 'CR'
c     --- NR - CR1E - CR - NH1 
c     --- (from NH1 - X - X - NR) [param19_eef1.inp]
      tnlocal(8,1) = 'NR'
      tnlocal(8,2) = 'CR1E'
      tnlocal(8,3) = 'CR'
      tnlocal(8,4) = 'NH1'
c     --- NR - CR1E - NH1 - CR1E 
c     --- (from NR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(9,1) = 'NR'
      tnlocal(9,2) = 'CR1E'
      tnlocal(9,3) = 'NH1'
      tnlocal(9,4) = 'CR1E'
c     --- NH1 - CR1E - CR - NR 
c     --- (from NH1 - X - X - NR) [param19_eef1.inp]
      tnlocal(10,1) = 'NH1'
      tnlocal(10,2) = 'CR1E'
      tnlocal(10,3) = 'CR'
      tnlocal(10,4) = 'NR'
c     --- CR1E - CR - NR - CR1E 
c     --- (from CR1E - X - X - CR1E) [param19_eef1.inp]
      tnlocal(11,1) = 'CR1E'
      tnlocal(11,2) = 'CR'
      tnlocal(11,3) = 'NR'
      tnlocal(11,4) = 'CR1E'
c     --- CR - CR1E - CR1E - CR1E 
c     --- (from CR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(12,1) = 'CR'
      tnlocal(12,2) = 'CR1E'
      tnlocal(12,3) = 'CR1E'
      tnlocal(12,4) = 'CR1E'
c     --- CR1E - CR1E - CR1E - CR1E 
c     --- (from CR1E - X - X - CR1E) [param19_eef1.inp]
      tnlocal(13,1) = 'CR1E'
      tnlocal(13,2) = 'CR1E'
      tnlocal(13,3) = 'CR1E'
      tnlocal(13,4) = 'CR1E'
c     --- CR1E - CR - CR1E - CR1E 
c     --- (from CR1E - X - X - CR1E) [param19_eef1.inp]
      tnlocal(14,1) = 'CR1E'
      tnlocal(14,2) = 'CR'
      tnlocal(14,3) = 'CR1E'
      tnlocal(14,4) = 'CR1E'
c     --- CR - CR1E - NH1 - CR 
c     --- (from CR - X - X - CR) [param19_eef1.inp]
      tnlocal(15,1) = 'CR'
      tnlocal(15,2) = 'CR1E'
      tnlocal(15,3) = 'NH1'
      tnlocal(15,4) = 'CR'
c     --- CR - CR - CR - CR1E 
c     --- (from CR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(16,1) = 'CR'
      tnlocal(16,2) = 'CR'
      tnlocal(16,3) = 'CR'
      tnlocal(16,4) = 'CR1E'
c     --- CR1E - CR1E - CR - CR 
c     --- (from CR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(17,1) = 'CR1E'
      tnlocal(17,2) = 'CR1E'
      tnlocal(17,3) = 'CR'
      tnlocal(17,4) = 'CR'
c     --- CR1E - CR - CR - CR1E 
c     --- (from CR1E - X - X - CR1E) [param19_eef1.inp]
      tnlocal(18,1) = 'CR1E'
      tnlocal(18,2) = 'CR'
      tnlocal(18,3) = 'CR'
      tnlocal(18,4) = 'CR1E'
c     --- C - CR1E - CR1E - CR 
c     --- (from CR - X - X - C) [param19_eef1.inp]
      tnlocal(19,1) = 'C'
      tnlocal(19,2) = 'CR1E'
      tnlocal(19,3) = 'CR1E'
      tnlocal(19,4) = 'CR'
c     --- CR1E - NH1 - CR - CR 
c     --- (from CR - X - X - CR1E) [param19_eef1.inp]
      tnlocal(20,1) = 'CR1E'
      tnlocal(20,2) = 'NH1'
      tnlocal(20,3) = 'CR'
      tnlocal(20,4) = 'CR'
c     --- CR1E - C - CR1E - CR1E 
c     --- (from CR1E - X - X - CR1E) [param19_eef1.inp]
      tnlocal(21,1) = 'CR1E'
      tnlocal(21,2) = 'C'
      tnlocal(21,3) = 'CR1E'
      tnlocal(21,4) = 'CR1E'
      call twh_set_tornames(type,ntslocal,tnlocal)
      
      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * C19eef1 improper harmonic torsion                             *
c     ******************************************************************

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (20.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 1
c     --- FE - X - X - NP [param19_eef1.inp]
      ninlocal(1,1) = 'FE'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'NP'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (90.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 9
c     --- C - C - CR1E - CH2E [param19_eef1.inp]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'C'
      ninlocal(1,3) = 'CR1E'
      ninlocal(1,4) = 'CH2E'
c     --- C - CR1E - C - CH2E [param19_eef1.inp]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'CR1E'
      ninlocal(2,3) = 'C'
      ninlocal(2,4) = 'CH2E'
c     --- C - CR1E - CR1E - CH2E [param19_eef1.inp]
      ninlocal(3,1) = 'C'
      ninlocal(3,2) = 'CR1E'
      ninlocal(3,3) = 'CR1E'
      ninlocal(3,4) = 'CH2E'
c     --- C - CR1E - NH1 - CH2E [param19_eef1.inp]
      ninlocal(4,1) = 'C'
      ninlocal(4,2) = 'CR1E'
      ninlocal(4,3) = 'NH1'
      ninlocal(4,4) = 'CH2E'
c     --- C - NH1 - CR1E - CH2E [param19_eef1.inp]
      ninlocal(5,1) = 'C'
      ninlocal(5,2) = 'NH1'
      ninlocal(5,3) = 'CR1E'
      ninlocal(5,4) = 'CH2E'
c     --- C - X - X - CH2E [param19_eef1.inp]
      ninlocal(6,1) = 'C'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'CH2E'
c     --- C - X - X - CH3E [param19_eef1.inp]
      ninlocal(7,1) = 'C'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'wild'
      ninlocal(7,4) = 'CH3E'
c     --- CR - CR1E - CR - CH2E [param19_eef1.inp]
      ninlocal(8,1) = 'CR'
      ninlocal(8,2) = 'CR1E'
      ninlocal(8,3) = 'CR'
      ninlocal(8,4) = 'CH2E'
c     --- CR - X - X - CH2E [param19_eef1.inp]
      ninlocal(9,1) = 'CR'
      ninlocal(9,2) = 'wild'
      ninlocal(9,3) = 'wild'
      ninlocal(9,4) = 'CH2E'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (150.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- C - CR1E - CR1E - OH1 [param19_eef1.inp]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'CR1E'
      ninlocal(1,3) = 'CR1E'
      ninlocal(1,4) = 'OH1'
c     --- C - X - X - OH1 [param19_eef1.inp]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'OH1'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (45.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 13
c     --- C - H - H - NH2 [param19_eef1.inp]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'H'
      ninlocal(1,3) = 'H'
      ninlocal(1,4) = 'NH2'
c     --- H - X - X - O [param19_eef1.inp]
      ninlocal(2,1) = 'H'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'O'
c     --- N - CH1E - CH2E - C [param19_eef1.inp]
      ninlocal(3,1) = 'N'
      ninlocal(3,2) = 'CH1E'
      ninlocal(3,3) = 'CH2E'
      ninlocal(3,4) = 'C'
c     --- N - X - X - CH2E [param19_eef1.inp]
      ninlocal(4,1) = 'N'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'CH2E'
c     --- NH1 - X - X - CH1E [param19_eef1.inp]
      ninlocal(5,1) = 'NH1'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'CH1E'
c     --- NH1 - X - X - CH2E [param19_eef1.inp]
      ninlocal(6,1) = 'NH1'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'CH2E'
c     --- NH1 - X - X - CH3E [param19_eef1.inp]
      ninlocal(7,1) = 'NH1'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'wild'
      ninlocal(7,4) = 'CH3E'
c     --- NH1 - X - X - H [param19_eef1.inp]
      ninlocal(8,1) = 'NH1'
      ninlocal(8,2) = 'wild'
      ninlocal(8,3) = 'wild'
      ninlocal(8,4) = 'H'
c     --- NH2 - X - X - H [param19_eef1.inp]
      ninlocal(9,1) = 'NH2'
      ninlocal(9,2) = 'wild'
      ninlocal(9,3) = 'wild'
      ninlocal(9,4) = 'H'
c     --- NC2 - X - X - HC [param19_eef1.inp]
      ninlocal(10,1) = 'NC2'
      ninlocal(10,2) = 'wild'
      ninlocal(10,3) = 'wild'
      ninlocal(10,4) = 'HC'
c     --- N - X - X - CT [param19_eef1.inp]
      ninlocal(11,1) = 'N'
      ninlocal(11,2) = 'wild'
      ninlocal(11,3) = 'wild'
      ninlocal(11,4) = 'CT'
c     --- NC2 - X - X - CT [param19_eef1.inp]
      ninlocal(12,1) = 'NC2'
      ninlocal(12,2) = 'wild'
      ninlocal(12,3) = 'wild'
      ninlocal(12,4) = 'CT'
c     --- NH1 - X - X - CT [param19_eef1.inp]
      ninlocal(13,1) = 'NH1'
      ninlocal(13,2) = 'wild'
      ninlocal(13,3) = 'wild'
      ninlocal(13,4) = 'CT'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (100.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 6
c     --- C - OC - OC - CH1E [param19_eef1.inp]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'OC'
      ninlocal(1,3) = 'OC'
      ninlocal(1,4) = 'CH1E'
c     --- C - OC - OC - CH2E [param19_eef1.inp]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'OC'
      ninlocal(2,3) = 'OC'
      ninlocal(2,4) = 'CH2E'
c     --- C - X - X - NH1 [param19_eef1.inp]
      ninlocal(3,1) = 'C'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'NH1'
c     --- C - X - X - O [param19_eef1.inp]
      ninlocal(4,1) = 'C'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'O'
c     --- C - X - X - OC [param19_eef1.inp]
      ninlocal(5,1) = 'C'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'OC'
c     --- CR - X - X - NH1 [param19_eef1.inp]
      ninlocal(6,1) = 'CR'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'NH1'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (25.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 14
c     --- C - X - X - C [param19_eef1.inp]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'C'
c     --- C - X - X - CR1E [param19_eef1.inp]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'CR1E'
c     --- CR1E - X - X - CR1E [param19_eef1.inp]
      ninlocal(3,1) = 'CR1E'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'CR1E'
c     --- CR1E - X - X - NH1 [param19_eef1.inp]
      ninlocal(4,1) = 'CR1E'
      ninlocal(4,2) = 'wild'
      ninlocal(4,3) = 'wild'
      ninlocal(4,4) = 'NH1'
c     --- NH1 - X - X - NH1 [param19_eef1.inp]
      ninlocal(5,1) = 'NH1'
      ninlocal(5,2) = 'wild'
      ninlocal(5,3) = 'wild'
      ninlocal(5,4) = 'NH1'
c     --- NH1 - X - X - NR [param19_eef1.inp]
      ninlocal(6,1) = 'NH1'
      ninlocal(6,2) = 'wild'
      ninlocal(6,3) = 'wild'
      ninlocal(6,4) = 'NR'
c     --- NR - X - X - C [param19_eef1.inp]
      ninlocal(7,1) = 'NR'
      ninlocal(7,2) = 'wild'
      ninlocal(7,3) = 'wild'
      ninlocal(7,4) = 'C'
c     --- NR - X - X - CR1E [param19_eef1.inp]
      ninlocal(8,1) = 'NR'
      ninlocal(8,2) = 'wild'
      ninlocal(8,3) = 'wild'
      ninlocal(8,4) = 'CR1E'
c     --- NR - X - X - CH3E [param19_eef1.inp]
      ninlocal(9,1) = 'NR'
      ninlocal(9,2) = 'wild'
      ninlocal(9,3) = 'wild'
      ninlocal(9,4) = 'CH3E'
c     --- CR - X - X - CR [param19_eef1.inp]
      ninlocal(10,1) = 'CR'
      ninlocal(10,2) = 'wild'
      ninlocal(10,3) = 'wild'
      ninlocal(10,4) = 'CR'
c     --- CR - X - X - C [param19_eef1.inp]
      ninlocal(11,1) = 'CR'
      ninlocal(11,2) = 'wild'
      ninlocal(11,3) = 'wild'
      ninlocal(11,4) = 'C'
c     --- CR - X - X - CR1E [param19_eef1.inp]
      ninlocal(12,1) = 'CR'
      ninlocal(12,2) = 'wild'
      ninlocal(12,3) = 'wild'
      ninlocal(12,4) = 'CR1E'
c     --- NR - X - X - CR [param19_eef1.inp]
      ninlocal(13,1) = 'NR'
      ninlocal(13,2) = 'wild'
      ninlocal(13,3) = 'wild'
      ninlocal(13,4) = 'CR'
c     --- NR - X - X - CT [param19_eef1.inp]
      ninlocal(14,1) = 'NR'
      ninlocal(14,2) = 'wild'
      ninlocal(14,3) = 'wild'
      ninlocal(14,4) = 'CT'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (75.00d0)*(CON_KCALTOK)
      iclocal(1) = (0.0d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 2
c     --- C - X - X - H [param19_eef1.inp]
      ninlocal(1,1) = 'C'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'H'
c     --- C - X - X - HA [param19_eef1.inp]
      ninlocal(2,1) = 'C'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'HA'
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (55.00d0)*(CON_KCALTOK)
      iclocal(1) = (35.26439d0)*(onepi/180.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nislocal = 3
c     --- CH1E - X - X - CH1E [param19_eef1.inp]
      ninlocal(1,1) = 'CH1E'
      ninlocal(1,2) = 'wild'
      ninlocal(1,3) = 'wild'
      ninlocal(1,4) = 'CH1E'
c     --- CH1E - X - X - CH2E [param19_eef1.inp]
      ninlocal(2,1) = 'CH1E'
      ninlocal(2,2) = 'wild'
      ninlocal(2,3) = 'wild'
      ninlocal(2,4) = 'CH2E'
c     --- CH1E - X - X - CH3E [param19_eef1.inp]
      ninlocal(3,1) = 'CH1E'
      ninlocal(3,2) = 'wild'
      ninlocal(3,3) = 'wild'
      ninlocal(3,4) = 'CH3E'
      call twh_set_impnames(type,nislocal,ninlocal)
      
      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
     
      return
      end
