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
      subroutine twh_ffgordon(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters to test the          *
c     * gordon n-6 forcefield                                          *
c     *                                                                *
c     * originally written 11-30-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer ii,type,isame
      character*(FFNAMELEN) ffname
c     --- nonbond
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      character*(FFNAMELEN) vnlocal
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
c     --- bond
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
c     --- angle
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
c     --- torsion
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)


      ffname = 'Gordon'
      call twh_classical_potential(GLB_SET,'Gordon n-6')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      type = 0

c     --- Reference
c     P.A. Gordon "Development of intermolecular potentials for
c     predicting transport properties of hydrocarbons
c     J. Chem. Phys. vol. 125 p. 014504 (2006)

c     --- (CH4n9) from Table 1
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = 16.0426d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.785d0
      nbclocal(2) = (0.232d0)*CON_KCALTOK
      nbclocal(3) = 9.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH4n9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH4n10) from Table 1
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = 16.0426d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.778d0
      nbclocal(2) = (0.256d0)*CON_KCALTOK
      nbclocal(3) = 10.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH4n10'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH4n11) from Table 1
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = 16.0426d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.736d0
      nbclocal(2) = (0.275d0)*CON_KCALTOK
      nbclocal(3) = 11.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH4n11'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH4n12) from Table 1
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = 16.0426d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.730d0
c      nbclocal(2) = (0.29389d0)*CON_KCALTOK
c     --- using the TraPPE value instead of the converted one
      nbclocal(2) = 148.0d0
      nbclocal(3) = 12.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH4n12'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH4n18) from Table 1
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = 16.0426d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.706d0
      nbclocal(2) = (0.374d0)*CON_KCALTOK
      nbclocal(3) = 18.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH4n18'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH4n40) from Table 1
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = 16.0426d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.690d0
      nbclocal(2) = (0.491d0)*CON_KCALTOK
      nbclocal(3) = 40.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH4n40'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH3n14) from Table 3
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = (12.011d0)+(3.0d0)*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.680d0
      nbclocal(2) = (0.220d0)*CON_KCALTOK
      nbclocal(3) = 14.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH3n14'
      cbnlocal = 'Csp3'
      canlocal = 'CH3sp3'
      ctnlocal = 'CH3sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- (CH2n20) from Table 3
      type = type + 1
      ii = twh_getnbtype(type,type)
      call twh_nbff(GLB_SET,type,ffname)
      masslocal = (12.011d0)+(2.0d0)*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      nbclocal(1) = 3.8750d0
      nbclocal(2) = (0.147d0)*CON_KCALTOK
      nbclocal(3) = 20.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      nbnlocal ='CH2n20'
      cbnlocal = 'Csp3'
      canlocal = 'CH2sp3'
      ctnlocal = 'CH2sp3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     --- uses TraPPE-UA intramolecular terms
c     TraPPE fixed bond length [TraPPE-1]
c     [CHx](sp3) -- [CHy](sp3,sp2,aro) single bonds
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
c     --- set nvsame
      call twh_set_vibnames(type,isame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     --- harmonic angles from TraPPE-UA
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
      call twh_set_bendnames(type,isame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     --- TraPPE-UA dihedral
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
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
