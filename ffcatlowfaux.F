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
      subroutine twh_ffcatlowfaux(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the CatlowFaux   *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer j,ij,type
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

      ffname = 'CatlowFaux'
      call twh_classical_potential(GLB_SET,'Exponential-12-6')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     --- combination of 12-6 Lennard-Jones and Buckingham
c     --- NOTE: must explicitly set up all cross interactions
c     --- the proper format for this is to set up all interactions of
c     --- the lower index atom with all of the higher index atoms.
      
c     --- D.A. Faux, W. Smith, T.R. Forester  
c     --- J. Phys. Chem. B (1997) 101, 1762-1768.
c     --- use the following rules to convert from their notation 
c     --- nbclocal(1) = -C
c     --- nbclocal(2) = B
c     --- nbclocal(3) = A
c     --- nbclocal(4) = -1/rho
      
c     --- this strange combination of LJ and exponential was first
c     --- used (I think) by Catlow's group and the parameters below 
c     --- come from a variety of places.  If there is no mention of 
c     --- where the parameters are from then I took them out of 
c     --- the Faux paper (listed above) but that does not mean the 
c     --- Faux paper is the original reference for each of the cross
c     --- terms as Faux leveraged the work that Catlow had already 
c     --- done
      
c     --- (OZ) the oxygen in zeolite 4A
      type = 1
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OZ'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of OZ with OZ (1-1)
      j = 1
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(21.250d0)*CON_EVTOK
      nbclocal(3) = (22764.000d0)*CON_EVTOK
      nbclocal(4) = -(1.0d0/0.1490d0)
      nbclocal(2) = 0.0d0*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with Si (1-2)
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(50.000d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (1150.000d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with Al (1-3)
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.0d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (1300.000d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with Na (1-4)
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.0d0)*CON_EVTOK
      nbclocal(3) = (5836.840d0)*CON_EVTOK
      nbclocal(4) = -(1.0d0/0.2387d0)
      nbclocal(2) = (0.0d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with OW (1-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(23.429d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (5645.172d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with HW (1-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(5.413d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (660.778d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with Cc (1-7)
c     --- N. Raj; G. Sastre; C.R.A. Catlow J. Phys. Chem. B
c     --- 103, 11007-11015 (1999)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(17.654d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (1.1d4)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of OZ with Hc (1-8)
c     --- N. Raj; G. Sastre; C.R.A. Catlow J. Phys. Chem. B
c     --- 103, 11007-11015 (1999)
      j = 8
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(5.5717d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (1.5564d3)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
     
c     --- (Si) Silicon in zeolite 4a
      type = 2
      masslocal = 28.086d0
      elelocal = 'Si'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Si'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of Si with OW (2-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.586d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (10.644d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with HW (2-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.089d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (0.552d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Al) Aluminum in zeolite 4a
      type = 3
      masslocal = 26.982d0
      elelocal = 'Al'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Al'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of Al with OW (3-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.500d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (11.758d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Al with HW (3-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.075d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (0.626d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Na) Sodium in zeolite 4a
      type = 4
      masslocal = 22.990d0
      elelocal = 'Na'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of Na with OW (4-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.0d0)*CON_EVTOK
      nbclocal(3) = (5836.840d0)*CON_EVTOK
      nbclocal(4) = -(1.0d0/0.2387d0)
      nbclocal(2) = (0.0d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Na with Cc (4-7)
      j = 7
      ij = twh_getnbtype(type,j)
c     --- extrapolated using the OPLS-AA Sodium-2, the Cc-Cc from
c     --- Raj, and geometric mixing rules
c     --- M.G. Martin 03-07-2001
      nbclocal(1) = -(3.44d0)*CON_EVTOK
      nbclocal(3) = (0.00d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (2095.7816d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Na with Hc (4-8)
      j = 8
      ij = twh_getnbtype(type,j)
c     --- extrapolated using the OPLS-AA Sodium-2, the Hc-Hc from
c     --- Raj, and geometric mixing rules
c     --- M.G. Martin 03-06-2001
      nbclocal(1) = -(1.218d0)*CON_EVTOK
      nbclocal(3) = (0.00d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (289.59d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (OW) Oxygen in SPC-E water
c     --- H.J.C. Berendsen; J.R. Grigera; T.P. Straatsma;
c     --- "The Missing Term in Effective Pair Potentials"
c     --- J. Phys. Chem. 91, 6269-6271 (1987).
c     --- plugged this into the Faux/Catlow equations as a water 
c     --- potential
      type = 5
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='OW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of OW with OW (5-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -((0.37122*10.0d0)**6)*CON_KJTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = ((0.3428d0*10.0d0)**12)*CON_KJTOK
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (HW) Hydrogen in SPC-E water
c     --- H.J.C. Berendsen; J.R. Grigera; T.P. Straatsma;
c     --- "The Missing Term in Effective Pair Potentials"
c     --- J. Phys. Chem. 91, 6269-6271 (1987).
c     --- plugged this into the Faux/Catlow equations as a water 
c     --- potential
      type = 6
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of HW with HW (6-6)
c     --- note: HW has no nonbond interaction with itself
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.0d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (0.0d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Cc) Carbon in an alkane
      type = 7
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cc'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of Cc with Cc (7-7)
c     --- from N. Raj, G. Sastre, C.R.A. Catlow
c     --- J. Phys Chem B 1999, 103, 11007-11015.
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(18.0933d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (19692.0d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Cc with Hc (7-8)
c     --- from N. Raj, G. Sastre, C.R.A. Catlow
c     --- J. Phys Chem B 1999, 103, 11007-11015.
      j = 8
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(5.8415d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (2800.0d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Cc with HZ (7-9)
c     --- J.O. Titiloye, S.C. Parker, F.S. Stone, C.R.A. Catlow
c     --- J. Phys Chem 1991 95, 4038-4044.
      j = 9
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(5.844d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (2854.571d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Hc) Hydrogen bonded to a carbon in an alkane
      type = 8
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Hc'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of Hc with Hc (8-8)
c     --- from N. Raj, G. Sastre, C.R.A. Catlow
c     --- J. Phys Chem B 1999, 103, 11007-11015.
      j = 8
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(1.9867d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (384.84d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Hc with HZ (8-9)
c     --- J.O. Titiloye, S.C. Parker, F.S. Stone, C.R.A. Catlow
c     --- J. Phys Chem 1991 95, 4038-4044.
      j = 9
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.110d0)*CON_EVTOK
      nbclocal(3) = (0.0d0)*CON_EVTOK
      nbclocal(4) = -(0.0d0)
      nbclocal(2) = (3.553d0)*CON_EVTOK
      call twh_set_nbcoeff(ij,nbclocal) 

c     --- (HZ) Hydrogen in an OH group of a zeolite 
      type = 9
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HZ'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Fixed bond length for SPC-E water and simple chain models      *
c     ******************************************************************

c     --- H.J.C. Berendsen; J.R. Grigera; T.P. Straatsma;
c     --- "The Missing Term in Effective Pair Potentials"
c     --- J. Phys. Chem. 91, 6269-6271 (1987).
      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'OW'
      vnlocal(1,2) = 'HW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ******************************************************************
c     * Nonbonded interactions for bond lengths                        *
c     ******************************************************************
c     --- Faux zeolite 4A has bond bending interactions, yet no bonds
c     --- need to use this bond type to make the logic come out OK

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_NONBOND)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 2
      vnlocal(1,1) = 'OZ'
      vnlocal(1,2) = 'Al'
      vnlocal(2,1) = 'OZ'
      vnlocal(2,2) = 'Si'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ******************************************************************
c     * Catlow Harmonic bond potentials                                *
c     ******************************************************************
c     --- C.R.A. Catlow, C.M. Freeman, B. Vessal, S.M. Tomlinson, 
c     --- M. Leslie J. Chem. Soc. Faraday Tran. 1999 87(13) 1947-1950.

c     --- Cc - Cc alkane carbon-carbon bond
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.4057d0
      vclocal(1) = 0.5d0*(60.0d0)*CON_EVTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'Cc'
      vnlocal(1,2) = 'Cc'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- Cc - Hc alkane carbon-hydrogen bond
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.054d0
      vclocal(1) = 0.5d0*(31.0d0)*CON_EVTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'Cc'
      vnlocal(1,2) = 'Hc'
      call twh_set_vibnames(type,nvslocal,vnlocal)


      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Nonbonded interaction between 1-3 beads                        *
c     ******************************************************************
c     --- only the angle type matters as this uses the nonbond 
c     --- interaction between whichever two beads make up the ends of 
c     --- the angle
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONBOND)
      bclocal(0) = 90.0d0
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'Al'
      bnlocal(1,2) = 'OZ'
      bnlocal(1,3) = 'Al'
      bnlocal(2,1) = 'Al'
      bnlocal(2,2) = 'OZ'
      bnlocal(2,3) = 'Si'
      bnlocal(3,1) = 'Si'
      bnlocal(3,2) = 'OZ'
      bnlocal(3,3) = 'Si'
      call twh_set_bendnames(type,nbslocal,bnlocal)
      
c     ******************************************************************
c     * Nonbonded interaction between 1-3 beads                        *
c     *       plus a simple harmonic term                              *
c     ******************************************************************
c     --- D.A. Faux, W. Smith, T.R. Forester J. Phys. Chem. B (1997)
c     --- 101, 1762-1768.
c     --- personal communication 02-12-2001 M.G. Martin with D.A. Faux
c     --- these are listed for angles Si-OZ-Al and Al-Si-OZ, but that 
c     --- is a typo.  They are really for the OZ-Si-OZ and OZ-Al-OZ
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NB_HARM)
      bclocal(0) = 109.47d0
      bclocal(1) = (4.5815)*(0.5d0)*(CON_EVTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'OZ'
      bnlocal(1,2) = 'Al'
      bnlocal(1,3) = 'OZ'
      bnlocal(2,1) = 'OZ'
      bnlocal(2,2) = 'Si'
      bnlocal(2,3) = 'OZ'
      call twh_set_bendnames(type,nbslocal,bnlocal)


c     ******************************************************************
c     * Catlow hydrocarbon angles - simple harmonic                    *
c     ******************************************************************
c     --- N. Raj, G. Sastre, C.R.A. Catlow J. Phys. Chem. B
c     --- 103, 11007-11015 (1999).
c     --- they got these parameters from kiselev, but it can be pretty 
c     --- hard to figure out exactly what kiselev did so take a look 
c     --- at both papers and decide for yourself who to cite
c     --- A.V. Keselev, A.A. Lopatkin, A.A. Shulga Zeolites 5,261 (1985)

c     --- Cc-Cc-Cc hydrocarbon
c     --- Hc-Cc-Cc hydrocarbon
c     --- Hc-Cc-Hc hydrocarbon
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (3.39d0)*(0.5d0)*CON_EVTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'Cc'
      bnlocal(1,2) = 'Cc'
      bnlocal(1,3) = 'Cc'
      bnlocal(2,1) = 'Hc'
      bnlocal(2,2) = 'Cc'
      bnlocal(2,3) = 'Cc'
      bnlocal(3,1) = 'Hc'
      bnlocal(3,2) = 'Cc'
      bnlocal(3,3) = 'Hc'
      call twh_set_bendnames(type,nbslocal,bnlocal)


c     ******************************************************************
c     * Rigid angle for SPC-E water at the tetrahedral position        *
c     ******************************************************************
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 109.47d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HW'
      bnlocal(1,2) = 'OW'
      bnlocal(1,3) = 'HW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * nonbond vdw with no coulombic                                  *
c     * any scaling factors for the 1-4 are accounted for in ffnonbond *
c     ******************************************************************
      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.0d0
c     --- N. Raj, G. Sastre, C.R.A. Catlow J. Phys. Chem. B
c     --- 103, 11007-11015 (1999).
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- * - Cc - Cc - *
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'Cc'
      tnlocal(1,3) = 'Cc'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     ******************************************************************
c     * nonbond vdw with coulombic                                     *
c     * any scaling factors for the 1-4 vdw are accounted for in       *
c     * ffnonbond, the coulombic is not scaled                         *
c     ******************************************************************
      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- * - OZ - Al - *
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'OZ'
      tnlocal(1,3) = 'Al'
      tnlocal(1,4) = 'wild'
c     --- * - OZ - Si - *
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'OZ'
      tnlocal(2,3) = 'Si'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)


      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
