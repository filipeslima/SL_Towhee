#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_ffkfvbvs(lfailure,ntypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the KFvBvS       *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer j,ij,type
      character*(FFNAMELEN) ffname
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'KFvBvS'
      call twh_classical_potential(GLB_SET,'Exponential-6')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     --- exponential-6 potential 
c     --- NOTE: must explicitly set up all cross interactions
c     --- the proper format for this is to set up all interactions of
c     --- the lower index atom with all of the higher index atoms.
c     --- Explicit atom carbon Williams Exp-6 potential
c     --- U(r) = A*r^(-6) + B*exp[C*r]
      
c     --- The following parameters are for zeolites/silica and are
c     --- taken from
c     --- G.J. Kramer; N.P. Farragher; B.W.H. van Beest; 
c     --- R.A. van Santen; 
c     --- `Interatomic force fields for silicas, aluminophosphates,
c     --- and zeolites: Derivation based on ab initio calculations'
c     --- Phys. Rev. B 43, 5068-5080 (1991).
c     --- inevitably, this paper uses a different meaning for 
c     --- the A, B, and C parameters.  below is the conversion 
c     --- from their notation to the one used in this code
c     --- nbclocal(1) = -C
c     --- nbclocal(3) = A
c     --- nbclocal(4) = -b
      
c     --- (O) oxygen
      type = 1
      masslocal = 15.999d0
      elelocal = ' O'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
c     --- cross term of O with O (1-1)
      j = 1
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(175.0000d0)*CON_EVTOK
      nbclocal(2) = (1388.7730d0)*CON_EVTOK
      nbclocal(3) = -(2.76000d0)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Si (1-2)
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(133.5381d0)*CON_EVTOK
      nbclocal(2) = (18003.7572d0)*CON_EVTOK
      nbclocal(3) = -(4.87318d0)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Al (1-3)
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(130.5659d0)*CON_EVTOK
      nbclocal(2) = (16008.5345d0)*CON_EVTOK
      nbclocal(3) = -(4.79667d0)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with P (1-4)
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(19.8793d0)*CON_EVTOK
      nbclocal(2) = (9034.2080d0)*CON_EVTOK
      nbclocal(3) = -(5.19098d0)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Na (1-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(0.0000d0)*CON_EVTOK
      nbclocal(2) = (3542.2072d0)*CON_EVTOK
      nbclocal(3) = -(4.13455d0)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Si) silicon
c     --- note that in this force field silicon only has a nonzero
c     --- term with oxygen.  Otherwise the interactions are soley
c     --- coulombic.  Note that this is normally not a problem, but
c     --- will lead to disaster if used with their Cl parameters
      type = 2
      masslocal = 28.086d0
      elelocal = 'Si'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Si'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (Al) aluminum
c     --- note that in this force field aluminum only has a nonzero
c     --- term with oxygen.  Otherwise the interactions are soley
c     --- coulombic.  Note that this is normally not a problem, but
c     --- will lead to disaster if used with their Cl parameters
      type = 3
      masslocal = 26.982d0
      elelocal = 'Al'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Al'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (P) phosphorous
c     --- note that in this nbff phosphorous only has a nonzero
c     --- term with oxygen.  Otherwise the interactions are soley
c     --- coulombic.  Note that this is normally not a problem, but
c     --- will lead to disaster if used with their Cl parameters
      type = 4
      masslocal = 30.974d0
      elelocal = ' P'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     --- (Na) sodium
      type = 5
      masslocal = 22.990d0
      elelocal = 'Na'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
c     --- cross term of Na with Cl (5-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = -(521.3348d0)*CON_EVTOK
      nbclocal(2) = (5783.1124d0)*CON_EVTOK
      nbclocal(3) = -(3.11593d0)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Cl) chlorine
c     --- note that chlorine only has nonbonded parameters with 
c     --- sodium.  This will result in disaster if you use any of 
c     --- the other atoms in the force field as the coulombic terms
c     --- would go to infinity without any barrier
      type = 6
      masslocal = 35.453d0
      elelocal = 'Cl'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
