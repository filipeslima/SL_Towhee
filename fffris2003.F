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
      subroutine twh_fffris2003(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the fris2003     *
c     * forcefield                                                     *
c     * A.L. Frischnect; J.G. Curro; "Improved United Atom Force Field *
c     * for Poly(dimethylsiloxane), Macromolecules, 2003 36, 2122-2129 *
c     * the parameters listed here are those noted as hybrid/UA in that*
c     * paper                                                          *
c     *                                                                *
c     * originally written 10-13-2003 by M.G. Martin                   *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
     &        ,nimprotypes
c     --- local variables
      integer j,ij,type
      double precision eps96,sig96,eps126,sig126
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

      ffname = 'Fris2003'
c     --- 12-9-6 potential
      call twh_classical_potential(GLB_SET,'12-9-6')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     --- initialize types
      type = 0

c     --- They present the nonbond forcefield as a sum of a 12-6 
c     --- and a 9-6 potential.  These terms can be combined into 
c     --- the 12-9-6 potential.


c     --- (Si) silicon
      type = 1
      masslocal = 28.982d0
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
c     --- cross term of Si with Si (1-1)
      j = 1
      ij = twh_getnbtype(type,j)
      eps96 = (0.1310d0)*(CON_KCALTOK)
      sig96 = 4.29d0
      eps126 = (0.0d0)*(CON_KCALTOK)
      sig126 = 0.0d0
      nbclocal(1) = 4.0d0*(eps126)*(sig126**12)
      nbclocal(2) = 2.0d0*(eps96)*(sig96**9)
      nbclocal(3) = -3.0d0*(eps96)*(sig96**6) 
     &     - 4.0d0*eps126*(sig126**6)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with O (1-2)
      j = 2
      ij = twh_getnbtype(type,j)
      eps96 = (0.0772d0)*(CON_KCALTOK)
      sig96 = 3.94d0
      eps126 = (0.0d0)*(CON_KCALTOK)
      sig126 = 0.0d0
      nbclocal(1) = 4.0d0*(eps126)*(sig126**12)
      nbclocal(2) = 2.0d0*(eps96)*(sig96**9)
      nbclocal(3) = -3.0d0*(eps96)*(sig96**6) 
     &     - 4.0d0*eps126*(sig126**6)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with CH3 (1-3)
      j = 3
      ij = twh_getnbtype(type,j)
      eps96 = (0.0d0)*(CON_KCALTOK)
      sig96 = 0.0d0
      eps126 = (0.1596d0)*(CON_KCALTOK)
      sig126 = 3.83d0
      nbclocal(1) = 4.0d0*(eps126)*(sig126**12)
      nbclocal(2) = 2.0d0*(eps96)*(sig96**9)
      nbclocal(3) = -3.0d0*(eps96)*(sig96**6) 
     &     - 4.0d0*eps126*(sig126**6)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (O) oxygen
      type = 2
      masslocal = 15.999d0
      elelocal = 'O'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of O with O (2-2)
      j = 2
      ij = twh_getnbtype(type,j)
      eps96 = (0.0800d0)*(CON_KCALTOK)
      sig96 = 3.30d0
      eps126 = (0.0d0)*(CON_KCALTOK)
      sig126 = 0.0d0
      nbclocal(1) = 4.0d0*(eps126)*(sig126**12)
      nbclocal(2) = 2.0d0*(eps96)*(sig96**9)
      nbclocal(3) = -3.0d0*(eps96)*(sig96**6) 
     &     - 4.0d0*eps126*(sig126**6)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with CH3 (2-3)
      j = 3
      ij = twh_getnbtype(type,j)
      eps96 = (0.0d0)*(CON_KCALTOK)
      sig96 = 0.0d0
      eps126 = (0.1247d0)*(CON_KCALTOK)
      sig126 = 3.38d0
      nbclocal(1) = 4.0d0*(eps126)*(sig126**12)
      nbclocal(2) = 2.0d0*(eps96)*(sig96**9)
      nbclocal(3) = -3.0d0*(eps96)*(sig96**6) 
     &     - 4.0d0*eps126*(sig126**6)
      call twh_set_nbcoeff(ij,nbclocal) 

c     --- (CH3) methyl united-atom group
      type = 3
      masslocal = 15.0347d0
      elelocal = 'C'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of CH3 with CH3 (3-3)
      j = 3
      ij = twh_getnbtype(type,j)
      eps96 = (0.0d0)*(CON_KCALTOK)
      sig96 = 0.0d0
      eps126 = (0.1944d0)*(CON_KCALTOK)
      sig126 = 3.73d0
      nbclocal(1) = 4.0d0*(eps126)*(sig126**12)
      nbclocal(2) = 2.0d0*(eps96)*(sig96**9)
      nbclocal(3) = -3.0d0*(eps96)*(sig96**6) 
     &     - 4.0d0*eps126*(sig126**6)
      call twh_set_nbcoeff(ij,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     --- bond potentials
c     --- simple harmonic

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.64d0
      vclocal(1) = (350.12d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Si-O 
      vnlocal(1,1) = 'Si'
      vnlocal(1,2) = 'O'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.90d0
      vclocal(1) = (189.65d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Si-CH3
      vnlocal(1,1) = 'Si'
      vnlocal(1,2) = 'CH3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return
       
c     --- angle potentials
c     --- standard harmonic
c     --- Note: there is a typo in their paper.  The units for 
c     --- the angle potential are in fact kcal/mol rad^2

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 146.46d0
      bclocal(1) = (14.14d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- Si-O-Si
      nbslocal = 1
      bnlocal(1,1) = 'Si'
      bnlocal(1,2) = 'O'
      bnlocal(1,3) = 'Si'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.82d0
      bclocal(1) = (94.5d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- O-Si-O
      nbslocal = 1
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'Si'
      bnlocal(1,3) = 'O'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.24d0
      bclocal(1) = (49.97d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- CH3-Si-CH3
      nbslocal = 1
      bnlocal(1,1) = 'CH3'
      bnlocal(1,2) = 'Si'
      bnlocal(1,3) = 'CH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.69d0
      bclocal(1) = (49.97d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
c     --- O-Si-CH3
      nbslocal = 1
      bnlocal(1,1) = 'O'
      bnlocal(1,2) = 'Si'
      bnlocal(1,3) = 'CH3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     --- torsion potentials

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.225)*(CON_KCALTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- Si-O-Si-O
      tnlocal(1,1) = 'Si'
      tnlocal(1,2) = 'O'
      tnlocal(1,3) = 'Si'
      tnlocal(1,4) = 'O'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.01)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- Si-O-Si-CH3
      tnlocal(1,1) = 'Si'
      tnlocal(1,2) = 'O'
      tnlocal(1,3) = 'Si'
      tnlocal(1,4) = 'CH3'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     --- no impropers
      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

      return
      end
