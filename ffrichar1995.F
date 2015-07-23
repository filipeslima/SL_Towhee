#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2008-2014 Marcus G. Martin                       *
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
      subroutine twh_ffrichar1995(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the Richar1995   *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 08-11-2008 M.G. Martin                      *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nbitypes
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
      double precision charge
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
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)

      ffname = 'Richar1995'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     --- Reference:
c     A.J. Richards; K. Watanabe; N. Austin; M.R. Stapleton;
c     "Computer Simulation of the Gas Separation Properties of Zeolite
c     Li-X", Journal of Porous Materials, 2, 43-49 (1995).

c     --- NOTE: must explicitly set up all cross interactions
c     --- the proper format for this is to set up all interactions of
c     --- the lower index atom with all of the higher index atoms.
      
c     --- (O) the oxygen in O2
      type = 1
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      nbnlocal = 'O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of O with O (1-1)
      j = 1
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.09d0
      nbclocal(2) = (0.0884d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with N (1-2)
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.21d0
      nbclocal(2) = (0.0799d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Si (1-3)
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 1.88d0
      nbclocal(2) = (0.0702d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Al (1-4)
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.05d0
      nbclocal(2) = (0.0701d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Oz (1-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.90d0
      nbclocal(2) = (0.194d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with Li (1-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.08d0
      nbclocal(2) = (0.0291d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of O with mp (1-7)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
     
c     --- (N) Nitrogen in N2
      type = 2
      masslocal = 14.0067d0
      elelocal = ' N'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      nbnlocal = 'N'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of N with N (2-2)
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.32d0
      nbclocal(2) = (0.0723d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of N with Si (2-3)
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 1.99d0
      nbclocal(2) = (0.0563d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of N with Al (2-4)
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.16d0
      nbclocal(2) = (0.0578d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of N with Oz (2-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.01d0
      nbclocal(2) = (0.175d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of N with Li (2-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.19d0
      nbclocal(2) = (0.0240d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Si) Silicon in zeolite Li-X
      type = 3
      masslocal = 28.086d0
      elelocal = 'Si'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      charge = 304.0d0/192.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      nbnlocal = 'Si'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- Si has no interactions with atoms equal or greater index
c     --- cross term of Si with Si (3-3)
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with Al (3-4)
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with Oz (3-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with Li (3-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Si with mp (3-7)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Al) Aluminum in zeolite Li-X
      type = 4
      masslocal = 26.982d0
      elelocal = 'Al'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      charge = 304.0d0/192.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      nbnlocal = 'Al'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- Al has no interactions with atoms equal or greater index
c     --- cross term of Al with Al (4-4)
      j = 4
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Al with Oz (4-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Al with Li (4-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Al with mp (4-7)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (Oz) the oxygen in zeolite Li-X
      type = 5
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,-1.0d0)
      nbnlocal = 'Oz'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- Oz has no interactions with atoms equal or greater index
c     --- cross term of Oz with Oz (5-5)
      j = 5
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Oz with Li (5-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Oz with mp (5-7)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 

c     --- (Li) Lithium in zeolite Li-X
      type = 6
      masslocal = 6.941d0
      elelocal = 'Li'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,1.0d0)
      nbnlocal = 'Li'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- Li has no interactions with atoms equal or greater index
c     --- cross term of Li with Li (6-6)
      j = 6
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of Li with mp (6-7)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
      
c     --- (mp) midpoint in O2 or N2
      type = 7
      masslocal = 0.0d0
      elelocal = '-'
      bplocal = 'sp'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      nbnlocal = 'mp'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- mp has no interactions with atoms equal or greater index
c     --- cross term of mp with mp (7-7)
      j = 7
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*CON_KCALTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Fixed bond lengths                                             *
c     ******************************************************************

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.549d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Nitrogen with its midpoint
      vnlocal(1,1) = 'N'
      vnlocal(1,2) = 'mp'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.604d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
c     --- Oxygen with its midpoint
      vnlocal(1,1) = 'O'
      vnlocal(1,2) = 'mp'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ******************************************************************
c     * no interaction for bond lengths                                *
c     ******************************************************************
c     --- zeolite bonds that do not interact as the zeolite is rigid

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_NONE)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- everything else is no interaction
c     --- Si - Oz
      nvslocal = nvslocal + 1
      vnlocal(nvslocal,1) = 'Si'
      vnlocal(nvslocal,2) = 'Oz'
c     --- Al - Oz
      nvslocal = nvslocal + 1
      vnlocal(nvslocal,1) = 'Al'
      vnlocal(nvslocal,2) = 'Oz'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * rigid angles                                                   *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 180.0d0
      bclocal(1) = 1.0d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- O2 and N2 dimers are modeled as trimers with a midpoint
      bnlocal(1,1) = 'wild'
      bnlocal(1,2) = 'mp'
      bnlocal(1,3) = 'wild'
      call twh_set_bendnames(type,nbslocal,bnlocal)
      
c     ******************************************************************
c     * no interactions angles                                         *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONE)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     --- all other angles are the zeolite and have no interactions
      bnlocal(1,1) = 'wild'
      bnlocal(1,2) = 'wild'
      bnlocal(1,3) = 'wild'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * null torsions                                                  *
c     ******************************************************************
      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- all torsions are null in the zeolite
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'wild'
      tnlocal(1,3) = 'wild'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.405d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     N - mp in N2
      qbinlocal(1) = 'N'
      qbinlocal(2) = 'mp'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.112d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     O - mp in O2
      qbinlocal(1) = 'O'
      qbinlocal(2) = 'mp'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     everything else does not exchange charge
      qbinlocal(1) = 'wild'
      qbinlocal(2) = 'wild'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
