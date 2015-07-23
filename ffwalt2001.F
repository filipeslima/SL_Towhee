#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2014 Marcus G. Martin                       *
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
c     *                                                                *
c     * See the file towhee.F for more information about the code      *
c     ******************************************************************
      subroutine twh_ffwalt2001(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the walt2001     *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 02-02-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nbitypes
c     --- local variables
      integer type,jtype,ij,isame
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal
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
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)

      ffname = 'Walt2001'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     --- reference
c     J.H. Walther; R. Jaffe; T. Halicioglu; P. Koumoutsakos;
c     "Carbon Nanotubes in Water: Structural Characteristics and 
c     Energetics", J. Phys. Chem. B 105, (2001), 9980-9987.

c     Carbon in a nanotube
      type = 1
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- interaction with self (1-1)
      jtype = 1
      ij = twh_getnbtype(type,jtype)
      nbclocal(1) = 3.851d0
      nbclocal(2) = (0.4396d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- interaction of C with O (1-2)
      jtype = 2
      ij = twh_getnbtype(type,jtype)
      nbclocal(1) = 3.19d0
      nbclocal(2) = (0.3126d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- interaction of C with H (1-3)
      jtype = 3
      ij = twh_getnbtype(type,jtype)
      nbclocal(1) = 1.0d0
      nbclocal(2) = (0.0d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     Oxygen in water
      type = 2
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- interaction with self (2-2)
      jtype = 2
      ij = twh_getnbtype(type,jtype)
      nbclocal(1) = 3.166d0
      nbclocal(2) = (0.6502d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- interaction of O with H (2-3)
      jtype = 3
      ij = twh_getnbtype(type,jtype)
      nbclocal(1) = 1.0d0
      nbclocal(2) = (0.0d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ij,nbclocal) 

c     Hydrogen in water
      type = 3
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- interaction with self (3-3)
      jtype = 3
      ij = twh_getnbtype(type,jtype)
      nbclocal(1) = 1.0d0
      nbclocal(2) = (0.0d0)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ij,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond lengths - Morse                                           *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_MORSE)
      vclocal(0) = 1.418d0
      vclocal(1) = (478.9d0)*CON_KJTOK
      vclocal(2) = -(2.1867d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- carbon bonds in a nanotube
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     ******************************************************************
c     * bond lengths - Standard Harmonic                               *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0d0
      vclocal(1) = (0.5d0)*(4637.0d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- O-H bonds in water
      isame = isame + 1
      vnlocal(isame,1) = 'O'
      vnlocal(isame,2) = 'H'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * angles - harmonic cosine                                       *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
      bclocal(0) = 120.00d0
      bclocal(1) = (562.2d0)*(0.5d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C - C - C in a nanotube
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ******************************************************************
c     * angles - standard harmonic                                     *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (383.0d0)*(0.5d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H - O - H in flexible SPC water
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'O'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * torsions - 2-fold cosine                                       *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,16)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(1) = (0.5d0)*(25.12d0)*CON_KJTOK
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     C - C - C - C
      isame = isame + 1
      tnlocal(isame,1) = 'C'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'C'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'C'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.41d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'O'
      qbinlocal(2) = 'H'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
