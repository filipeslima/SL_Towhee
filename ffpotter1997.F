#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2014 Marcus G. Martin                       *
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
      subroutine twh_ffpotter1997(lfailure,ntypes,nbondtypes
     &     ,nangletypes,nbitypes)
c     ******************************************************************
c     * sets up parameters for potter19997                             *
c     *                                                                *
c     * originally written 04-26-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,nbitypes
c     --- local variables
      character*(FFNAMELEN) ffname

c     --- reference
c     S.C. Potter; D.J. Tildesley; A.N. Burgess; S.C. Rogers;
c     "A transferable potential model for the liquid-vapour equilibria
c     of fluoromethanes" Mol. Phys. 92, 825-833 (1997)

      ffname = 'Potter1997'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB plus manual')

      call twh_ffpotter1997_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffpotter1997_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffpotter1997_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffpotter1997_bond_increment(lfailure,ffname,nbitypes)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffpotter1997_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * sets up the nonbonded parameters for potter1997                *
c     *                                                                *
c     * originally written 04-26-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ntypes
c     --- local variables
      integer type,ii
      double precision masslocal
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      character*(FFNAMELEN) nbnlocal,cbnlocal,canlocal,ctnlocal

      type = 0

c     --- C
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.15d0
      nbclocal(2) = 54.6d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F(2)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F(2)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.975d0
      nbclocal(2) = 40.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F(3)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F(3)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.975d0
      nbclocal(2) = 40.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F(4)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F(4)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.975d0
      nbclocal(2) = 40.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H(1)
      type = type + 1
      masslocal = 1.00794d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H(1)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.17d0
      nbclocal(2) = 10.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H(2)
      type = type + 1
      masslocal = 1.00794d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H(2)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.17d0
      nbclocal(2) = 10.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffpotter1997_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * bond lengths for the potter1997 forcefield                     *
c     * all are rigid                                                  *
c     *                                                                *
c     * originally written 04-26-2006 by M.G. Martin                   *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nbondtypes
c     --- local variables
      integer type,nsame
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)

c     --- initialize types
      type = 0

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.360d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C-F(2)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C'
      vnlocal(nsame,2) = 'F(2)'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.340d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C-F(3)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C'
      vnlocal(nsame,2) = 'F(3)'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.335d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C-F(4)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C'
      vnlocal(nsame,2) = 'F(4)'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.090d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C-H(2)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C'
      vnlocal(nsame,2) = 'H(2)'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.080d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C-H(1)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C'
      vnlocal(nsame,2) = 'H(1)'
      call twh_set_vibnames(type,nsame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffpotter1997_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * sets up the angle terms for Potter1997                         *
c     *                                                                *
c     * originally written 04-26-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nangletypes
c     --- local variables
      integer type,nsame
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      character*(FFNAMELEN) bnlocal
      dimension bnlocal(MAXSAMEB,MAXBNAME)

c     --- initialize type count
      type = 0

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 113.61d0
      bclocal(1) = 232.3d0*(0.5d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H-C-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H(2)'
      bnlocal(nsame,2) = 'C'
      bnlocal(nsame,3) = 'H(2)'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.63d0
      bclocal(1) = 478.4d0*(0.5d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F-C-F
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F(4)'
      bnlocal(nsame,2) = 'C'
      bnlocal(nsame,3) = 'F(4)'
c     --- F-C-F
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F(3)'
      bnlocal(nsame,2) = 'C'
      bnlocal(nsame,3) = 'F(3)'
c     --- F-C-F
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F(2)'
      bnlocal(nsame,2) = 'C'
      bnlocal(nsame,3) = 'F(2)'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.64d0
      bclocal(1) = 323.1d0*(0.5d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H-C-F
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H(2)'
      bnlocal(nsame,2) = 'C'
      bnlocal(nsame,3) = 'F(2)'
c     --- H-C-F
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H(1)'
      bnlocal(nsame,2) = 'C'
      bnlocal(nsame,3) = 'F(3)'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffpotter1997_bond_increment(lfailure,ffname
     &     ,nbitypes)
c     ******************************************************************
c     * bond increments created based upon the charge assignments      *
c     *                                                                *
c     * originally written 04-26-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nbitypes
c     --- local variables
      integer type
      double precision qbivlocal
      character*(FFNAMELEN) qbinlocal
      dimension qbinlocal(2)

      type = 0

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.204d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- F(2)-C
      qbinlocal(1) = 'F(2)'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.182d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- F(3)-C
      qbinlocal(1) = 'F(3)'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -0.15d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- F(4)-C
      qbinlocal(1) = 'F(4)'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.097d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- H(2)-C
      qbinlocal(1) = 'H(2)'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.106d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- H(1)-C
      qbinlocal(1) = 'H(1)'
      qbinlocal(2) = 'C'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
