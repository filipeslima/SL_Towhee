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
      subroutine twh_ffcui1998(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nbitypes)
c     ******************************************************************
c     * sets up parameters for Cui1998                                 *
c     *                                                                *
c     * originally written 05-02-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
     &        ,nbitypes
c     --- local variables
      character*(FFNAMELEN) ffname

c     --- reference
c     S.T. Cui; J.I. Siepmann; H.D. Cochran; P.T. Cummings;
c     "Intermolecular potentials and vapor-liquid phase equilibria
c     of perfluorinated alkanes",
c     Fluid Phase Equilib. 146, 15-61 (1998)

      ffname = 'Cui1998'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

      call twh_ffcui1998_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffcui1998_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffcui1998_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffcui1998_dihedral(lfailure,ffname,ndihedtypes)
      if ( lfailure ) return
      call twh_ffcui1998_bond_increment(lfailure,ffname,nbitypes)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcui1998_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * sets up the nonbonded parameters for Cui1998                   *
c     *                                                                *
c     * originally written 05-02-2006 by M.G. Martin                   *
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

c     --- CF3(T)
      type = type + 1
      masslocal = 12.011d0 + 3.0d0*(18.998d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CF3(T)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.60d0
      nbclocal(2) = 79.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- CF2(T)
      type = type + 1
      masslocal = 12.011d0 + 2.0d0*(18.998d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CF2(T)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.60d0
      nbclocal(2) = 30.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- CF3(M)
      type = type + 1
      masslocal = 12.011d0 + 3.0d0*(18.998d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CF3(M)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.65d0
      nbclocal(2) = 100.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- CF2(M)
      type = type + 1
      masslocal = 12.011d0 + 2.0d0*(18.998d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CF2(M)'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.65d0
      nbclocal(2) = 30.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcui1998_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * bond lengths for the cui1998 forcefield                        *
c     * rigid C-C bond lengths                                         *
c     *                                                                *
c     * originally written 05-02-2006 by M.G. Martin                   *
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
      vclocal(0) = 1.54d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C-C
      nsame = nsame + 1
      vnlocal(nsame,1) = 'Cwild'
      vnlocal(nsame,2) = 'Cwild'
      call twh_set_vibnames(type,nsame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcui1998_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * sets up the angle terms for Cui1998                           *
c     *                                                                *
c     * originally written 05-02-2006 by M.G. Martin                   *
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
      bclocal(0) = 114.6d0
      bclocal(1) = 62500d0*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- CFx - CF2(T) - CFy
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Cwild'
      bnlocal(nsame,2) = 'CF2(T)'
      bnlocal(nsame,3) = 'Cwild'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 114.0d0
      bclocal(1) = 62500d0*(0.5d0)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- CFx - CF2(M) - CFy
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Cwild'
      bnlocal(nsame,2) = 'CF2(M)'
      bnlocal(nsame,3) = 'Cwild'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcui1998_dihedral(lfailure,ffname,ndihedtypes)
c     ******************************************************************
c     * dihedral terms for Cui1998                                     *
c     *                                                                *
c     * originally written 05-02-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ndihedtypes
c     --- local variables
      integer type,nsame
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)

      type = 0

      type = type + 1
      call twh_torstyle(GLB_SET,type,10)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,7)
      tclocal(0) = 595.4d0
      tclocal(1) = -282.7d0
      tclocal(2) = 1355.2d0
      tclocal(3) = 6800.0d0
      tclocal(4) = -7875.3d0
      tclocal(5) = -14168.0d0
      tclocal(6) = 9213.7d0
      tclocal(7) = 4123.7d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- CFx - CF2(T) - CF2(T) - CFy
      nsame = nsame + 1
      tnlocal(nsame,1) = 'Cwild'
      tnlocal(nsame,2) = 'CF2(T)'
      tnlocal(nsame,3) = 'CF2(T)'
      tnlocal(nsame,4) = 'Cwild'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,22)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(0) = 269.6d0
      tclocal(1) = 503.6d0
      tclocal(2) = 679.9d0
      tclocal(3) = 3.0d0
      tclocal(4) = 420.6d0
      tclocal(5) = 30.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- CFx - CF2(M) - CF2(M) - CFy
      nsame = nsame + 1
      tnlocal(nsame,1) = 'Cwild'
      tnlocal(nsame,2) = 'CF2(M)'
      tnlocal(nsame,3) = 'CF2(M)'
      tnlocal(nsame,4) = 'Cwild'
      call twh_set_tornames(type,nsame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffcui1998_bond_increment(lfailure,ffname,nbitypes)
c     ******************************************************************
c     * all of the charges are assigned in the qbasevalue in nonbond so*
c     * the only thing here is a global wildcard 0.0                   *
c     *                                                                *
c     * originally written 05-02-2006 by M.G. Martin                   *
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
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- wild
      qbinlocal(1) = 'wild'
      qbinlocal(2) = 'wild'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
