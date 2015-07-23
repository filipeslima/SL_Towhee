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
      subroutine twh_ffshah2004(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nbitypes)
c     ******************************************************************
c     * sets up parameters for shah2004                                *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nbitypes
c     --- local variables
      character*(FFNAMELEN) ffname

c     --- reference
c     J.K. Shah; E.J. Maginn; "A Monte Carlo simulation study of the
c     ionic liquid 1-n-butyl-3-methylimidazolium hexafluorophosphate:
c     liquid structure, volumetric properties and infinite dilution
c     solution thermodynamics of CO2", Fluid Phase Equilib. 222-223
c     (2004) 195-203.

      ffname = 'Shah2004'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

      call twh_ffshah2004_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffshah2004_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffshah2004_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffshah2004_dihedral(lfailure,ffname,ndihedtypes)
      if ( lfailure ) return
      call twh_ffshah2004_bond_increment(lfailure,ffname,nbitypes)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffshah2004_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * sets up the nonbonded parameters for shah2004                  *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
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

c     --- uses the full van der waals 1-4 (personal communication 
c     --- between E.J. Maginn and M.G. Martin 04-21-2006)

      type = 0

c     --- N1
      type = type + 1
      masslocal = 14.0070d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' N')
      call twh_qbasevalue(GLB_SET,type,0.111d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = 0.710d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- N3
      type = type + 1
      masslocal = 14.0070d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' N')
      call twh_qbasevalue(GLB_SET,type,0.133d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = 0.710d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C2
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 1.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.233d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.880d0
      nbclocal(2) = 0.443d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C4
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 1.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.040d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.880d0
      nbclocal(2) = 0.443d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C5
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 1.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.010d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.880d0
      nbclocal(2) = 0.443d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C6
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 3.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.183d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = 0.865d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C7
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 2.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.195d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = 0.493d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C8
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 2.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.066d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = 0.493d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C9
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 2.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.128d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = 0.493d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C10
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 3.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.043d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C10'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = 0.732d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- P
      type = type + 1
      masslocal = 30.9740d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' P')
      call twh_qbasevalue(GLB_SET,type,1.460d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.740d0
      nbclocal(2) = 0.836d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F (generic F term that uses an average bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F1 (specific F term that uses an individual bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F2 (specific F term that uses an individual bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F3 (specific F term that uses an individual bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F4 (specific F term that uses an individual bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F5 (specific F term that uses an individual bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F6 (specific F term that uses an individual bond length)
      type = type + 1
      masslocal = 18.9980d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,-0.394d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.118d0
      nbclocal(2) = 0.255d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffshah2004_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * bond lengths for the shah2004 forcefield                       *
c     * these are not explicitly specified and instead must be         *
c     * inferred from their example structure                          *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
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
      double precision dist,apos,bpos
      dimension apos(3),bpos(3)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)

c     --- initialize types
      type = 0

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -1.542d0
      apos(2) = -2.310d0
      apos(3) = -1.640d0
      bpos(1) = -2.041d0
      bpos(2) = -1.177d0
      bpos(3) = -0.852d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6 - N3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C6'
      vnlocal(nsame,2) = 'N3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -1.322d0
      apos(2) = -0.100d0
      apos(3) = -0.530d0
      bpos(1) = -2.041d0
      bpos(2) = -1.177d0
      bpos(3) = -0.852d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2 - N3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C2'
      vnlocal(nsame,2) = 'N3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -3.274d0
      apos(2) = -1.083d0
      apos(3) = -0.232d0
      bpos(1) = -2.041d0
      bpos(2) = -1.177d0
      bpos(3) = -0.852d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4 - N3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C4'
      vnlocal(nsame,2) = 'N3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -3.274d0
      apos(2) = -1.083d0
      apos(3) = -0.232d0
      bpos(1) = -3.279d0
      bpos(2) = 0.082d0
      bpos(3) = 0.469d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4 - C5
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C4'
      vnlocal(nsame,2) = 'C5'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -2.0478d0
      apos(2) = 0.682d0
      apos(3) = 0.272d0
      bpos(1) = -3.279d0
      bpos(2) = 0.082d0
      bpos(3) = 0.469d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1 - C5
      nsame = nsame + 1
      vnlocal(nsame,1) = 'N1'
      vnlocal(nsame,2) = 'C5'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -2.0478d0
      apos(2) = 0.682d0
      apos(3) = 0.272d0
      bpos(1) = -1.322d0
      bpos(2) = -0.100d0
      bpos(3) = -0.530d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1 - C2
      nsame = nsame + 1
      vnlocal(nsame,1) = 'N1'
      vnlocal(nsame,2) = 'C2'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -2.0478d0
      apos(2) = 0.682d0
      apos(3) = 0.272d0
      bpos(1) = -1.574d0
      bpos(2) = 1.940d0
      bpos(3) = 0.888d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1 - C7
      nsame = nsame + 1
      vnlocal(nsame,1) = 'N1'
      vnlocal(nsame,2) = 'C7'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -1.993d0
      apos(2) = 3.187d0
      apos(3) = 0.103d0
      bpos(1) = -1.574d0
      bpos(2) = 1.940d0
      bpos(3) = 0.888d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8 - C7
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C8'
      vnlocal(nsame,2) = 'C7'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = -1.993d0
      apos(2) = 3.187d0
      apos(3) = 0.103d0
      bpos(1) = -1.340d0
      bpos(2) = 3.348d0
      bpos(3) = -1.278d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8 - C9
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C8'
      vnlocal(nsame,2) = 'C9'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 0.173d0
      apos(2) = 3.582d0
      apos(3) = -1.238d0
      bpos(1) = -1.340d0
      bpos(2) = 3.348d0
      bpos(3) = -1.278d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C10 - C9
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C10'
      vnlocal(nsame,2) = 'C9'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 1.995d0
      apos(2) = -1.001d0
      apos(3) = 0.512d0
      bpos(1) = 1.428d0
      bpos(2) = -1.098d0
      bpos(3) = -1.066d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F1
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F1'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 1.995d0
      apos(2) = -1.001d0
      apos(3) = 0.512d0
      bpos(1) = 0.502d0
      bpos(2) = -1.550d0
      bpos(3) = 1.005d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F2
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F2'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 1.995d0
      apos(2) = -1.001d0
      apos(3) = 0.512d0
      bpos(1) = 2.529d0
      bpos(2) = -2.527d0
      bpos(3) = 0.414d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 1.995d0
      apos(2) = -1.001d0
      apos(3) = 0.512d0
      bpos(1) = 2.462d0
      bpos(2) = -0.861d0
      bpos(3) = 2.055d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F4
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F4'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 1.995d0
      apos(2) = -1.001d0
      apos(3) = 0.512d0
      bpos(1) = 3.400d0
      bpos(2) = -0.413d0
      bpos(3) = -0.031d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F5
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F5'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c     --- compute the bond length from the coordinates
      apos(1) = 1.995d0
      apos(2) = -1.001d0
      apos(3) = 0.512d0
      bpos(1) = 1.366d0
      bpos(2) = 0.553d0
      bpos(3) = 0.557d0
      call twh_simple_distance(dist,apos,bpos)
      vclocal(0) = dist
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F6
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F6'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
c     --- generic bond length from Morrow2002, but rigid
c     --- suggested for use by E.J. Maginn (pers. comm with M.G. Martin
c     --- 04-19-2006)
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.646d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F (generic)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F'
      call twh_set_vibnames(type,nsame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffshah2004_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * sets up the angle terms for Shah2004                           *
c     *                                                                *
c     * originally written 04-20-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from this subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nangletypes
c     --- local variables
      integer type,nsame,index
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      character*(FFNAMELEN) bnlocal
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      double precision xcord,ycord,zcord
      dimension xcord(3),ycord(3),zcord(3)
      double precision angle,radtodeg

c     --- retrieve constant
      call twh_constant_radtodeg(GLB_GET,radtodeg)

c     --- initialize type count
      type = 0

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C6
      index = 1
      xcord(index) = -1.542d0
      ycord(index) = -2.310d0
      zcord(index) = -1.640d0
c     N3
      index = 2
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 3
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6-N3-C4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C6'
      bnlocal(nsame,2) = 'N3'
      bnlocal(nsame,3) = 'C4'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C6
      index = 1
      xcord(index) = -1.542d0
      ycord(index) = -2.310d0
      zcord(index) = -1.640d0
c     N3
      index = 2
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C2
      index = 3
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6-N3-C2
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C6'
      bnlocal(nsame,2) = 'N3'
      bnlocal(nsame,3) = 'C2'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     N3
      index = 1
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C2
      index = 2
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N1
      index = 3
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N3-C2-N1
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N3'
      bnlocal(nsame,2) = 'C2'
      bnlocal(nsame,3) = 'N1'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C2
      index = 1
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N3
      index = 2
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 3
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N3-C4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C2'
      bnlocal(nsame,2) = 'N3'
      bnlocal(nsame,3) = 'C4'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     N3
      index = 1
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 2
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 3
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N3-C4-C5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N3'
      bnlocal(nsame,2) = 'C4'
      bnlocal(nsame,3) = 'C5'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C4
      index = 1
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 2
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
c     N1
      index = 3
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4-C5-N1
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C4'
      bnlocal(nsame,2) = 'C5'
      bnlocal(nsame,3) = 'N1'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C2
      index = 1
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N1
      index = 2
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C5
      index = 3
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N1-C5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C2'
      bnlocal(nsame,2) = 'N1'
      bnlocal(nsame,3) = 'C5'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C2
      index = 1
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N1
      index = 2
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C7
      index = 3
      xcord(index) = -1.574d0
      ycord(index) = 1.940d0
      zcord(index) = 0.888d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N1-C7
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C2'
      bnlocal(nsame,2) = 'N1'
      bnlocal(nsame,3) = 'C7'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C5
      index = 1
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
c     N1
      index = 2
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C7
      index = 3
      xcord(index) = -1.574d0
      ycord(index) = 1.940d0
      zcord(index) = 0.888d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C5-N1-C7
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C5'
      bnlocal(nsame,2) = 'N1'
      bnlocal(nsame,3) = 'C7'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     N1
      index = 1
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C7
      index = 2
      xcord(index) = -1.574d0
      ycord(index) = 1.940d0
      zcord(index) = 0.888d0
c     C8
      index = 3
      xcord(index) = -1.993d0
      ycord(index) = 3.187d0
      zcord(index) = 0.103d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1-C7-C8
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N1'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'C8'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C7
      index = 1
      xcord(index) = -1.574d0
      ycord(index) = 1.940d0
      zcord(index) = 0.888d0
c     C8
      index = 2
      xcord(index) = -1.993d0
      ycord(index) = 3.187d0
      zcord(index) = 0.103d0
c     C9
      index = 3
      xcord(index) = -1.340d0
      ycord(index) = 3.348d0
      zcord(index) = -1.278d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7-C8-C9
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C7'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'C9'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     C8
      index = 1
      xcord(index) = -1.993d0
      ycord(index) = 3.187d0
      zcord(index) = 0.103d0
c     C9
      index = 2
      xcord(index) = -1.340d0
      ycord(index) = 3.348d0
      zcord(index) = -1.278d0
c     C10
      index = 3
      xcord(index) = 0.173d0
      ycord(index) = 3.582d0
      zcord(index) = -1.238d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8-C9-C10
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C8'
      bnlocal(nsame,2) = 'C9'
      bnlocal(nsame,3) = 'C10'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F1
      index = 1
      xcord(index) = 1.428d0
      ycord(index) = -1.098d0
      zcord(index) = -1.066d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F2
      index = 3
      xcord(index) = 0.502d0
      ycord(index) = -1.550d0
      zcord(index) = 1.005d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F1-P-F2
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F1'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F2'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F1
      index = 1
      xcord(index) = 1.428d0
      ycord(index) = -1.098d0
      zcord(index) = -1.066d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F3
      index = 3
      xcord(index) = 2.529d0
      ycord(index) = -2.527d0
      zcord(index) = 0.414d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F1-P-F3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F1'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F3'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F1
      index = 1
      xcord(index) = 1.428d0
      ycord(index) = -1.098d0
      zcord(index) = -1.066d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F4
      index = 3
      xcord(index) = 2.462d0
      ycord(index) = -0.861d0
      zcord(index) = 2.055d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F1-P-F4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F1'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F4'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F1
      index = 1
      xcord(index) = 1.428d0
      ycord(index) = -1.098d0
      zcord(index) = -1.066d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F5
      index = 3
      xcord(index) = 3.400d0
      ycord(index) = -0.413d0
      zcord(index) = -0.031d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F1-P-F5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F1'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F5'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F1
      index = 1
      xcord(index) = 1.428d0
      ycord(index) = -1.098d0
      zcord(index) = -1.066d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F6
      index = 3
      xcord(index) = 1.366d0
      ycord(index) = 0.553d0
      zcord(index) = 0.557d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F1-P-F6
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F1'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F6'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F2
      index = 1
      xcord(index) = 0.502d0
      ycord(index) = -1.550d0
      zcord(index) = 1.005d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F3
      index = 3
      xcord(index) = 2.529d0
      ycord(index) = -2.527d0
      zcord(index) = 0.414d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F2-P-F3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F2'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F3'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F2
      index = 1
      xcord(index) = 0.502d0
      ycord(index) = -1.550d0
      zcord(index) = 1.005d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F4
      index = 3
      xcord(index) = 2.462d0
      ycord(index) = -0.861d0
      zcord(index) = 2.055d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F2-P-F4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F2'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F4'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F2
      index = 1
      xcord(index) = 0.502d0
      ycord(index) = -1.550d0
      zcord(index) = 1.005d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F5
      index = 3
      xcord(index) = 3.400d0
      ycord(index) = -0.413d0
      zcord(index) = -0.031d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F2-P-F5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F2'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F5'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F2
      index = 1
      xcord(index) = 0.502d0
      ycord(index) = -1.550d0
      zcord(index) = 1.005d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F6
      index = 3
      xcord(index) = 1.366d0
      ycord(index) = 0.553d0
      zcord(index) = 0.557d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F2-P-F6
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F2'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F6'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F3
      index = 1
      xcord(index) = 2.529d0
      ycord(index) = -2.527d0
      zcord(index) = 0.414d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F4
      index = 3
      xcord(index) = 2.462d0
      ycord(index) = -0.861d0
      zcord(index) = 2.055d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F3-P-F4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F3'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F4'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F3
      index = 1
      xcord(index) = 2.529d0
      ycord(index) = -2.527d0
      zcord(index) = 0.414d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F5
      index = 3
      xcord(index) = 3.400d0
      ycord(index) = -0.413d0
      zcord(index) = -0.031d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F3-P-F5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F3'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F5'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F3
      index = 1
      xcord(index) = 2.529d0
      ycord(index) = -2.527d0
      zcord(index) = 0.414d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F6
      index = 3
      xcord(index) = 1.366d0
      ycord(index) = 0.553d0
      zcord(index) = 0.557d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F3-P-F6
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F3'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F6'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F4
      index = 1
      xcord(index) = 2.462d0
      ycord(index) = -0.861d0
      zcord(index) = 2.055d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F5
      index = 3
      xcord(index) = 3.400d0
      ycord(index) = -0.413d0
      zcord(index) = -0.031d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F4-P-F5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F4'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F5'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F4
      index = 1
      xcord(index) = 2.462d0
      ycord(index) = -0.861d0
      zcord(index) = 2.055d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F6
      index = 3
      xcord(index) = 1.366d0
      ycord(index) = 0.553d0
      zcord(index) = 0.557d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F4-P-F6
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F4'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F6'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
c     F5
      index = 1
      xcord(index) = 3.400d0
      ycord(index) = -0.413d0
      zcord(index) = -0.031d0
c     P
      index = 2
      xcord(index) = 1.995d0
      ycord(index) = -1.001d0
      zcord(index) = 0.512d0
c     F6
      index = 3
      xcord(index) = 1.366d0
      ycord(index) = 0.553d0
      zcord(index) = 0.557d0
      call twh_simple_angle(angle,xcord,ycord,zcord)
      bclocal(0) = angle*radtodeg
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F5-P-F6
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F5'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F6'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_MULTI_RIGID)
      bclocal(0) = 90.0d0
      bclocal(1) = 180.0d0
      bclocal(2) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F-P-F for symmetric PF6
      nsame = nsame + 1
      bnlocal(nsame,1) = 'F'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'F'
c     --- set nbsame
      call twh_set_bendnames(type,nsame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffshah2004_dihedral(lfailure,ffname,ndihedtypes)
c     ******************************************************************
c     * dihedral terms for Shah2004                                    *
c     *                                                                *
c     * originally written 04-24-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ndihedtypes
c     --- local variables
      integer type,nsame,index
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision xcord,ycord,zcord
      dimension xcord(4),ycord(4),zcord(4)
      double precision dihedral

c     --- uses 1-4 vdw plus 1-4 coulomb scaled by 1/2 (personal 
c     --- communication between E.J. Maginn and M.G. Martin 04-21-2006)

      type = 0

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (5.89d0)*CON_KJTOK*(0.5d0)
      tclocal(2) =  (-1.13d0)*CON_KJTOK*(0.5d0)
      tclocal(3) =  (13.17d0)*CON_KJTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7-C8-C9-C10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C7'
      tnlocal(nsame,2) = 'C8'
      tnlocal(nsame,3) = 'C9'
      tnlocal(nsame,4) = 'C10'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (11.41d0)*CON_KJTOK*(0.5d0)
      tclocal(2) =  (0.96d0)*CON_KJTOK*(0.5d0)
      tclocal(3) =  (2.05d0)*CON_KJTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1-C7-C8-C9
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'C9'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (-5.85d0)*CON_KJTOK*(0.5d0)
      tclocal(2) =  (-1.80d0)*CON_KJTOK*(0.5d0)
      tclocal(3) =  (0.00d0)*CON_KJTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N1-C7-C8
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'C8'
c     --- C5-N1-C7-C8
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C5'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'C8'
      call twh_set_tornames(type,nsame,tnlocal)

c     --- now need to assign the torsions for the ring - this was held
c     --- rigid in their work.

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     N3
      index = 1
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 2
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 3
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
c     N1
      index = 4
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- N3-C4-C5-N1
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N3'
      tnlocal(nsame,2) = 'C4'
      tnlocal(nsame,3) = 'C5'
      tnlocal(nsame,4) = 'N1'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     C4
      index = 1
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 2
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
c     N1
      index = 3
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C2
      index = 4
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4-C5-N1-C2
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C4'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C2'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     C4
      index = 1
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 2
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
c     N1
      index = 3
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C7
      index = 4
      xcord(index) = -1.574d0
      ycord(index) = 1.940d0
      zcord(index) = 0.888d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4-C5-N1-C7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C4'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C7'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     C5
      index = 1
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
c     N1
      index = 2
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C2
      index = 3
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N3
      index = 4
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C5-N1-C2-N3
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C5'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C2'
      tnlocal(nsame,4) = 'N3'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     C7
      index = 1
      xcord(index) = -1.574d0
      ycord(index) = 1.940d0
      zcord(index) = 0.888d0
c     N1
      index = 2
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C2
      index = 3
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N3
      index = 4
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7-N1-C2-N3
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C7'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C2'
      tnlocal(nsame,4) = 'N3'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     N1
      index = 1
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C2
      index = 2
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N3
      index = 3
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 4
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1-C2-N3-C4
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C4'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     N1
      index = 1
      xcord(index) = -2.0478d0
      ycord(index) = 0.682d0
      zcord(index) = 0.272d0
c     C2
      index = 2
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N3
      index = 3
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C6
      index = 4
      xcord(index) = -1.542d0
      ycord(index) = -2.310d0
      zcord(index) = -1.640d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1-C2-N3-C6
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C6'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     C2
      index = 1
      xcord(index) = -1.322d0
      ycord(index) = -0.100d0
      zcord(index) = -0.530d0
c     N3
      index = 2
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 3
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 4
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N3-C4-C5
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N3'
      tnlocal(nsame,3) = 'C4'
      tnlocal(nsame,4) = 'C5'
      call twh_set_tornames(type,nsame,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,21)
      call twh_loftor(GLB_SET,type,.false.)
      call twh_ntorloop(GLB_SET,type,1)
c     C6
      index = 1
      xcord(index) = -1.542d0
      ycord(index) = -2.310d0
      zcord(index) = -1.640d0
c     N3
      index = 2
      xcord(index) = -2.041d0
      ycord(index) = -1.177d0
      zcord(index) = -0.852d0
c     C4
      index = 3
      xcord(index) = -3.274d0
      ycord(index) = -1.083d0
      zcord(index) = -0.232d0
c     C5
      index = 4
      xcord(index) = -3.279d0
      ycord(index) = 0.082d0
      zcord(index) = 0.469d0
      call twh_simple_dihedral(dihedral,xcord,ycord,zcord)
      tclocal(0) = 1.0d-05
      tclocal(1) = dihedral
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6-N3-C4-C5
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C6'
      tnlocal(nsame,2) = 'N3'
      tnlocal(nsame,3) = 'C4'
      tnlocal(nsame,4) = 'C5'
      call twh_set_tornames(type,nsame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffshah2004_bond_increment(lfailure,ffname,nbitypes)
c     ******************************************************************
c     * all of the charges are assigned in the qbasevalue in nonbond so*
c     * the only thing here is a global wildcard 0.0                   *
c     *                                                                *
c     * originally written 04-24-2006 by M.G. Martin                   *
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
