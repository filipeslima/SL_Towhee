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
      subroutine twh_ffopls1996(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the OPLS-1996    *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 06-29-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      character*(FFNAMELEN) ffname

      ffname = 'OPLS-1996'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Geometric')

      call twh_ffopls1996_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffopls1996_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffopls1996_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffopls1996_dihedral(lfailure,ffname,ndihedtypes)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffopls1996_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * sets up the nonbonded portion of the OPLS-1996 forcefield      *
c     *                                                                *
c     * originally written 06-29-2006 by M.G. Martin                   *
c     * last modified 10-14-2008 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      character*(FFNAMELEN) ffname
      integer ntypes
c     --- local variables
      integer type,ii,atomnum
      double precision charge
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      type = 0

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.570d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.370d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.910d0
      nbclocal(2) = (0.160d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.850d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-10'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.175d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CD-11'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.110d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -12'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.850d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -13'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.425d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-14'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.285d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-15'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.285d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-16'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -17'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O2-18'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.800d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-19'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.310d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3-20'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-21'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.330d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-22'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-23'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.070d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-24'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-25'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.850d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -26'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.110d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-27'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.310d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-28'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-29'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.310d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-30'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-31'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.180d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'SH-32'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.450d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HS-33'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.270d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-34'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.235d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S -35'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.470d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-36'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.235d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-37'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S -38'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-39'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-40'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.570d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -41'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.420d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NB-42'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.490d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CP-43'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.410d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CF-44'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CC-45'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.130d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-46'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.540d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -47'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.460d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CP-48'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CG-49'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.330d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CB-50'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.055d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-51'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.800d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-52'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.460d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-53'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.640d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.250d0
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-54'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-55'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.440d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-56'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.310d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-57'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.070d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -58'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.550d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -59'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.450d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-60'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.250d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-61'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.250d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OS-62'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.400d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.000d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-63'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.250d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-64'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-65'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.960d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C4-66'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.730d0
      nbclocal(2) = (0.294d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-67'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = (0.207d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-68'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.175d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-69'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.910d0
      nbclocal(2) = (0.160d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-70'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.960d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-71'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C9-72'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.850d0
      nbclocal(2) = (0.140d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-73'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.850d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.115d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CD-75'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.110d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-76'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C7-77'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-78'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.070d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-79'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-80'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = (0.207d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-81'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'SH-82'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.470d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.700d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'SH-83'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.450d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S -84'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.470d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S -85'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HS-86'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.235d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HS-87'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.270d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-88'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.180d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = (0.207d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-89'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.180d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.905d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-90'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.235d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-91'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.235d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-92'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-93'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-96'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.150d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = (0.207d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ar-97'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 18
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.401d0
      nbclocal(2) = (0.2339d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Kr-98'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 36
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.624d0
      nbclocal(2) = (0.3170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Xe-99'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 54
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.935d0
      nbclocal(2) = (0.4330d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3-101'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.40d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3-102'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.30d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-104'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.35d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-105'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.33d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-106'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.31d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.775d0
      nbclocal(2) = (0.207d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-107'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.25d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.960d0
      nbclocal(2) = (0.145d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OS-108'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.50d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.000d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-109'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.25d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-110'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.25d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2-118'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.118d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl-119'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 17
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.250d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0
      nbclocal(2) = (0.300d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH-120'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.420d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.080d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl-121'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 17
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.140d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.470d0
      nbclocal(2) = (0.300d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-122'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.248d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl-123'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 17
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.062d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.470d0
      nbclocal(2) = (0.266d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ne-129'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 10
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.780d0
      nbclocal(2) = (0.069d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'He-130'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 2
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.556d0
      nbclocal(2) = (0.020d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -131'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.800d0
      nbclocal(2) = (0.115d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C3-132'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.285d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.80d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-135'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.18d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-136'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.12d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-137'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.06d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-138'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.24d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-139'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-140'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.06d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-141'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.076d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-142'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.115d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.076d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-143'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.230d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.076d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-144'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.115d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.42d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-145'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.115d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HA-146'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.115d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.42d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-147'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-148'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.065d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-149'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.005d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-154'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.683d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.120d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-155'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.418d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-156'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.040d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-157'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.145d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-158'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.205d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-159'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-160'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.1263d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-161'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.5323d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.062d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-162'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.6351d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.07d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-163'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.4286d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F -164'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 9
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.2057d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.94d0
      nbclocal(2) = (0.061d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-165'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0825d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-166'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.150d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-167'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.585d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.07d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-168'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-169'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.070d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-170'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OS-180'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.400d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.90d0
      nbclocal(2) = (0.140d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-181'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.110d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-182'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.140d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-183'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.170d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-184'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-185'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.030d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OS-186'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.400d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.00d0
      nbclocal(2) = (0.140d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-187'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.070d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-188'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO-189'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-190'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO-191'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.265d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-192'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO-193'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-194'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO-195'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.365d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-196'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO-197'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.400d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CO-198'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.465d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'SH-200'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'SH-201'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.470d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.700d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S -202'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.435d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'S -203'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 16
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.2175d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HS-204'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.255d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HS-205'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.235d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-206'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.060d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-207'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.120d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-208'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.180d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-209'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0375d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-210'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0975d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-211'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.1575d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-212'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.2175d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-213'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0375d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-214'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0975d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-215'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.1575d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-216'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.2175d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NT-220'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.900d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-221'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.350d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-222'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.020d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-223'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.080d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-224'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.140d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-225'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-226'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.55d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NT-227'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.900d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-228'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.400d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NT-229'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.660d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-230'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.66d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-231'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.040d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -235'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -236'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -237'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.760d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -238'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -239'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.140d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -240'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.380d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -241'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-242'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.020d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-243'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.110d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-244'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.080d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-245'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.050d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-246'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.010d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -247'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.142d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -248'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.390d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -249'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.542d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -250'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.333d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N -251'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.490d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -252'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.420d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -253'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.420d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -254'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.370d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-255'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.060d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.500d0
      nbclocal(2) = (0.020d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-256'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.120d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-257'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.060d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-258'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-259'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.060d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -267'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.520d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OH-268'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.530d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.000d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -269'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.440d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-270'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.450d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -271'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O2-272'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.800d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-273'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.280d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-274'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.220d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-275'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.160d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-276'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.100d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -277'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.450d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -278'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.450d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-279'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.420d0
      nbclocal(2) = (0.015d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -280'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.470d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.750d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -281'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.470d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-282'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.060d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.420d0
      nbclocal(2) = (0.015d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3-286'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.40d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N3-287'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.30d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-289'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.35d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-290'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.33d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-291'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.130d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-292'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.190d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-293'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.250d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-294'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.310d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-300'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.800d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-301'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.460d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-302'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.640d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.250d0
      nbclocal(2) = (0.050d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-303'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.700d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3-304'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.440d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-305'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.200d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-306'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.110d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-307'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.190d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-308'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.050d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-311'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.46d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-312'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.36d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-313'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.85d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-314'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.37d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-315'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.15d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-316'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-317'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.04d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-318'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-319'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.60d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -320'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.50d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-321'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.51d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -322'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.45d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-323'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.07d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-324'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.08d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -325'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.41d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -326'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.40d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -327'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.36d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -328'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.42d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-329'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-330'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-331'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.14d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-332'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.08d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-333'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.56d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -334'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.55d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NC-335'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.54d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-336'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.46d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-337'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.06d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-338'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -339'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.38d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -340'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.48d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96d0
      nbclocal(2) = (0.21d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-341'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.79d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-342'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.385d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-343'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.355d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-344'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-345'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NC-346'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.53d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CQ-347'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.22d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NC-348'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.55d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CB-349'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.38d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CB-350'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.15d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-351'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.44d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NB-352'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.49d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CK-353'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.20d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-354'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.50d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-355'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.20d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-356'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.81d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-357'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.385d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-358'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.355d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-359'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.20d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -360'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.35d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-361'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.56d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-362'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.46d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NC-363'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.51d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CB-364'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.34d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CB-365'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.12d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -366'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.52d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -367'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.38d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-368'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.80d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-369'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.40d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -370'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.51d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96d0
      nbclocal(2) = (0.21d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-371'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.01d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-372'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.12d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-373'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.01d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-374'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.14d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-375'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.01d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-376'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.13d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-377'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.64d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C -378'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.65d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.75d0
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-379'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.74d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-380'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.66d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-381'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.06d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CM-382'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.10d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -383'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.49d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -384'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.30d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.96d0
      nbclocal(2) = (0.21d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -385'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.48d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2-386'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.81d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.25d0
      nbclocal(2) = (0.17d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-387'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.46d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2-388'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.43d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-389'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.14d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-390'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.14d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-391'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.01d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.50d0
      nbclocal(2) = (0.08d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HC-392'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.16d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.50d0
      nbclocal(2) = (0.05d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'P -393'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 15
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.780d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.740d0
      nbclocal(2) = (0.200d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O2-393'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.660d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.960d0
      nbclocal(2) = (0.210d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OS-395'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.430d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.000d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-396'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.020d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F -400'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 9
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -1.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.73295d0
      nbclocal(2) = (0.72000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl-401'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 17
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -1.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.41724d0
      nbclocal(2) = (0.11779d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Br-402'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 35
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -1.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.62376d0
      nbclocal(2) = (0.09000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Li-404'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 3
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.25992d0
      nbclocal(2) = (6.25000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Na-405'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 11
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.0d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.89744d0
      nbclocal(2) = (1.60714d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Li-406'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 3
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.126452d0
      nbclocal(2) = (0.018279d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Na-407'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 11
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.330445d0
      nbclocal(2) = (0.002772d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'K -408'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 19
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.934628d0
      nbclocal(2) = (0.000328d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Rb-409'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 37
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 5.621773d0
      nbclocal(2) = (0.000171d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cs-410'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 55
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 1.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 6.715999d0
      nbclocal(2) = (0.000081d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Mg-411'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 12
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 2.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.644471d0
      nbclocal(2) = (0.875044d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ca-412'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 20
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 2.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.412031d0
      nbclocal(2) = (0.449657d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Sr-413'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 38
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 2.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.102688d0
      nbclocal(2) = (0.118226d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ba-414'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 56
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 2.00d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.816610d0
      nbclocal(2) = (0.047096d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O -434'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -1.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.200d0
      nbclocal(2) = (0.250d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'HO-435'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.300d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'U -436'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 92
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 2.500d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.81524d0
      nbclocal(2) = (0.400d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'OU-437'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 8
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.250d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.11815d0
      nbclocal(2) = (0.200d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-500'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.075d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-501'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.055d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-502'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.130d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-503'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.570d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -504'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.420d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CT-505'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.005d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.500d0
      nbclocal(2) = (0.066d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-506'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.295d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-507'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.015d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-508'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.015d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-509'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.385d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CA-510'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 6
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.215d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.550d0
      nbclocal(2) = (0.070d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NB-511'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.490d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'NA-512'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 7
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = -0.540d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.250d0
      nbclocal(2) = (0.170d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H -513'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      atomnum = 1
      call twh_getelement(lfailure,'get e and m',masslocal,elelocal
     &     ,atomnum)
      if ( lfailure ) return
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      charge = 0.460d0
      call twh_qbasevalue(GLB_SET,type,charge)
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffopls1996_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * sets up the bond potentials for OPLS-1996                      *
c     *                                                                *
c     * originally written 06-30-2006 by M.G. Martin                   *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nbondtypes
c     --- local variables
      integer type
      integer nvslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)

c     --- initialize type count
      type = 0

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.9572d0
      vclocal(1) = (529.60d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- OW-HW
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'OW-wild'
      vnlocal(1,2) = 'HW-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.5139d0
      vclocal(1) = (38.25d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- HW-HW
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'HW-wild'
      vnlocal(1,2) = 'HW-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-HC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'HC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.419d0
      vclocal(1) = (447.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CD
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CD-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.444d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CJ
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CJ-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.444d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CM
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CM-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.522d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (490.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.383d0
      vclocal(1) = (424.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.388d0
      vclocal(1) = (418.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.358d0
      vclocal(1) = (457.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -NC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'NC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.229d0
      vclocal(1) = (570.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -O
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'O -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.25d0
      vclocal(1) = (656.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -O2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'O2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.364d0
      vclocal(1) = (450.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -OH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'OH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.327d0
      vclocal(1) = (214.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C -OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C -wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.495d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.459d0
      vclocal(1) = (388.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-CB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'CB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.352d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-CG
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'CG-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.352d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-CC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'CC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.495d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.352d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C*-CW
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C*-wild'
      vnlocal(1,2) = 'CW-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.382d0
      vclocal(1) = (546.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CW
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CW-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-CA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'CA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-CD
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'CD-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.504d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-CC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'CC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-CH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'CH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.463d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-N3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'N3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.448d0
      vclocal(1) = (382.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.448d0
      vclocal(1) = (382.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.448d0
      vclocal(1) = (382.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.448d0
      vclocal(1) = (382.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (386.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-OH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'OH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CO-OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CO-wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CO-C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CO-wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CO-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CO-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-S
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'S -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-SH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'SH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-CH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'CH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-CM
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'CM-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.463d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-N3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'N3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (386.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-OH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'OH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C3-S
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C3-wild'
      vnlocal(1,2) = 'S -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C2-SH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C2-wild'
      vnlocal(1,2) = 'SH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.404d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CD
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CD-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.433d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CJ
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CJ-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.433d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CM
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CM-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CN
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CN-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (481.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.339d0
      vclocal(1) = (483.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-NC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'NC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.320d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- NC-NC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'NC-wild'
      vnlocal(1,2) = 'NC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.080d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-HA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'HA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.370d0
      vclocal(1) = (520.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-CB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'CB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.370d0
      vclocal(1) = (520.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CV-CW
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CV-wild'
      vnlocal(1,2) = 'CW-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-CD
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'CD-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.419d0
      vclocal(1) = (447.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-CN
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'CN-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.374d0
      vclocal(1) = (436.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.391d0
      vclocal(1) = (414.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.354d0
      vclocal(1) = (461.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-NC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'NC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (512.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-CF
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'CF-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (518.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-CG
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'CG-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.504d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.375d0
      vclocal(1) = (512.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-CV
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'CV-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (518.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-CW
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'CW-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.385d0
      vclocal(1) = (422.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.394d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CD-CD
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CD-wild'
      vnlocal(1,2) = 'CD-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CD-CN
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CD-wild'
      vnlocal(1,2) = 'CN-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-CC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'CC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.40d0
      vclocal(1) = (469.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CD-CC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CD-wild'
      vnlocal(1,2) = 'CC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CE-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CE-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.304d0
      vclocal(1) = (529.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CE-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CE-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.394d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CF-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CF-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CG-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CG-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-CH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'CH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (386.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-OH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'OH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.425d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CH-OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CH-wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.324d0
      vclocal(1) = (502.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CI-NC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CI-wild'
      vnlocal(1,2) = 'NC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.350d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CJ-CJ
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CJ-wild'
      vnlocal(1,2) = 'CJ-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.350d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CJ-CM
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CJ-wild'
      vnlocal(1,2) = 'CM-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (448.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CJ-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CJ-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CK-HA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CK-wild'
      vnlocal(1,2) = 'HA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CK-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CK-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.371d0
      vclocal(1) = (440.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CK-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CK-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.304d0
      vclocal(1) = (529.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CK-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CK-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (549.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CM-CM
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CM-wild'
      vnlocal(1,2) = 'CM-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.460d0
      vclocal(1) = (532.8d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-C+
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'C+-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.084d0
      vclocal(1) = (532.8d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C+-HC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C+-wild'
      vnlocal(1,2) = 'HC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CM-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CM-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CM-HC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CM-wild'
      vnlocal(1,2) = 'HC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.09d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- HC-C
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'HC-wild'
      vnlocal(1,2) = 'C -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (448.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CM-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CM-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.365d0
      vclocal(1) = (448.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CM-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CM-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.38d0
      vclocal(1) = (428.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CN-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CN-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (447.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CP-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CP-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (488.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CP-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CP-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CQ-HA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CQ-wild'
      vnlocal(1,2) = 'HA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CQ-H5
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CQ-wild'
      vnlocal(1,2) = 'H5-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.324d0
      vclocal(1) = (502.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CQ-NC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CQ-wild'
      vnlocal(1,2) = 'NC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CR-HA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CR-wild'
      vnlocal(1,2) = 'HA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CR-H5
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CR-wild'
      vnlocal(1,2) = 'H5-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.343d0
      vclocal(1) = (477.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CR-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CR-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.335d0
      vclocal(1) = (488.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CR-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CR-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.529d0
      vclocal(1) = (268.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.09d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-HC
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'HC-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.449d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-N*
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'N*-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.475d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.463d0
      vclocal(1) = (337.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-N3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'N3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.471d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.41d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-OH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'OH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.41d0
      vclocal(1) = (320.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-S
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'S -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-SH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'SH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.781d0
      vclocal(1) = (245.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-Cl
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'Cl-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CV-HA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CV-wild'
      vnlocal(1,2) = 'HA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CV-H4
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CV-wild'
      vnlocal(1,2) = 'H4-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.394d0
      vclocal(1) = (410.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CV-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CV-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CW-H4
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CW-wild'
      vnlocal(1,2) = 'H4-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.08d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CW-HA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CW-wild'
      vnlocal(1,2) = 'HA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CW-NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CW-wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H -N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H -wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H -N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H -wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H -NA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H -wild'
      vnlocal(1,2) = 'NA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H2-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H2-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H2-N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H2-wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H2-NT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H2-wild'
      vnlocal(1,2) = 'NT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H3-N2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H3-wild'
      vnlocal(1,2) = 'N2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.01d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- H3-N3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'H3-wild'
      vnlocal(1,2) = 'N3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.945d0
      vclocal(1) = (553.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- HO-OH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'HO-wild'
      vnlocal(1,2) = 'OH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 0.945d0
      vclocal(1) = (553.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- HO-OS
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'HO-wild'
      vnlocal(1,2) = 'OS-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.336d0
      vclocal(1) = (274.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- HS-SH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'HS-wild'
      vnlocal(1,2) = 'SH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.48d0
      vclocal(1) = (525.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- O2-P
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'O2-wild'
      vnlocal(1,2) = 'P -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.61d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- OH-P
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'OH-wild'
      vnlocal(1,2) = 'P -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.61d0
      vclocal(1) = (230.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- OS-P
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'OS-wild'
      vnlocal(1,2) = 'P -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.843d0
      vclocal(1) = (212.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-P
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'P -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 2.038d0
      vclocal(1) = (166.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- S -S
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'S -wild'
      vnlocal(1,2) = 'S -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (530.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C9-C9
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C9-wild'
      vnlocal(1,2) = 'C9-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (530.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C9-C8
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C9-wild'
      vnlocal(1,2) = 'C8-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (530.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C9-C7
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C9-wild'
      vnlocal(1,2) = 'C7-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (530.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C8-C8
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C8-wild'
      vnlocal(1,2) = 'C8-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (530.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C8-C7
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C8-wild'
      vnlocal(1,2) = 'C7-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.34d0
      vclocal(1) = (530.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C7-C7
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C7-wild'
      vnlocal(1,2) = 'C7-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C8-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C8-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C8-CH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C8-wild'
      vnlocal(1,2) = 'CH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C8-C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C8-wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C8-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C8-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C7-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C7-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C7-CH
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C7-wild'
      vnlocal(1,2) = 'CH-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C7-C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C7-wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.50d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- C7-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'C7-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (260.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-C2
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'C2-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.391d0
      vclocal(1) = (414.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-NB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'NB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.381d0
      vclocal(1) = (427.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CA-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CA-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.51d0
      vclocal(1) = (317.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CB-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CB-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.370d0
      vclocal(1) = (520.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CC-CB
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CC-wild'
      vnlocal(1,2) = 'CB-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.332d0
      vclocal(1) = (367.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- CT-F
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'CT-wild'
      vnlocal(1,2) = 'F -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.81d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- SY-C3
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'SY-wild'
      vnlocal(1,2) = 'C3-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.77d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- SY-CA
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'SY-wild'
      vnlocal(1,2) = 'CA-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.76d0
      vclocal(1) = (222.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- SY-C8
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'SY-wild'
      vnlocal(1,2) = 'C8-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.44d0
      vclocal(1) = (700.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- SY-OY
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'SY-wild'
      vnlocal(1,2) = 'OY-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.67d0
      vclocal(1) = (434.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- SY-N
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'SY-wild'
      vnlocal(1,2) = 'N -wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.77d0
      vclocal(1) = (340.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- SY-CT
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'SY-wild'
      vnlocal(1,2) = 'CT-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.80d0
      vclocal(1) = (500.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 0
c     --- U -OU
      nvslocal = nvslocal + 1
      vnlocal(1,1) = 'U -wild'
      vnlocal(1,2) = 'OU-wild'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffopls1996_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * sets up the angle terms for OPLS-1996                          *
c     *                                                                *
c     * originally written 07-01-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nangletypes
c     --- local variables
      integer type,isame,nbslocal
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)

c     --- initialize type count
      type = 0

      type = type + 1
c     --- standard harmonic
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 150.0d0
      bclocal(1) = (180.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
C     --- OU-U OU
      isame = isame + 1
      bnlocal(isame,1) = 'OU-wild'
      bnlocal(isame,2) = 'U -wild'
      bnlocal(isame,3) = 'OU-wild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffopls1996_dihedral(lfailure,ffname,ndihedtypes)
c     ******************************************************************
c     * OPLS-1996 dihedral                                             *
c     *                                                                *
c     * originally written 07-01-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ndihedtypes
c     --- local variables
      integer type,isame
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision onepi

      call twh_constant_pi(GLB_GET,onepi)

c     --- initialize type count
      type = 0

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      tclocal(1) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_KCALTOK)*(0.5d0)
      tclocal(3) =  (0.318d0)*(CON_KCALTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame = 0
c     --- H-C-C-H hydrocarbon all-atom
      isame = isame + 1
      tnlocal(isame,1) = 'HC-wild'
      tnlocal(isame,2) = 'CT-wild'
      tnlocal(isame,3) = 'CT-wild'
      tnlocal(isame,4) = 'HC-wild'
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
