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
      subroutine twh_ffmorrow2002(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes,nimprotypes,nbitypes)
c     ******************************************************************
c     * sets up parameters for morrow2002                              *
c     *                                                                *
c     * originally written 04-19-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
     &        ,nimprotypes,nbitypes
c     --- local variables
      character*(FFNAMELEN) ffname

c     --- references
c     T.I. Morrow; E.J. Maginn;
c     "Molecular Dynamics Study of the Ionic Liquid 1-n-Butyl-3-
c     methylimidazolium Hexafluorophosphate" J. Phys. Chem. B
c     106, 12807-12813 (2002).
c     T.I. Morrow; E.J. Maginn;
c     "Erratum 2002, Vol. 106B" J. Phys. Chem. B 107, 9160 (2003)

      ffname = 'Morrow2002'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

      call twh_ffmorrow2002_nonbond(lfailure,ffname,ntypes)
      if ( lfailure ) return
      call twh_ffmorrow2002_bond(lfailure,ffname,nbondtypes)
      if ( lfailure ) return
      call twh_ffmorrow2002_angle(lfailure,ffname,nangletypes)
      if ( lfailure ) return
      call twh_ffmorrow2002_dihedral(lfailure,ffname,ndihedtypes)
      if ( lfailure ) return
      call twh_ffmorrow2002_improper(lfailure,ffname,nimprotypes)
      if ( lfailure ) return
      call twh_ffmorrow2002_bond_increment(lfailure,ffname,nbitypes)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffmorrow2002_nonbond(lfailure,ffname,ntypes)
c     ******************************************************************
c     * sets up the nonbonded parameters for morrow2002                *
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

c     --- full vdw nonbonded interactions are included in the 1-4 terms
c     --- personal communication from E.J. Maginn to M.G. Martin on 
c     --- 04-21-2006.

      type = 0

c     --- C2
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.056d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.600d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.209d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C4
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.141d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.600d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.209d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C5
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.217d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.600d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.209d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C6
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.157d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.550d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.084d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C7
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.095d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.550d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.084d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C8
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.122d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.020d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.234d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C9
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,0.256d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.020d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.234d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- C10
      type = type + 1
      masslocal = 12.011d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_qbasevalue(GLB_SET,type,-0.209d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C10'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.080d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.326d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F1
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F2
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F3
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F4
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F5
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- F6
      type = type + 1
      masslocal = 18.998d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' F')
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.400d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.377d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H1
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.177d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.800d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.192d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H2
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.181d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.936d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.033d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H3
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.207d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.936d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.033d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H4
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.125d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.640d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.092d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H5
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.073d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.640d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.092d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H6
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.142d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.640d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.092d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H7
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.055d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.092d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H8
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.045d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H8'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.092d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H9
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.055d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H9'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.117d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H10
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.001d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H10'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.117d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H11
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,-0.029d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H11'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.117d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H12
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,-0.099d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H12'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.117d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H13
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.051d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H13'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.100d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H14
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.040d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H14'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.100d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- H15
      type = type + 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' H')
      call twh_qbasevalue(GLB_SET,type,0.075d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H15'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.680d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.100d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- N1
      type = type + 1
      masslocal = 14.007d0
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
      nbclocal(1) = 3.700d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.837d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- N3
      type = type + 1
      masslocal = 14.007d0
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
      nbclocal(1) = 3.700d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 0.837d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- P1
      type = type + 1
      masslocal = 30.974d0
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' P')
      call twh_qbasevalue(GLB_SET,type,-0.904d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'P'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.30d0/((2.0d0)**(1.0d0/6.0d0))
      nbclocal(2) = 2.448d0*CON_KJTOK
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffmorrow2002_bond(lfailure,ffname,nbondtypes)
c     ******************************************************************
c     * bond lengths for the morrow2002 forcefield                     *
c     * these are not explicitly specified and instead must be         *
c     * inferred from their example structure                          *
c     *                                                                *
c     * originally written 04-19-2006 by M.G. Martin                   *
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
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.470d0
      vclocal(1) = 921.1d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6 - N3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C6'
      vnlocal(nsame,2) = 'N3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.483d0
      vclocal(1) = 921.1d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7 - N1
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C7'
      vnlocal(nsame,2) = 'N1'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.382d0
      vclocal(1) = 1674.7d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C5 - N1
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C5'
      vnlocal(nsame,2) = 'N1'
c     --- C4 - N3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C4'
      vnlocal(nsame,2) = 'N3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.337d0
      vclocal(1) = 1674.7d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2 - N1
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C2'
      vnlocal(nsame,2) = 'N1'
c     --- C2 - N3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C2'
      vnlocal(nsame,2) = 'N3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.361d0
      vclocal(1) = 1716.6d0*CON_KJTOK
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
      vclocal(0) = 1.078d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2 - H1 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C2'
      vnlocal(nsame,2) = 'H1'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.078d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4 - H2 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C4'
      vnlocal(nsame,2) = 'H2'
c     --- C5 - H3 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C5'
      vnlocal(nsame,2) = 'H3'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.089d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6 - H4 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C6'
      vnlocal(nsame,2) = 'H4'
c     --- C6 - H5 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C6'
      vnlocal(nsame,2) = 'H5'
c     --- C6 - H6 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C6'
      vnlocal(nsame,2) = 'H6'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.091d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7 - H7 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C7'
      vnlocal(nsame,2) = 'H7'
c     --- C7 - H8 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C7'
      vnlocal(nsame,2) = 'H8'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.096d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8 - H9 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C8'
      vnlocal(nsame,2) = 'H9'
c     --- C8 - H10 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C8'
      vnlocal(nsame,2) = 'H10'
c     --- C9 - H11 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C9'
      vnlocal(nsame,2) = 'H11'
c     --- C9 - H12 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C9'
      vnlocal(nsame,2) = 'H12'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.093d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C10 - H13 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C10'
      vnlocal(nsame,2) = 'H13'
c     --- C10 - H14 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C10'
      vnlocal(nsame,2) = 'H14'
c     --- C10 - H15 (all H bonds were held rigid)
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C10'
      vnlocal(nsame,2) = 'H15'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.530d0
      vclocal(1) = 837.4d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7 - C8
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C7'
      vnlocal(nsame,2) = 'C8'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.534d0
      vclocal(1) = 931.6d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8 - C9
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C8'
      vnlocal(nsame,2) = 'C9'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.530d0
      vclocal(1) = 931.6d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- C9 - C10
      nsame = nsame + 1
      vnlocal(nsame,1) = 'C9'
      vnlocal(nsame,2) = 'C10'
      call twh_set_vibnames(type,nsame,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.646d0
      vclocal(1) = 1089.8d0*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nsame = 0
c     --- P - F
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F'
c     --- P - F1
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F1'
c     --- P - F2
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F2'
c     --- P - F3
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F3'
c     --- P - F4
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F4'
c     --- P - F5
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F5'
c     --- P - F6
      nsame = nsame + 1
      vnlocal(nsame,1) = 'P'
      vnlocal(nsame,2) = 'F6'
      call twh_set_vibnames(type,nsame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffmorrow2002_angle(lfailure,ffname,nangletypes)
c     ******************************************************************
c     * sets up the angle terms for Morrow2002                         *
c     *                                                                *
c     * originally written 04-19-2006 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
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
      bclocal(0) = 112.6d0
      bclocal(1) = 586.2d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8-C7-N1
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C8'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'N1'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 108.3d0
      bclocal(1) = 544.3d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C5,4-N1,3-C2
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C5'
      bnlocal(nsame,2) = 'N1'
      bnlocal(nsame,3) = 'C2'
c     --- C5,4-N1,3-C2
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C4'
      bnlocal(nsame,2) = 'N3'
      bnlocal(nsame,3) = 'C2'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.9d0
      bclocal(1) = 544.3d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7-N1-C2,5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C7'
      bnlocal(nsame,2) = 'N1'
      bnlocal(nsame,3) = 'C2'
c     --- C7-N1-C2,5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C7'
      bnlocal(nsame,2) = 'N1'
      bnlocal(nsame,3) = 'C5'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.8d0
      bclocal(1) = 544.3d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C6-N3-C2,4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C6'
      bnlocal(nsame,2) = 'N3'
      bnlocal(nsame,3) = 'C2'
c     --- C6-N3-C2,4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C6'
      bnlocal(nsame,2) = 'N3'
      bnlocal(nsame,3) = 'C4'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.6d0
      bclocal(1) = 125.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H4,5,6-C6-N3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H4'
      bnlocal(nsame,2) = 'C6'
      bnlocal(nsame,3) = 'N3'
c     --- H4,5,6-C6-N3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H5'
      bnlocal(nsame,2) = 'C6'
      bnlocal(nsame,3) = 'N3'
c     --- H4,5,6-C6-N3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H6'
      bnlocal(nsame,2) = 'C6'
      bnlocal(nsame,3) = 'N3'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 125.5d0
      bclocal(1) = 104.7d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H1-C2-N1,3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H1'
      bnlocal(nsame,2) = 'C2'
      bnlocal(nsame,3) = 'N1'
c     --- H1-C2-N1,3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H1'
      bnlocal(nsame,2) = 'C2'
      bnlocal(nsame,3) = 'N3'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.2d0
      bclocal(1) = 544.3d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1,3-C5,4-C4,5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N1'
      bnlocal(nsame,2) = 'C5'
      bnlocal(nsame,3) = 'C4'
c     --- N1,3-C5,4-C4,5
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N3'
      bnlocal(nsame,2) = 'C4'
      bnlocal(nsame,3) = 'C5'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.1d0
      bclocal(1) = 544.3d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1-C2-N3
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N1'
      bnlocal(nsame,2) = 'C2'
      bnlocal(nsame,3) = 'N3'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 130.8d0
      bclocal(1) = 104.7d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H2,3-C4,5-C5,4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H2'
      bnlocal(nsame,2) = 'C4'
      bnlocal(nsame,3) = 'C5'
c     --- H2,3-C4,5-C5,4
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H3'
      bnlocal(nsame,2) = 'C5'
      bnlocal(nsame,3) = 'C4'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.3d0
      bclocal(1) = 148.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H4,5,6-C6-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H4'
      bnlocal(nsame,2) = 'C6'
      bnlocal(nsame,3) = 'Hwild'
c     --- H4,5,6-C6-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H5'
      bnlocal(nsame,2) = 'C6'
      bnlocal(nsame,3) = 'Hwild'
c     --- H4,5,6-C6-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H6'
      bnlocal(nsame,2) = 'C6'
      bnlocal(nsame,3) = 'Hwild'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 122.1d0
      bclocal(1) = 104.7d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1,3-C5-4-H3,2
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N1'
      bnlocal(nsame,2) = 'C5'
      bnlocal(nsame,3) = 'H3'
c     --- N1,3-C5-4-H3,2
      nsame = nsame + 1
      bnlocal(nsame,1) = 'N3'
      bnlocal(nsame,2) = 'C4'
      bnlocal(nsame,3) = 'H2'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.2d0
      bclocal(1) = 148.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H7-C7-H8
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H7'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'H8'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.5d0
      bclocal(1) = 139.8d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H7,8-C7-C8
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H7'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'C8'
c     --- H7,8-C7-C8
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H8'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'C8'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = 139.8d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H9,10-C8-C7
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H9'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'C7'
c     --- H9,10-C8-C7
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H10'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'C7'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 111.6d0
      bclocal(1) = 244.5d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7,8-C8,9-C9,10
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C7'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'C9'
c     --- C7,8-C8,9-C9,10
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C8'
      bnlocal(nsame,2) = 'C9'
      bnlocal(nsame,3) = 'C10'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.4d0
      bclocal(1) = 144.4d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H-C8,9-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Hwild'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'Hwild'
c     --- H-C8,9-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Hwild'
      bnlocal(nsame,2) = 'C9'
      bnlocal(nsame,3) = 'Hwild'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.7d0
      bclocal(1) = 144.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8,9-C9,10-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C8'
      bnlocal(nsame,2) = 'C9'
      bnlocal(nsame,3) = 'Hwild'
c     --- C8,9-C9,10-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'C9'
      bnlocal(nsame,2) = 'C10'
      bnlocal(nsame,3) = 'Hwild'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.7d0
      bclocal(1) = 144.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H9,10-C8-C9
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H9'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'C9'
c     --- H9,10-C8-C9
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H10'
      bnlocal(nsame,2) = 'C8'
      bnlocal(nsame,3) = 'C9'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.7d0
      bclocal(1) = 144.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H11,12-C9-C10
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H11'
      bnlocal(nsame,2) = 'C9'
      bnlocal(nsame,3) = 'C10'
c     --- H11,12-C9-C10
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H12'
      bnlocal(nsame,2) = 'C9'
      bnlocal(nsame,3) = 'C10'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.6d0
      bclocal(1) = 148.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- H13,14,15-C10-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H13'
      bnlocal(nsame,2) = 'C10'
      bnlocal(nsame,3) = 'Hwild'
c     --- H13,14,15-C10-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H14'
      bnlocal(nsame,2) = 'C10'
      bnlocal(nsame,3) = 'Hwild'
c     --- H13,14,15-C10-H
      nsame = nsame + 1
      bnlocal(nsame,1) = 'H15'
      bnlocal(nsame,2) = 'C10'
      bnlocal(nsame,3) = 'Hwild'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.8d0
      bclocal(1) = 125.6d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- N-C7-H7,8 (incorrectly listed as N-C7-H8,9 in the erratum)
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Nwild'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'H7'
c     --- N-C7-H7,8 (incorrectly listed as N-C7-H8,9 in the erratum)
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Nwild'
      bnlocal(nsame,2) = 'C7'
      bnlocal(nsame,3) = 'H8'
      call twh_set_bendnames(type,nsame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 90.0d0
      bclocal(1) = 812.7d0*(CON_KJTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nsame = 0
c     --- F-P-F
      nsame = nsame + 1
      bnlocal(nsame,1) = 'Fwild'
      bnlocal(nsame,2) = 'P'
      bnlocal(nsame,3) = 'Fwild'
      call twh_set_bendnames(type,nsame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffmorrow2002_dihedral(lfailure,ffname,ndihedtypes)
c     ******************************************************************
c     * Morrow2002 cosine series dihedral                              *
c     *                                                                *
c     * originally written 04-19-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer ndihedtypes
c     --- local variables
      integer type,nsame
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision onepi,degtorad

c     --- coulombic 1-4 interactions are scaled by 0.5
c     --- personal communication between E.J. Maginn and M.G. Martin
c     --- 04-21-2006

      call twh_constant_pi(GLB_GET,onepi)
      degtorad = onepi/180.0d0

c     --- initialize type count
      type = 0

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (58.6d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N1,3-C5,4-C4,5
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C5'
      tnlocal(nsame,4) = 'C4'
c     --- C2-N1,3-C5,4-C4,5
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N3'
      tnlocal(nsame,3) = 'C4'
      tnlocal(nsame,4) = 'C5'
c     --- N1-C5-C4-N3
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'C4'
      tnlocal(nsame,4) = 'N3'
c     --- N1,3-C2-N3,1-C4,5
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C4'
c     --- N1,3-C2-N3,1-C4,5
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N3'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C5'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (12.6d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- H1-C2-N1,3-C5,4
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C5'
c     --- H1-C2-N1,3-C5,4
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C4'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (8.37d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- H2-C4-C5-H3
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H2'
      tnlocal(nsame,2) = 'C4'
      tnlocal(nsame,3) = 'C5'
      tnlocal(nsame,4) = 'H3'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KJTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4,5-C5,4-N1,3-C7,6
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C4'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C7'
c     --- C4,5-C5,4-N1,3-C7,6
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C5'
      tnlocal(nsame,2) = 'C4'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C6'
c     --- N1-C7-C8-C9
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'C9'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C4,5-N3,1-C6,7-H
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C4'
      tnlocal(nsame,2) = 'N3'
      tnlocal(nsame,3) = 'C6'
      tnlocal(nsame,4) = 'Hwild'
c     --- C4,5-N3,1-C6,7-H
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C5'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'Hwild'
c     --- N1-C7-C8-H9,10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'H9'
c     --- N1-C7-C8-H9,10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'H10'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.42d0)*(CON_KJTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (180.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C2-N3,1-C7-C8
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N3'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'C8'
c     --- C2-N3,1-C7-C8
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'C8'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.84d0)*(CON_KJTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C5-N1-C7-C8
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C5'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'C8'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.67d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- H-C9-C10-H
      nsame = nsame + 1
      tnlocal(nsame,1) = 'Hwild'
      tnlocal(nsame,2) = 'C9'
      tnlocal(nsame,3) = 'C10'
      tnlocal(nsame,4) = 'Hwild'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.82d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- H7,8-C7-C8-C9
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H7'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'C9'
c     --- H7,8-C7-C8-C9
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H8'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'C9'
c     --- C2-N1,3-C7,6-H
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N1'
      tnlocal(nsame,3) = 'C7'
      tnlocal(nsame,4) = 'Hwild'
c     --- C2-N1,3-C7,6-H
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C2'
      tnlocal(nsame,2) = 'N3'
      tnlocal(nsame,3) = 'C6'
      tnlocal(nsame,4) = 'Hwild'
c     --- C7-C8-C9-H11,12
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C7'
      tnlocal(nsame,2) = 'C8'
      tnlocal(nsame,3) = 'C9'
      tnlocal(nsame,4) = 'H11'
c     --- C7-C8-C9-H11,12
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C7'
      tnlocal(nsame,2) = 'C8'
      tnlocal(nsame,3) = 'C9'
      tnlocal(nsame,4) = 'H12'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.67d0)*(CON_KJTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C8-C9-C10-H
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C8'
      tnlocal(nsame,2) = 'C9'
      tnlocal(nsame,3) = 'C10'
      tnlocal(nsame,4) = 'Hwild'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (12.6d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- H2,3-C4,5-N3,1-C2
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H2'
      tnlocal(nsame,2) = 'C4'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C2'
c     --- H2,3-C4,5-N3,1-C2
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H3'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C2'
c     --- N1,3-C5,4-C4,5-H2,3
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'C4'
      tnlocal(nsame,4) = 'H2'
c     --- N1,3-C5,4-C4,5-H2,3
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N3'
      tnlocal(nsame,2) = 'C4'
      tnlocal(nsame,3) = 'C5'
      tnlocal(nsame,4) = 'H3'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.0d0)*(CON_KJTOK)
      tclocal(2) = 2.0d0
      tclocal(3) = (180.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1,3-C2-N3,1-C6,7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C6'
c     --- N1,3-C2-N3,1-C6,7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'N3'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C7'
c     --- H1-C2-N3,1-C6,7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C6'
c     --- H1-C2-N3,1-C6,7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H1'
      tnlocal(nsame,2) = 'C2'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C7'
c     --- H2,3-C4,5-N3,1-C6,7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H2'
      tnlocal(nsame,2) = 'C4'
      tnlocal(nsame,3) = 'N3'
      tnlocal(nsame,4) = 'C6'
c     --- H2,3-C4,5-N3,1-C6,7
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H3'
      tnlocal(nsame,2) = 'C5'
      tnlocal(nsame,3) = 'N1'
      tnlocal(nsame,4) = 'C7'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.63d0)*(CON_KJTOK)
      tclocal(2) = 1.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- C7-C8-C9-C10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'C7'
      tnlocal(nsame,2) = 'C8'
      tnlocal(nsame,3) = 'C9'
      tnlocal(nsame,4) = 'C10'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.5d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.82d0)*(CON_KJTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*degtorad
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      nsame = 0
c     --- H7,8-C7-C8-H9,10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H7'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'H9'
c     --- H7,8-C7-C8-H9,10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H7'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'H10'
c     --- H7,8-C7-C8-H9,10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H8'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'H9'
c     --- H7,8-C7-C8-H9,10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'H8'
      tnlocal(nsame,2) = 'C7'
      tnlocal(nsame,3) = 'C8'
      tnlocal(nsame,4) = 'H10'
c     --- H-C8-C9-C10
      nsame = nsame + 1
      tnlocal(nsame,1) = 'Hwild'
      tnlocal(nsame,2) = 'C8'
      tnlocal(nsame,3) = 'C9'
      tnlocal(nsame,4) = 'C10'
c     --- H-C8-C9-H [added by M.G. Martin assuming it is the same as
c     --- the published H-C7-C8-H term
      nsame = nsame + 1
      tnlocal(nsame,1) = 'Hwild'
      tnlocal(nsame,2) = 'C8'
      tnlocal(nsame,3) = 'C9'
      tnlocal(nsame,4) = 'Hwild'
c     --- set ntsame
      call twh_set_tornames(type,nsame,tnlocal)
      
      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end

      subroutine twh_ffmorrow2002_improper(lfailure,ffname,nimprotypes)
c     ******************************************************************
c     * Morrow 2002 improper torsions                                  *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
c     --- variables passed to/from the subroutine
      character*(FFNAMELEN) ffname
      logical lfailure
      integer nimprotypes
c     --- local variables
      integer type,nsame
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)

c     --- initialize type count
      type = 0

c     --- the naming in their paper looks inconsistent, I have
c     --- reordered so the central atom is listed first

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.50d0)*(CON_KJTOK)
      iclocal(1) = (0.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nsame = 0
c     --- H1-N1-N3-C2 (reordering so C2 is the center atom)
      nsame = nsame + 1
      ninlocal(nsame,1) = 'C2'
      ninlocal(nsame,2) = 'N3'
      ninlocal(nsame,3) = 'N1'
      ninlocal(nsame,4) = 'H1'
c     --- H2,3-N3,1-C4,5-C (reordering with C4 in the center)
      nsame = nsame + 1
      ninlocal(nsame,1) = 'C4'
      ninlocal(nsame,2) = 'N3'
      ninlocal(nsame,3) = 'C5'
      ninlocal(nsame,4) = 'H2'
c     --- H2,3-N3,1-C4,5-C (reordering with C5 in the center)
      nsame = nsame + 1
      ninlocal(nsame,1) = 'C5'
      ninlocal(nsame,2) = 'N1'
      ninlocal(nsame,3) = 'C4'
      ninlocal(nsame,4) = 'H3'
      call twh_set_impnames(type,nsame,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,2)
      call twh_impstyle(GLB_SET,type,1)
      iclocal(0) = (0.60d0)*(CON_KJTOK)
      iclocal(1) = (0.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      nsame = 0
c     --- N1,3-C4,5-C2-C6,7
      nsame = nsame + 1
      ninlocal(nsame,1) = 'N1'
      ninlocal(nsame,2) = 'C5'
      ninlocal(nsame,3) = 'C2'
      ninlocal(nsame,4) = 'C7'
c     --- N1,3-C4,5-C2-C6,7
      nsame = nsame + 1
      ninlocal(nsame,1) = 'N3'
      ninlocal(nsame,2) = 'C4'
      ninlocal(nsame,3) = 'C2'
      ninlocal(nsame,4) = 'C6'
      call twh_set_impnames(type,nsame,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end

      subroutine twh_ffmorrow2002_bond_increment(lfailure,ffname
     &     ,nbitypes)
c     ******************************************************************
c     * the bond increments for the Morrow2002 force field             *
c     * mostly set to null has we have directly assigned the charges   *
c     * for bmim in the qbase so a zero default wildcard, but using the*
c     * bond increments for PF6-                                       *
c     *                                                                *
c     * originally written 04-19-2006 by M.G. Martin                   *
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
      qbivlocal = 0.3935d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- generic P-F for the all equivalent model
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.421d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- assymmetric P-F1
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F1'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.426d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- assymmetric P-F2
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F2'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.368d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- assymmetric P-F3
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F3'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.368d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- assymmetric P-F4
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F4'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.364d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- assymmetric P-F5
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F5'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.414d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- assymmetric P-F6
      qbinlocal(1) = 'P'
      qbinlocal(2) = 'F6'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     --- default zero wildcards
      qbinlocal(1) = 'wild'
      qbinlocal(2) = 'wild'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
