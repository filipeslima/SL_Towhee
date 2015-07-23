#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2006-2011 Marcus G. Martin                       *
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
      subroutine twh_ffshukla1987(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for shukla1987                              *
c     *                                                                *
c     * originally written 04-18-2006 by M.G. Martin                   *
c     * last modified 08-04-2011 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Shukla1987'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'wild')
      type = 0
c     --- reference
c     K.P. Shukla "Thermodynamic properties of simple fluid mixtures
c     from perturbation theory", Mol. Phys. 62, 5, 1143-1163 (1987)
c     Parameters implemented here are listed in Table 6 of that work

c     --- Argon
      type = type + 1
      call twh_mass(GLB_SET,type,39.948d0) 
      call twh_elementname(GLB_SET,type,'Ar')
      call twh_polarizability(GLB_SET,type,1.642d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ar'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.404d0
      nbclocal(2) = 117.70d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- Krypton
      type = type + 1
      call twh_mass(GLB_SET,type,83.80d0)
      call twh_elementname(GLB_SET,type,'Kr')
      call twh_polarizability(GLB_SET,type,2.484d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Kr'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.622d0
      nbclocal(2) = 164.11d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- Xenon
      type = type + 1
      call twh_mass(GLB_SET,type,131.29d0)
      call twh_elementname(GLB_SET,type,'Xe')
      call twh_polarizability(GLB_SET,type,4.017d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Xe'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.948d0
      nbclocal(2) = 227.86d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- methane
      type = type + 1
      masslocal = 1.0d0*(12.011d0) + 4.0d0*(1.0079d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' C')
      call twh_polarizability(GLB_SET,type,2.60d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH4'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.733d0
      nbclocal(2) = 149.92d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- nitrogen dimer
      type = type + 1
      masslocal = 2.0d0*(14.007d0)
      call twh_mass(GLB_SET,type,masslocal)
      call twh_elementname(GLB_SET,type,' N')
      call twh_polarizability(GLB_SET,type,1.740d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N2'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.607d0
      nbclocal(2) = 100.31d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
