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
      subroutine twh_ffpana1989(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for pana1989                                *
c     *                                                                *
c     * originally written 03-02-2006 by M.G. Martin                   *
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
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Pana1989'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'wild')
      type = 0
c     --- references
c     A. Z. Panagiotopoulos
c     "Exact Calculations of Fluid-Phase Equilibria by Monte Carlo
c     Simulation in a New Statistical Ensemble"
c     Int. J. Thermophys. 10, 447-457 (1989)
c     --- polarizabilities are the entries marked DC from Table D.3
c     --- page 577 of
c     "Theory of molecular fluids Volume 1: Fundamentals"
c     C.G. Gray and K.E. Gubbins Clarendon Press Oxford 1984.

c     --- Neon
      type = type + 1
      masslocal = 20.180d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'Ne'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_polarizability(GLB_SET,type,0.396d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ne'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.81d0
      nbclocal(2) = 34.3d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- Argon
      type = type + 1
      masslocal = 39.948d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = 'Ar'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_polarizability(GLB_SET,type,1.642d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Ar'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.390d0
      nbclocal(2) = 117.5d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- Krypton
      type = type + 1
      masslocal = 83.80d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'Kr'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_polarizability(GLB_SET,type,2.484d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Kr'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.607d0
      nbclocal(2) = 161.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- Xenon
      type = type + 1
      masslocal = 131.29d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal =  'Xe'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_polarizability(GLB_SET,type,4.02d0)
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Xe'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.91d0
      nbclocal(2) = 227.0d0
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
