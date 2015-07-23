#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2011 Marcus G. Martin                       *
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
      subroutine twh_fflgm(lfailure,ntypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the aqvist       *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
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

      ffname = 'LGM'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     (Br-) Bromine (ionic -1)
c     --- T.P. Lybrand, I. Ghosh, J.A. McCammon
c     `Hydration of chloride and bromide anions: determination
c     of relative free energy by computer simulation'
c     J. Am. Chem. Soc. 107, 7793-7794 (1985)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.5950d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0900d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 79.904d0
      elelocal = 'Br'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Br-'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Cl-) Chlorine (ionic -1)
c     --- T.P. Lybrand, I. Ghosh, J.A. McCammon
c     `Hydration of chloride and bromide anions: determination 
c     of relative free energy by computer simulation'
c     J. Am. Chem. Soc. 107, 7793-7794 (1985)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.4954d0)*(2.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1070d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 79.904d0
      elelocal = 'Cl'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl-'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
