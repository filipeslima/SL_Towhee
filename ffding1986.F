#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2011 Marcus G. Martin                       *
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
      subroutine twh_ffding1986(lfailure,ntypes)
c     ******************************************************************
c     * sets up parameters for the Ding1986 forcefield                 *
c     *                                                                *
c     * originally written 09-30-2004 by M.G. Martin                   *
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

      ffname = 'Ding1986'
      call twh_classical_potential(GLB_SET,'Stillinger-Weber')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     ******************************************************************
c     This forcefield is a Stillinger-Weber Potential published in   
c     K. Ding; H.C. Andersen; "Molecular-dynamics simulation of 
c     amorphous germanium", Phys. Rev. B 34 6987-6991 (1986). 
c     ******************************************************************

c     --- initialize type
      type = 0
c     --- germanium
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.93d0*CON_EVTOK
      nbclocal(2) = 7.049556277d0
      nbclocal(3) = 0.6022245584d0
      nbclocal(4) = 2.181d0
      nbclocal(5) = 4.0d0
      nbclocal(6) = 1.8d0
      nbclocal(7) = 31.0d0
      nbclocal(8) = 1.2d0
      nbclocal(9) = -1.0d0/3.0d0
      masslocal = 72.61d0
      elelocal = 'Ge'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ge'
      cbnlocal = 'Ge'
      canlocal = 'Ge'
      ctnlocal = 'Ge'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
