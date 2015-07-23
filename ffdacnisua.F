#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2002-2014 Marcus G. Martin                       *
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
      subroutine twh_ffdacnisua(lfailure,ntypes,nbondtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the dacnis-ua    *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-20-2002 by M.G. Martin  *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      integer nvslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'DACNIS-UA'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      type = 0

c     **************************************************************
c     * DACNIS-UA force field                                      *
c     **************************************************************
c     force field developed for alkane absorption in silicalite
c     not clear if any more United atom parameters will be developed
c     for this force field, but here it is if you want to look at 
c     methane and ethane in silicalite
c     References:
c     [mart2001] M.G. Martin, A.P. Thompson, T.M. Nenoff
c     `Effect of pressure, membrane thickness, and 
c     placement of control volumes on the flux of methane
c     through thin silicalite membranes: A dual control
c     volume grand canonical molecular dynamics study'
c     J. Chem. Phys. 114, 7174-7181 (2001)
      
c     (CH4sp3) Methane sp3 (final 12-08-1999 MGM)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.72d0
      nbclocal(2) = 158.5d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 16.043d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH4sp3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      
c     (CH3sp3) alkane methyl group sp3 (final 1-17-2000 MGM)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.76d0
      nbclocal(2) = 108.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.0347d0
      elelocal = ' C'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CH3sp3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
          
c     (Ozeo) zeolite Oxygen [includes 0.5 Si mass] 
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.20d0
      nbclocal(2) = 82.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 30.042d0
      elelocal = ' O'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ozeo'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      bplocal = 'null'
      call twh_bondpatt(GLB_SET,type,bplocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     [CHx](sp3) -- [CHy](sp3,sp2,aro) single bonds
c     --- from table 1.6 of McMurry Organic Chemistry 3rd edition
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.54d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CH3sp3'
      vnlocal(1,2) = 'CH3sp3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

      return
      end
