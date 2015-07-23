#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2005-2014 Marcus G. Martin                       *
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
      subroutine twh_ffgala1994(lfailure,ntypes,nbondtypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the gala1994     *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 02-08-2005 by M.G. Martin                   *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nbitypes
c     --- local variables
      integer type,ii,isame
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
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)

      ffname = 'Gala1994'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      type = 0

c     --- reference
c     G. Galassi; D.J. Tildesley;
c     "Phase Diagrams of Diatomic Molecules Using the Gibbs Ensemble
c     Monte Carlo Method", Mol. Sim. 13 (1994) 11-24.

c     Nitrogen in N2
      type = type + 1
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'N'
      cbnlocal = nbnlocal
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.310d0
      nbclocal(2) = 37.30d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     Fluorine in F2
      type = type + 1
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F'
      cbnlocal = nbnlocal
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.825d0
      nbclocal(2) = 52.8d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     Chlorine in Cl2
      type = type + 1
      masslocal = 35.453d0
      elelocal = 'Cl'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl'
      cbnlocal = nbnlocal
      canlocal = 'none'
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.332d0
      nbclocal(2) = 178.3d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond lengths - rigid                                           *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.329d0)*(3.310d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- dinitrogen
      isame = isame + 1
      vnlocal(isame,1) = 'N'
      vnlocal(isame,2) = 'N'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.505d0)*(2.825d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- difluorine
      isame = isame + 1
      vnlocal(isame,1) = 'F'
      vnlocal(isame,2) = 'F'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.630d0)*(3.332d0)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- dichlorine
      isame = isame + 1
      vnlocal(isame,1) = 'Cl'
      vnlocal(isame,2) = 'Cl'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'N'
      qbinlocal(2) = 'N'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'F'
      qbinlocal(2) = 'F'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
      qbinlocal(1) = 'Cl'
      qbinlocal(2) = 'Cl'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
