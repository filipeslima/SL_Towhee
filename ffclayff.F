#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2004-2014 Marcus G. Martin                       *
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
      subroutine twh_ffclayff(lfailure,ntypes,nbondtypes,nangletypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the clayff       *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 06-16-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes
c     --- local variables
      integer ii,type
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'ClayFF'
      type = 0
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
c     ******************************************************************
c     * Clayff references                                              *
c     ******************************************************************
c     Randall T. Cygan; Jian-Jie Liang; Andrey G. Kalinichev;
c     `Molecular Models of Hydroxide, Oxyhydroxide, and Clay Phases and
c     the Development of a General Force Field'  J. Phys. Chem. B
c     2004, 108, 1255-1266.

c     -- Hydrogen
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     water hydrogen (standard charge 0.4100)
      nbnlocal = 'h*'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (0.0d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     hydroxyl hydrogen (standard charge 0.4250)
      nbnlocal = 'ho'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     -- Oxygen --

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     water oxygen (standard charge -0.8200)
      nbnlocal = 'o*'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     hydroxyl oxygen (standard charge -0.9500)
      nbnlocal = 'oh'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     bridging oxygen (standard charge -1.0500)
      nbnlocal = 'ob'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     bridging oxygen with octahedral substitution (charge -1.1808)
      nbnlocal = 'obos'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     bridging oxygen with tetrahedral substitution (charge -1.1688)
      nbnlocal = 'obts'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     bridging oxygen with double substitution (charge -1.2996)
      nbnlocal = 'obss'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.5532d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1554d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     hydroxyl oxygen with substitution (charge -1.0808)
      nbnlocal = 'ohs'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     -- Silicon --

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.7064d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (1.8405d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 28.086d0
      elelocal = 'Si'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     tetrahedral silicon (charge 2.1000)
      nbnlocal = 'st'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     -- Aluminum --

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.7943d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (1.3298d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 26.982d0
      elelocal = 'Al'
      bplocal = 'p3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     octahedral aluminum (charge 1.5750)
      nbnlocal = 'ao'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.7064d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (1.8405d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 26.982d0
      elelocal = 'Al'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     tetrahedral aluminum (charge 1.5750)
      nbnlocal = 'at'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- Magnesium

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (5.9090d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (9.0298d-7)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 24.305d0
      elelocal = 'Mg'
      bplocal = 'p3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     octahedral magnesium (charge 1.3600)
      nbnlocal = 'mgo'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (5.9090d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (9.0298d-7)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 24.305d0
      elelocal = 'Mg'
      bplocal = 's?'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     hydroxide magnesium (charge 1.0500)
      nbnlocal = 'mgh'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (6.2484d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (5.0298d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 40.08d0
      elelocal = 'Ca'
      bplocal = 'p3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     octahedral calcium (charge 1.3600)
      nbnlocal = 'cao'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (6.2484d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (5.0298d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 40.08d0
      elelocal = 'Ca'
      bplocal = 's?'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     hydroxide calcium (charge 1.0500)
      nbnlocal = 'cah'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (5.5070d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (9.0298d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 55.847d0
      elelocal = 'Fe'
      bplocal = 'p3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     octahedral iron (charge 1.5750)
      nbnlocal = 'feo'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.7257d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (9.0298d-6)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 6.941d0
      elelocal = 'Li'
      bplocal = 'p3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     octahedral lithium (charge 0.5250)
      nbnlocal = 'lio'
      cbnlocal = 'metal'
      canlocal = 'metal'
      ctnlocal = 'generic'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.6378d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1301d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 22.990d0
      elelocal = 'Na'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     aqueous sodium ion (charge 1.0)
      nbnlocal = 'Na'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.7423d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 39.098d0
      elelocal = ' K'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     aqueous potassium ion (charge 1.0)
      nbnlocal = 'K'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.3002d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 132.90d0
      elelocal = 'Cs'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     aqueous cesium ion (charge 1.0)
      nbnlocal = 'Cs'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.2237d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1000d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 40.08d0
      elelocal = 'Ca'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     aqueous calcium ion (charge 2.0)
      nbnlocal = 'Ca'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.2840d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.0470d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 137.33d0
      elelocal = 'Ba'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     aqueous barium ion (charge 2.0)
      nbnlocal = 'Ba'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.9388d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.1001d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 35.453d0
      elelocal = 'Cl'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     aqueous chloride ion (charge -1.0)
      nbnlocal = 'Cl'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * ClayFF bond potentials                                         *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0000d0
      vclocal(1) = (554.1349d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     o* - h*
      ii = ii + 1
      vnlocal(ii,1) = 'o*'
      vnlocal(ii,2) = 'h*'
c     oh - ho
      ii = ii + 1
      vnlocal(ii,1) = 'oh'
      vnlocal(ii,2) = 'ho'
c     ohs - ho
      ii = ii + 1
      vnlocal(ii,1) = 'ohs'
      vnlocal(ii,2) = 'ho'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- everything else just interacts via nonbond
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_NONBOND)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     metal - oh
      ii = ii + 1
      vnlocal(ii,1) = 'metal'
      vnlocal(ii,2) = 'oh'
c     metal - ob
      ii = ii + 1
      vnlocal(ii,1) = 'metal'
      vnlocal(ii,2) = 'ob'
c     metal - obts
      ii = ii + 1
      vnlocal(ii,1) = 'metal'
      vnlocal(ii,2) = 'obts'
c     metal - obss
      ii = ii + 1
      vnlocal(ii,1) = 'metal'
      vnlocal(ii,2) = 'obss'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)
     

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * ClayFF harmonic angles                                         *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (45.7696d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
c     h* - o* - h*
      bnlocal(1,1) = 'h*'
      bnlocal(1,2) = 'o*'
      bnlocal(1,3) = 'h*'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (30.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     Metal - oh - ho 
      ii = ii + 1
      bnlocal(ii,1) = 'metal'
      bnlocal(ii,2) = 'oh'
      bnlocal(ii,3) = 'ho'
c     Metal - ohs - ho 
      ii = ii + 1
      bnlocal(ii,1) = 'metal'
      bnlocal(ii,2) = 'ohs'
      bnlocal(ii,3) = 'ho'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)


      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end
