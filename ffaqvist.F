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
      subroutine twh_ffaqvist(lfailure,ntypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the aqvist       *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
c     * last modified 08-04-2011 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes
c     --- local variables
      integer ii,type
      double precision adum, bdum
      character*(FFNAMELEN) ffname
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'Aqvist'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     (Ba+2) Barium (ionic +2) 
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 1341.5d0
      bdum = 24.13
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))
     &     *CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 137.33d0
      elelocal = 'Ba'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ba+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Ca+2) Calcium (ionic +2)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 264.1d0
      bdum = 18.82d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))
     &     *CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 40.08d0
      elelocal = 'Ca'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ca+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Cs+) Cesium (ionic +1)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 1647.9d0
      bdum = 5.44d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))
     &     *CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 132.90d0
      elelocal = 'Cs'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cs+'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (K+) Potasium (ionic +1)
c     --- J. Aqvist `Ion-water interaction potentials derived from
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 522.7d0
      bdum = 4.35d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))
     &     *CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 39.098d0
      elelocal = ' K'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='K+'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Li+) Lithium (ionic +1)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 25.0d0
      bdum = 2.60d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))
     &     *CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 6.941d0
      elelocal = 'Li'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Li+'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Mg+2) Magnesium (ionic +2)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 37.0d0
      bdum = 8.32d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))*CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 24.305d0
      elelocal = 'Mg'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Mg+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Na+) Sodium (ionic +1)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 143.7d0
      bdum = 3.89d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))*CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 22.990d0
      elelocal = 'Na'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na+'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Rb+) Rubidium (ionic +1)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 824.4d0
      bdum = 4.64d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))*CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 85.47d0
      elelocal = 'Rb'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Rb+'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (Sr+2) Strontium (ionic +2)
c     --- J. Aqvist `Ion-water interaction potentials derived from 
c     free energy perturbation simulations' J. Phys. Chem. 94, 
c     8021-8024 (1990)
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = 613.5d0
      bdum = 20.54d0
      nbclocal(1) = (adum/bdum)**(1.0/3.0)
      nbclocal(2) =(0.25d0*bdum*bdum/(nbclocal(1)**6))*CON_KCALTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 87.62d0
      elelocal = 'Sr'
      bplocal = 'ion'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sr+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

      return
      end
