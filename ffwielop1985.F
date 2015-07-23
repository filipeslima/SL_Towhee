#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2007-2014 Marcus G. Martin                       *
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
      subroutine twh_ffwielop1985(lfailure,ntypes,nbondtypes
     &     ,nangletypes,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for Wielop1985       *
c     *                                                                *
c     * originally written 01-12-2007 by M.G. Martin                   *
c     * last modified 07-24-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_GETNBTYPE
#define FUNCTION_INVERSELAWOFCOSINE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,nbitypes
c     --- local variables
c     --- integer scalars
      integer ii,type,isame
c     --- character string scalars
      character*(FFNAMELEN) ffname
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
c     --- character string arrays
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)
c     --- double precision scalars
      double precision masslocal,aaa,bbb,ccc,cosangle,onepi
      double precision qbivlocal
c     --- double precision arrays
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)

      ffname = 'Wielop1985'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     --- reference
c     P. A. Wielopolski; E. R. Smith;
c     "Molecular dynamics studies of dielectric behaviour and
c     orientational correlations of liquid ethylene oxide (oxirane)
c     Mol. Phys. 54, 467-478 (1985)
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.6666d0
      nbclocal(2) = 73.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      masslocal = 15.999d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' O'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O:A'
      cbnlocal = nbnlocal
      canlocal = 'O'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.800d0
      nbclocal(2) = 90.58d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      masslocal = 15.999d0
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' O'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O:B'
      cbnlocal = nbnlocal
      canlocal = 'O'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.7143d0
      nbclocal(2) = 90.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      masslocal = 12.011d0 + 2.0d0*(1.00794d0)
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2:A'
      cbnlocal = nbnlocal
      canlocal = 'C'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.595d0
      nbclocal(2) = 72.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      masslocal = 12.011d0 + 2.0d0*(1.00794d0)
      call twh_mass(GLB_SET,type,masslocal) 
      elelocal = ' C'
      call twh_elementname(GLB_SET,type,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2:B'
      cbnlocal = nbnlocal
      canlocal = 'C'
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Fixed bond length
c     ******************************************************************

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.431d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
      isame = isame + 1
      vnlocal(isame,1) = 'O:A'
      vnlocal(isame,2) = 'CH2:A'
      isame = isame + 1
      vnlocal(isame,1) = 'O:B'
      vnlocal(isame,2) = 'CH2:B'
      call twh_set_vibnames(type,isame,vnlocal)

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.466d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
      isame = isame + 1
      vnlocal(isame,1) = 'CH2:A'
      vnlocal(isame,2) = 'CH2:A'
      isame = isame + 1
      vnlocal(isame,1) = 'CH2:B'
      vnlocal(isame,2) = 'CH2:B'
      call twh_set_vibnames(type,isame,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Rigid angle
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      aaa = 1.431d0
      bbb = 1.466d0
      ccc = 1.431d0
      cosangle = twh_inverselawofcosine(aaa,bbb,ccc)
      bclocal(0) = twh_arccos(cosangle)*180.0d0/onepi
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
      isame = isame + 1
      bnlocal(isame,1) = 'O'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
      call twh_set_bendnames(type,isame,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      aaa = 1.431d0
      bbb = 1.431d0
      ccc = 1.466d0
      cosangle = twh_inverselawofcosine(aaa,bbb,ccc)
      bclocal(0) = twh_arccos(cosangle)*180.0d0/onepi
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'O'
      bnlocal(isame,3) = 'C'
      call twh_set_bendnames(type,isame,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     --- bond increment

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.1608d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C-O (A model)
      qbinlocal(1) = 'CH2:A'
      qbinlocal(2) = 'O:A'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.1608d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C-O (B model)
      qbinlocal(1) = 'CH2:B'
      qbinlocal(2) = 'O:B'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C-C (A model)
      qbinlocal(1) = 'CH2:A'
      qbinlocal(2) = 'CH2:A'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.0d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C-C (B model)
      qbinlocal(1) = 'CH2:B'
      qbinlocal(2) = 'CH2:B'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
