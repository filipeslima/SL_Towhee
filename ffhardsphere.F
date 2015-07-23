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
      subroutine twh_ffhardsphere(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes)
c     ******************************************************************
c     * sets up parameters for a generic Hard Sphere forcefield        *
c     *                                                                *
c     * originally written 08-11-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer ii,type,isame
      character*(FFNAMELEN) ffname
      integer nvslocal,nbslocal,ntslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'HardSphere'
c     --- set classical potential
      call twh_classical_potential(GLB_SET,'Hard Sphere')
      call twh_classical_mixrule(GLB_SET,'Arithmetic')

c     --- initialize type
      type = 0

c     --- default hard sphere (1.0 sigma)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
      masslocal = 0.60220457d0
      elelocal = 'hs'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='hs1.0'
      cbnlocal = 'hs-fix-1.0'
      canlocal = 'hs'
      ctnlocal = 'hs'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- ideal particle (0.0 sigma)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      masslocal = 0.60220457d0
      elelocal = 'hs'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='hs0.0'
      cbnlocal = 'hs-flex'
      canlocal = 'hs'
      ctnlocal = 'hs'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- half size sphere (0.5 sigma)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.5d0
      masslocal = 0.60220457d0
      elelocal = 'hs'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='hs-fix-0.5'
      cbnlocal = 'hs'
      canlocal = 'hs'
      ctnlocal = 'hs'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of 3-site water O
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.25d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O'
      cbnlocal = 'O'
      canlocal = 'O'
      ctnlocal = 'O'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of 3-site water H
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.650d0
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H'
      cbnlocal = 'H'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of clayff Cl
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 4.4d0
      masslocal = 35.453d0
      elelocal =  'Cl'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Cl'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of clayff Na
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.35d0
      masslocal = 22.990d0
      elelocal =  'Na'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'Na'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of freely jointed 3-site water O
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.25d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'FJ-O'
      cbnlocal = 'FJ-O'
      canlocal = 'FJ-O'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of 3-site water H
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.650d0
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'FJ-H'
      cbnlocal = 'FJ-H'
      canlocal = 'FJ-H'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- parameters for a hard sphere trimer
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
      masslocal = 0.60220457d0
      elelocal =  'tr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'trimer'
      cbnlocal = 'trimer'
      canlocal = 'trimer'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of rigid tangent sphere "water" 
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.25d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'R-tan-O'
      cbnlocal = 'R-tan-O'
      canlocal = 'R-tan-O'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of rigid tangent sphere "water"
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.650d0
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'R-tan-H'
      cbnlocal = 'R-tan-H'
      canlocal = 'R-tan-H'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of freely jointed tangent sphere "water"
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.25d0
      masslocal = 15.999d0
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'FJ-tan-O'
      cbnlocal = 'FJ-tan-O'
      canlocal = 'FJ-tan-O'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- a hard sphere version of freely jointed tangent sphere "water"
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.650d0
      masslocal = 1.0079d0
      elelocal =  ' H'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'FJ-tan-H'
      cbnlocal = 'FJ-tan-H'
      canlocal = 'FJ-tan-H'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * vibrations                                                     *
c     ******************************************************************

c     --- fixed bond length of 1.0
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- tangent hard spheres for regular 1.0
      isame = isame + 1
      vnlocal(isame,1) = 'hs-fix-1.0'
      vnlocal(isame,2) = 'hs-fix-1.0'
c     --- tangent hard spheres for trimer
      isame = isame + 1
      vnlocal(isame,1) = 'trimer'
      vnlocal(isame,2) = 'trimer'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- fixed bond length of 0.75
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.75d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'hs-fix-1.0'
      vnlocal(1,2) = 'hs-fix-0.5'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- fixed bond length of 0.5
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.5d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'hs-fix-0.5'
      vnlocal(1,2) = 'hs-fix-0.5'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- harmonic bond lengths for idea chain example
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0d0
      vclocal(1) = 100643.323d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'hs-flex'
      vnlocal(1,2) = 'hs-flex'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- O-H for our hard sphere version of TIP3P-like water
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.9572d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     O-H [jcp1983]
      ii = ii + 1
      vnlocal(ii,1) = 'O'
      vnlocal(ii,2) = 'H'
c     O-H for the freely jointed water
      ii = ii + 1
      vnlocal(ii,1) = 'FJ-O'
      vnlocal(ii,2) = 'FJ-H'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- fixed bond length to make the atoms tangent
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.45d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     --- rigid tangent water
      ii = ii + 1
      vnlocal(ii,1) = 'R-tan-O'
      vnlocal(ii,2) = 'R-tan-H'
c     --- freely jointed tangent water
      ii = ii + 1
      vnlocal(ii,1) = 'FJ-tan-O'
      vnlocal(ii,2) = 'FJ-tan-H'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * nonbonded bending interactions                                 *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONBOND)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'hs'
      bnlocal(1,2) = 'hs'
      bnlocal(1,3) = 'hs'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONBOND)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     --- for our freely jointed water
      ii = ii + 1
      bnlocal(ii,1) = 'FJ-H'
      bnlocal(ii,2) = 'FJ-O'
      bnlocal(ii,3) = 'FJ-H'
c     --- freely jointed tangent water
      ii = ii + 1
      bnlocal(ii,1) = 'FJ-tan-H'
      bnlocal(ii,2) = 'FJ-tan-O'
      bnlocal(ii,3) = 'FJ-tan-H'
c     --- set nbsame
      nbslocal = 1
      call twh_set_bendnames(type,nbslocal,bnlocal)


c     ******************************************************************
c     * rigid angles                                                   *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 104.52d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     H-O-H [jcp1983]
      ii = ii + 1
      bnlocal(ii,1) = 'H'
      bnlocal(ii,2) = 'O'
      bnlocal(ii,3) = 'H'
c     rigid tangent "water"
      ii = ii + 1
      bnlocal(ii,1) = 'R-tan-H'
      bnlocal(ii,2) = 'R-tan-O'
      bnlocal(ii,3) = 'R-tan-H'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 60.0d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     rigid trimer
      ii = ii + 1
      bnlocal(ii,1) = 'trimer'
      bnlocal(ii,2) = 'trimer'
      bnlocal(ii,3) = 'trimer'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * nonbonded dihedral interactions                                *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
      tnlocal(1,1) = 'hs'
      tnlocal(1,2) = 'hs'
      tnlocal(1,3) = 'hs'
      tnlocal(1,4) = 'hs'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
