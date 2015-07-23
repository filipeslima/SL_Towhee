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
      subroutine twh_ffljium(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up parameters for a generic Lennard-Jones forcefield      *
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
      integer ii,type
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
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision nbclocal
      dimension nbclocal(MINNBCOEFF:MAXNBCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal

      ffname = 'LJium'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     --- default lennard-jonesium (1,1)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
      nbclocal(2) = 1.0d0
      nbclocal(3) = 1.0d0
      nbclocal(4) = 1.0d0
      masslocal = 0.60220457d0
      elelocal = 'lj'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lj1.0'
      cbnlocal = 'lj'
      canlocal = 'lj'
      ctnlocal = 'lj'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- FENE lennard-jonesium (1,1)
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
      nbclocal(2) = 1.0d0
      nbclocal(3) = 1.0d0
      nbclocal(4) = 1.0d0
      masslocal = 0.60220457d0
      elelocal = 'lj'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lj1.0_FENE'
      cbnlocal = 'ljFENE'
      canlocal = 'lj'
      ctnlocal = 'lj'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * vibrations                                                     *
c     ******************************************************************

c     --- fixed bond length of 1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'lj'
      vnlocal(1,2) = 'lj'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- FENE [Kremer & Grest J. Chem. Phys. 92, 5057-5086 (1990)]
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FENE)
      vclocal(1) = 1.5d0
      vclocal(2) = 30.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'ljFENE'
      vnlocal(1,2) = 'ljFENE'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * nonbonded bending interactions                                 *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONBOND)
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'lj'
      bnlocal(1,2) = 'lj'
      bnlocal(1,3) = 'lj'
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
      tnlocal(1,1) = 'lj'
      tnlocal(1,2) = 'lj'
      tnlocal(1,3) = 'lj'
      tnlocal(1,4) = 'lj'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
