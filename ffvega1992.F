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
      subroutine twh_ffvega1992(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up a square well potential for the vega1992 force field   *
c     * this force field was used in a study of square well monomers   *
c     * ranging in lamda values from 1.25 to 2.0                       *
c     * I added some default values for the bonded terms for those who *
c     * might want to study square well chains                         *
c     *                                                                *
c     * originally written 06-24-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer ii,jj,itype,jtype,type
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

      ffname = 'Vega1992'
c     --- square well potential
      call twh_classical_potential(GLB_SET,'Square Well')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     --- initialize types
      type = 0

c     --- type 1.25
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
c     --- set the mass so that the g/ml becomes a reduced number density
      masslocal = 0.60220457d0/(nbclocal(1)**3)
      call twh_mass(GLB_SET,type,masslocal) 
      nbclocal(2) = 1.25d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lam_1.25'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- type 1.375
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
c     --- set the mass so that the g/ml becomes a reduced number density
      masslocal = 0.60220457d0/(nbclocal(1)**3)
      call twh_mass(GLB_SET,type,masslocal) 
      nbclocal(2) = 1.375d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lam_1.375'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- type 1.5
      type = type + 1
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lam_1.5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
c     --- set the mass so that the g/ml becomes a reduced number density
      masslocal = 0.60220457d0/(nbclocal(1)**3)
      nbclocal(2) = 1.5d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- type 1.75
      type = type + 1
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lam_1.75'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
c     --- set the mass so that the g/ml becomes a reduced number density
      masslocal = 0.60220457d0/(nbclocal(1)**3)
      nbclocal(2) = 1.75d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     --- type 1.2
      type = type + 1
      elelocal = 'Sw'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='lam_1.2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 1.0d0
c     --- set the mass so that the g/ml becomes a reduced number density
      masslocal = 0.60220457d0/(nbclocal(1)**3)
      nbclocal(2) = 1.2d0
      nbclocal(3) = 1.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     --- bond potentials
c     --- fixed bond length that is the average of the inner hard 
c     --- sphere radii

      do itype = 1,ntypes
         ii = twh_getnbtype(itype,itype)
         do jtype = itype,ntypes
            jj = twh_getnbtype(jtype,jtype)
c           --- increment bond type
            type = type + 1
            call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
c           --- set to the arithmetic average of the hard radii
            vclocal(0) = 0.5d0*(twh_get_nbcoeff(1,ii)
     &           + twh_get_nbcoeff(1,jj))
            call twh_set_vibcoeffs(type,vclocal)
            call twh_vibff(GLB_SET,type,ffname)
            nvslocal = 1
c           --- use the bond names of those types
            vnlocal(1,1) = twh_get_cbname(itype)
            vnlocal(1,2) = twh_get_cbname(jtype)
            call twh_set_vibnames(type,nvslocal,vnlocal)
         enddo
      enddo

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return
       
c     --- angle potentials nonbond interaction for all
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_NONBOND)
      call twh_bendff(GLB_SET,type,ffname)
c     wild - wild - wild
      nbslocal = 1
      bnlocal(1,1) = 'wild'
      bnlocal(1,2) = 'wild'
      bnlocal(1,3) = 'wild'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     --- torsion potentials nonbond interaction for all
      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- wild - wild - wild - wild
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'wild'
      tnlocal(1,3) = 'wild'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
