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
      subroutine twh_fftele1987(lfailure,ntypes,nbondtypes,nangletypes
     &     ,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the tele1987     *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 02-09-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,nbitypes
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
      double precision qbivlocal 
      character*(FFNAMELEN) qbinlocal 
      dimension qbinlocal(2)

      ffname = 'Tele1987'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     reference
c     O. Teleman; B. Jonsson; S. Engstrom;
c     "A molecular dynamics simulation of a water model with 
c     intramolecular degrees of freedom", Mol. Phys. 60 193-203 (1987)

c     (O) oxygen
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.166d0
      nbclocal(2) = (0.65017)*CON_KJTOK
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (H) hydrogen
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * harmonic bond length                                           *
c     ******************************************************************

      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0d0
      vclocal(1) = (4637.0d0)*(0.5d0)*CON_KJTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'O'
      vnlocal(1,2) = 'H'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * harmonic angles                                                *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.47d0
      bclocal(1) = (383.0d0)*(0.5d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'H'
      bnlocal(1,2) = 'O'
      bnlocal(1,3) = 'H'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = 0.41d0
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - O in water
      qbinlocal(1) = 'H'
      qbinlocal(2) = 'O'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
