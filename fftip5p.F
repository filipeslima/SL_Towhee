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
      subroutine twh_fftip5p(lfailure,ntypes,nbondtypes,nangletypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the TIP5P        *
c     * forcefield for water                                           *
c     *                                                                *
c     * originally written 08-27-2004 by M.G. Martin                   *
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

c     --- set variables for this forcefield
      ffname = 'TIP5P'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     ******************************************************************
c     * TIP5P force field for water                                    *
c     ******************************************************************
c     [jcp2000] M.W. Mahoney; W.L. Jorgensen; 'A five-site model for
c     liquid water and the reprodution of the density anomaly by rigid,
c     nonpolarizable potential functions', J. Chem. Phys. 112, 
c     8910-8922 (2000).


c     ** Hydrogen **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = (0.5d0)*nbclocal(2)
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     H in tip5p water [jcp2000]
      nbnlocal = 'H'
      cbnlocal = 'H'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     ** Oxygen **

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.12d0
      nbclocal(2) = (0.16d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     water oxygen in the TIP5P model [jcp2000]
      nbnlocal = 'O'
      cbnlocal = 'O'
      canlocal = 'O'
      ctnlocal = 'O'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     ** other non-atoms **
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 0.0
      nbclocal(2) = (0.0d0)*(CON_KCALTOK)
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 0.0d0
      elelocal = '  '
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     L negative charge sites on the TIP5P water [jcp2000]
      nbnlocal = 'L'
      cbnlocal = 'L'
      canlocal = 'L'
      ctnlocal = 'L'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed bond lengths                                             *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.9572d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     O-H [jcp2000]
      ii = ii + 1
      vnlocal(ii,1) = 'O'
      vnlocal(ii,2) = 'H'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.70d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     O-L distance between oxygen and the lone pairs [jcp2000]
      ii = ii + 1
      vnlocal(ii,1) = 'O'
      vnlocal(ii,2) = 'L'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed angles                                                   *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 104.52d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     H-O-H [jcp2000]
      ii = ii + 1
      bnlocal(ii,1) = 'H'
      bnlocal(ii,2) = 'O'
      bnlocal(ii,3) = 'H'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 109.47d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     L-O-L [jcp2000]
      ii = ii + 1
      bnlocal(ii,1) = 'L'
      bnlocal(ii,2) = 'O'
      bnlocal(ii,3) = 'L'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 110.69478d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
C     --- L - O - H by geometry [jcp2000]
      ii = ii + 1
      bnlocal(ii,1) = 'L'
      bnlocal(ii,2) = 'O'
      bnlocal(ii,3) = 'H'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     --- safety check angles
      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end
