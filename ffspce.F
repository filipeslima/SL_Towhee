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
      subroutine twh_ffspce(lfailure,ntypes,nbondtypes,nangletypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the SPC/E        *
c     * forcefield                                                     *
c     *                                                                *
c     * rewritten from several old routines 12-24-2002 by M.G. Martin  *
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
      double precision adum, bdum
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

      ffname = 'SPC-E'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
      type = 0

c     **************************************************************
c     * SPC-E potential for water (Lennard-Jones 12-6)             *
c     * Note: this is a widely used water force field and I have   *
c     * also implemented it for use with the combined exponential  *
c     * and Lennard Jones 12-6 down in the potentyp 4 section      *
c     **************************************************************
c     --- H.J.C. Berendsen; J.R. Grigera; T.P. Straatsma;
c     --- "The Missing Term in Effective Pair Potentials"
c     --- J. Phys. Chem. 91, 6269-6271 (1987).
      
c     (OW) oxygen in SPC-E water
c     --- native units in the paper are A = (kJ/mol)^(1/6)nm
c     --- and B = (kJ/mol)^(1/12)nm
      type = type + 1
      ii = twh_getnbtype(type,type)
      adum = (((0.37122d0)*(10.0d0))**6)
      bdum = (((0.3428d0)*(10.0d0))**12)
      nbclocal(1) = (bdum/adum)**(1.0/6.0)
      nbclocal(2) = (0.25d0*adum/(nbclocal(1)**6))*CON_KJTOK
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
      nbnlocal = 'OW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     (HW) hydrogen in SPC-E water
c     --- SPC-E hydrogen only has a charge and mass site
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
      nbnlocal = 'HW'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Fixed bond length for SPC-E water and simple chain models      *
c     ******************************************************************

c     --- H.J.C. Berendsen; J.R. Grigera; T.P. Straatsma;
c     --- "The Missing Term in Effective Pair Potentials"
c     --- J. Phys. Chem. 91, 6269-6271 (1987).
      type =  type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'OW'
      vnlocal(1,2) = 'HW'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Rigid angle for SPC-E water at the tetrahedral position        *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 109.47d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HW'
      bnlocal(1,2) = 'OW'
      bnlocal(1,3) = 'HW'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

      return
      end
