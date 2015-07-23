#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2010-2014 Marcus G. Martin                       *
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
      subroutine twh_ffvahid2010(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * a SPEADMD model implementation based upon the work published in*
c     * A. Vahid; J. R. Elliott; "Transferable Intermolecular          *
c     * Potentials for Carboxylic Acids and Their Phase Behavior"      *
c     * AIChE Journal, 56, 2, 485-505 (2010)                           *
c     *                                                                *
c     * originally written 12-03-2005 by M.G. Martin                   *
c     * last modified 07-18-2014 M.G. Martin                           *
c     ******************************************************************
      implicit none
#define FUNCTION_ARCCOS
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
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
      integer isame,type,ii
      double precision center,bondlengtha,bondlengthb
      double precision lawcosa,lawcosb,angcos
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
      double precision onepi

      ffname = 'Cui2002'
c     --- multi-well potential
      call twh_classical_potential(GLB_SET,'Hard 2580 Multistep')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     --- all values implemented here are directly from the paper
c     J. Cui, J. R. Elliott Jr., "Phase diagrams for a multistep
c     potential model of n-alkanes by discontinuous molecular dynamics
c     and thermodynamic perturbation theory
c     J. Chem Phys., 116, 8625-8631 (2002)

c     --- set the nonbonded parameters for this force field

c     --- initialize type
      type = 0

      type = type + 1
      masslocal =   3.0d0*1.00794d0 + 12.0107d0
      elelocal =  ' C'
      bplocal = 's'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH3'
      cbnlocal = 'C'
      canlocal = 'C'
      ctnlocal = 'C'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   (0.363d0)*(10.0d0)
      nbclocal(2) =  -(84.0d0)
      nbclocal(3) =  -(73.9d0)
      nbclocal(4) =  -(35.9d0)
      nbclocal(5) =  -(17.5d0)
      nbclocal(6) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      type = type + 1
      masslocal =   2.0d0*1.00794d0 + 12.0107d0
      elelocal =  ' C'
      bplocal = 'sp3'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'CH2'
      cbnlocal = 'C'
      canlocal = 'C'
      ctnlocal = 'C'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      nbclocal(1) =   (0.357d0)*(10.0d0)
      nbclocal(2) =  -(30.1d0)
      nbclocal(3) =  -(25.8d0)
      nbclocal(4) =  -(25.8d0)
      nbclocal(5) =  -(22.9d0)
      nbclocal(6) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well bonds                                              *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_SQUARE_WELL)
      vclocal(1) = (0.137d0)*(10.0d0)
      vclocal(2) = (0.173d0)*(10.0d0)
      vclocal(3) = 0.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- C - C
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well angles                                             *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)
       
      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SQUARE_WELL)
c     --- set the walls
      bclocal(1) = (0.236d0)*(10.0d0)
      bclocal(2) = (0.272d0)*(10.0d0)
      bclocal(3) = 0.0d0
c     --- set a rough equilibrium angle for CBMC
      bondlengtha = 1.55d0
      bondlengthb = 1.55d0
      lawcosa = bondlengtha*bondlengtha + bondlengthb*bondlengthb
      lawcosb = 2.0d0*bondlengtha*bondlengthb
      center = 0.5d0*(bclocal(1) + bclocal(2))
      angcos = (center**2 - lawcosa)/(-lawcosb)
      bclocal(0) = twh_arccos(angcos)*180.0d0/onepi
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- C-C-C
      isame = isame + 1
      bnlocal(isame,1) = 'C'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * nonbonded interactions for dihedrals                           *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- C-C-C-C
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
