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
      subroutine twh_ffelli2002(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up generic force field parameters for the multi-well      *
c     * potential                                                      *
c     * has now been modified so adding a new nonbonded type is enough *
c     * to set up all of the other interactions                        *
c     * bonds are set to the arithmetic average of the hard sphere of  *
c     *       the nonbonded types                                      *
c     * angles and torsions are strictly nonbonded interactions        *
c     *                                                                *
c     * originally written 04-11-2005 by M.G. Martin                   *
c     * last modified 07-21-2014 M.G. Martin                           *
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
      integer nvslocal,nbslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      integer isame,type,ii,ntslocal
      double precision sigma,center,bondlengtha,bondlengthb
      double precision lawcosa,lawcosb,angcos
      double precision bclocal
      dimension bclocal(MINBENDCOEFF:MAXBENDCOEFF)
      double precision tclocal
      dimension tclocal(MINTORCOEFF:MAXTORCOEFF)
      double precision masslocal
      character*2 elelocal 
      character*5 bplocal
      character*(FFNAMELEN) nbnlocal,canlocal,cbnlocal,ctnlocal
      double precision onepi

      ffname = 'Elli2002'
c     --- multi-well potential
      call twh_classical_potential(GLB_SET,'Multiwell')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')

c     --- initialize types
      type = 0

c     --- methane united-atom
      type = type + 1
      masslocal = 12.0107d0 + 4.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'methane'
      cbnlocal = 'null'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.674d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-138.6d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-114.1d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-42.8d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-8.5d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)

c     --- ethane united-atom
      type = type + 1
      masslocal = 12.0107d0 + 3.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'ethane'
      cbnlocal = 'alkane'
      canlocal = 'null'
      ctnlocal = 'null'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.630d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-92.5d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-62.4d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-37.2d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-24.4d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- n-butane united-atom (CH3 mass)
      type = type + 1
      masslocal = 12.0107d0 + 3.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'butane-CH3'
      cbnlocal = 'alkane'
      canlocal = 'alkane'
      ctnlocal = 'alkane'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.570d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-63.4d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-37.8d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-28.8d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-28.7d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- n-butane united-atom (CH2 mass)
      type = type + 1
      masslocal = 12.0107d0 + 2.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'butane-CH2'
      cbnlocal = 'alkane'
      canlocal = 'alkane'
      ctnlocal = 'alkane'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.570d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-63.4d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-37.8d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-28.8d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-28.7d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- n-hexane united-atom (CH3 mass)
      type = type + 1
      masslocal = 12.0107d0 + 3.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'hexane-CH3'
      cbnlocal = 'alkane'
      canlocal = 'alkane'
      ctnlocal = 'alkane'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.570d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-41.6d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-39.0d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-27.7d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-25.1d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- n-hexane united-atom (CH2 mass)
      type = type + 1
      masslocal = 12.0107d0 + 2.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'hexane-CH2'
      cbnlocal = 'alkane'
      canlocal = 'alkane'
      ctnlocal = 'alkane'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.570d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-41.6d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-39.0d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-27.7d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-25.1d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- n-octane united-atom (CH3 mass)
      type = type + 1
      masslocal = 12.0107d0 + 3.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'octane-CH3'
      cbnlocal = 'alkane'
      canlocal = 'alkane'
      ctnlocal = 'alkane'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.570d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-39.2d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-32.6d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-26.9d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-25.7d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- n-octane united-atom (CH2 mass)
      type = type + 1
      masslocal = 12.0107d0 + 2.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'octane-CH2'
      cbnlocal = 'alkane'
      canlocal = 'alkane'
      ctnlocal = 'alkane'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.570d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-39.2d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-32.6d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-26.9d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-25.7d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
c     --- benzene united-atom
      type = type + 1
      masslocal = 12.0107d0 + 1.0d0*(1.00794d0)
      elelocal = ' C'
      bplocal = 'null'
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'benzene'
      cbnlocal = 'benzene'
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      ii = twh_getnbtype(type,type)
      sigma = 3.425d0
      isame = 0
c     --- hard core
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,1.0d10)
c     --- well 1
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.2d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-66.3d0)
c     --- well 2
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.5d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-51.4d0)
c     --- well 3
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(1.8d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-51.2d0)
c     --- well 4
      isame = isame + 1
      call twh_table_pair(GLB_SET,1,ii,isame,(2.0d0)*sigma)
      call twh_table_pair(GLB_SET,2,ii,isame,-37.4d0)
c     --- set npair
      call twh_table_npair(GLB_SET,ii,isame)
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well bonds                                              *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_SQUARE_WELL)
      vclocal(1) = 1.36d0
      vclocal(2) = 1.72d0
      vclocal(3) = 0.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- alkanes [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      vnlocal(isame,1) = 'alkane'
      vnlocal(isame,2) = 'alkane'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_SQUARE_WELL)
      vclocal(1) = 1.29d0
      vclocal(2) = 1.65d0
      vclocal(3) = 0.0d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- benzene [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      vnlocal(isame,1) = 'benzene'
      vnlocal(isame,2) = 'benzene'
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
      bclocal(1) = 2.3356d0
      bclocal(2) = 2.7156d0
      bclocal(3) = 0.0d0
c     --- set a rough equilibrium angle for CBMC
      bondlengtha = 1.54d0
      bondlengthb = 1.54d0
      lawcosa = bondlengtha*bondlengtha + bondlengthb*bondlengthb
      lawcosb = 2.0d0*bondlengtha*bondlengthb
      center = 0.5d0*(bclocal(1)+bclocal(2))
      angcos = (center**2 - lawcosa)/(-lawcosb)
      bclocal(0) = twh_arccos(angcos)*180.0d0/onepi
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- alkane [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      bnlocal(isame,1) = 'alkane'
      bnlocal(isame,2) = 'alkane'
      bnlocal(isame,3) = 'alkane'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SQUARE_WELL)
c     --- set the walls
      bclocal(1) = 2.36d0
      bclocal(2) = 2.74d0
      bclocal(3) = 0.0d0
c     --- set a rough equilibrium angle for CBMC
      bondlengtha = 1.47d0
      bondlengthb = 1.47d0
      lawcosa = bondlengtha*bondlengtha + bondlengthb*bondlengthb
      lawcosb = 2.0d0*bondlengtha*bondlengthb
      center = 0.5d0*(bclocal(1)+bclocal(2))
      angcos = (center**2 - lawcosa)/(-lawcosb)
      bclocal(0) = twh_arccos(angcos)*180.0d0/onepi
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     --- benzene [personal communication J.R. Elliott 4-18-2005]
      isame = isame + 1
      bnlocal(isame,1) = 'benzene'
      bnlocal(isame,2) = 'benzene'
      bnlocal(isame,3) = 'benzene'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * square well dihedrals                                          *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,18)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      tclocal(1) = 3.57d0 - 0.15d0
      tclocal(2) = 3.57d0 + 0.15d0
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     * - alkane - alkane - * [personal communication J.R. Elliott 
c     4-18-2005]
      tnlocal(1,1) = 'alkane'
      tnlocal(1,2) = 'alkane'
      tnlocal(1,3) = 'alkane'
      tnlocal(1,4) = 'alkane'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,18)
      call twh_loftor(GLB_SET,type,.false.)
      tclocal(-1) = 0.0d0
      tclocal(1) = 3.14
      tclocal(2) = 3.425
      tclocal(3) = 0.0d0
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     * - benzene - benzene - * [personal communication J.R. Elliott 
c     4-18-2005]
      tnlocal(1,1) = 'benzene'
      tnlocal(1,2) = 'benzene'
      tnlocal(1,3) = 'benzene'
      tnlocal(1,4) = 'benzene'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
