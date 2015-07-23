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
      subroutine twh_ffdreiding(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the dreiding     *
c     * forcefield                                                     *
c     * Note: Dreiding force field transforms its basic                *
c     * parameters into use with a variety of functional forms.        *
c     *                                                                *
c     * originally written 12-08-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_CANAME
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      logical lgroupi,lgroupj,lfound
      character*1 lettera,letterb,letteri,letterj,letterk
      character*(FFNAMELEN) code,ffname
      integer nvslocal,nbslocal,ntslocal
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      integer nallatom,ii,jj,itype,jtype,type,isame,ktype,index,icoeff
      double precision sixroottwo,bondorder,theta
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
      double precision iclocal
      dimension iclocal(MINIMPCOEFF:MAXIMPCOEFF)
      integer nislocal 
      character*(FFNAMELEN) ninlocal 
      dimension ninlocal(MAXSAMEIMP,MAXIMPNAME)
      double precision onepi

      ffname = 'DREIDING'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'LB or Geometric')
c     --- initialize types
      type = 0

c     ******************************************************************
c     * DREIDING                                                       *
c     ******************************************************************
c     References:
c     [jpc1990]
c     S.L. Mayo; B.D. Olafson; W.A. Goddard III;
c     "DREIDING: A Generic Force Field for Molecular Simulations"
c     J. Phys. Chem. 94, 8897-8909 (1990)

c     --- setup all of the basic parameters for later conversion
c     --- into the actual forcefield parameters
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.330d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.195d0
      nbclocal(8) = 0.0152d0
      nbclocal(9) = 12.382d0
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H_'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.510d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.195d0
      nbclocal(8) = 0.0001d0
      nbclocal(9) = 12.0d0
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='H_b'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.880d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.02d0
      nbclocal(8) = 0.095d0
      nbclocal(9) = 14.23d0
      masslocal = 10.811d0
      elelocal = ' B'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='B_3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.790d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 4.02d0
      nbclocal(8) = 0.095d0
      nbclocal(9) = 14.23d0
      masslocal = 10.811d0
      elelocal = ' B'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='B_2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.770d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 3.8983d0
      nbclocal(8) = 0.0951d0
      nbclocal(9) = 14.034d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.700d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.8983d0
      nbclocal(8) = 0.0951d0
      nbclocal(9) = 14.034d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.670d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.8983d0
      nbclocal(8) = 0.0951d0
      nbclocal(9) = 14.034d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.602d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.8983d0
      nbclocal(8) = 0.0951d0
      nbclocal(9) = 14.034d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.702d0
      nbclocal(6) = 106.7d0
      nbclocal(7) = 3.6621d0
      nbclocal(8) = 0.0774d0
      nbclocal(9) = 13.843d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N_3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.650d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.6621d0
      nbclocal(8) = 0.0774d0
      nbclocal(9) = 13.843d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N_R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.615d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.6621d0
      nbclocal(8) = 0.0774d0
      nbclocal(9) = 13.843d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N_2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.556d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.6621d0
      nbclocal(8) = 0.0774d0
      nbclocal(9) = 13.843d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'sp'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='N_1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.660d0
      nbclocal(6) = 104.51d0
      nbclocal(7) = 3.4046d0
      nbclocal(8) = 0.0957d0
      nbclocal(9) = 13.483d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O_3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.660d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.4046d0
      nbclocal(8) = 0.0957d0
      nbclocal(9) = 13.483d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O_R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.560d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.4046d0
      nbclocal(8) = 0.0957d0
      nbclocal(9) = 13.483d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O_2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.528d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.4046d0
      nbclocal(8) = 0.0957d0
      nbclocal(9) = 13.483d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O_1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.611d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.4720d0
      nbclocal(8) = 0.0725d0
      nbclocal(9) = 14.444d0
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F_'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.047d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.39d0
      nbclocal(8) = 0.31d0
      nbclocal(9) = 12.0d0
      masslocal = 26.982d0
      elelocal = 'Al'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Al3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.937d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.27d0
      nbclocal(8) = 0.31d0
      nbclocal(9) = 12.0d0
      masslocal = 28.086d0
      elelocal = 'Si'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Si3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.890d0
      nbclocal(6) = 93.3d0
      nbclocal(7) = 4.1500d0
      nbclocal(8) = 0.3200d0
      nbclocal(9) = 12.0d0
      masslocal = 30.974d0
      elelocal = ' P'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P_3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.040d0
      nbclocal(6) = 92.1d0
      nbclocal(7) = 4.0300d0
      nbclocal(8) = 0.3440d0
      nbclocal(9) = 12.0d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S_3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.997d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.9503d0
      nbclocal(8) = 0.2833d0
      nbclocal(9) = 13.861d0
      masslocal = 35.453d0
      elelocal = 'Cl'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cl'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.210d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.39d0
      nbclocal(8) = 0.40d0
      nbclocal(9) = 12.0d0
      masslocal = 69.723d0
      elelocal = 'Ga'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ga3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.210d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.27d0
      nbclocal(8) = 0.40d0
      nbclocal(9) = 12.0d0
      masslocal = 72.61d0
      elelocal = 'Ge'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ge3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.210d0
      nbclocal(6) = 92.1d0
      nbclocal(7) = 4.15d0
      nbclocal(8) = 0.41d0
      nbclocal(9) = 12.0d0
      masslocal = 74.922d0
      elelocal = 'As'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='As3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.210d0
      nbclocal(6) = 90.6d0
      nbclocal(7) = 4.03d0
      nbclocal(8) = 0.43d0
      nbclocal(9) = 12.0d0
      masslocal = 78.96d0
      elelocal = 'Se'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Se3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.167d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.95d0
      nbclocal(8) = 0.37d0
      nbclocal(9) = 12.0d0
      masslocal = 79.904d0
      elelocal = 'Br'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Br'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.390d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.59d0
      nbclocal(8) = 0.55d0
      nbclocal(9) = 12.0d0
      masslocal = 114.82d0
      elelocal = 'In'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='In3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.373d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.47d0
      nbclocal(8) = 0.55d0
      nbclocal(9) = 12.0d0
      masslocal = 118.71d0
      elelocal = 'Sn'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sn3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.432d0
      nbclocal(6) = 91.6d0
      nbclocal(7) = 4.35d0
      nbclocal(8) = 0.55d0
      nbclocal(9) = 12.0d0
      masslocal = 121.75d0
      elelocal = 'Sb'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sb3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.280d0
      nbclocal(6) = 90.3d0
      nbclocal(7) = 4.23d0
      nbclocal(8) = 0.57d0
      nbclocal(9) = 12.0d0
      masslocal = 127.60d0
      elelocal = 'Te'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Te3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.360d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.15d0
      nbclocal(8) = 0.51d0
      nbclocal(9) = 12.0d0
      masslocal = 126.90d0
      elelocal = ' I'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='I_'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.860d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.144d0
      nbclocal(8) = 0.5d0
      nbclocal(9) = 12.0d0
      masslocal = 22.990d0
      elelocal = 'Na'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na_+'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.940d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.472d0
      nbclocal(8) = 0.05d0
      nbclocal(9) = 12.0d0
      masslocal = 40.08d0
      elelocal = 'Ca'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ca_+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.285d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 4.54d0
      nbclocal(8) = 0.055d0
      nbclocal(9) = 12.0d0
      masslocal = 55.847d0
      elelocal = 'Fe'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fe_+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.330d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.54d0
      nbclocal(8) = 0.055d0
      nbclocal(9) = 12.0d0
      masslocal = 65.39d0
      elelocal = 'Zn'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Zn_+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     --- set the number of all-atom types as these are the only
c     --- ones that need intramolecular parameters set since the UA
c     --- use the same intramolecular parameters as the parent type
      nallatom = type

c     --- from here on down are all the united-atoms

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.770d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.2370d0
      nbclocal(8) = 0.3016d0
      nbclocal(9) = 12.0d0
      masslocal = 12.011d0 + 4.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C_34'
      cbnlocal = 'C_3'
      canlocal = 'C_3'
      ctnlocal = 'C_3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.770d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.1524d0
      nbclocal(8) = 0.2500d0
      nbclocal(9) = 12.0d0
      masslocal = 12.011d0 + 3.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_33'
      cbnlocal = 'C_3'
      canlocal = 'C_3'
      ctnlocal = 'C_3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.770d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 4.0677d0
      nbclocal(8) = 0.1984d0
      nbclocal(9) = 12.0d0
      masslocal = 12.011d0 + 2.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_32'
      cbnlocal = 'C_3'
      canlocal = 'C_3'
      ctnlocal = 'C_3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.770d0
      nbclocal(6) = 109.471d0
      nbclocal(7) = 3.9830d0
      nbclocal(8) = 0.1467d0
      nbclocal(9) = 12.0d0
      masslocal = 12.011d0 + 1.0d0*1.0079d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='C_31'
      cbnlocal = 'C_3'
      canlocal = 'C_3'
      ctnlocal = 'C_3'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.700d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 4.23d0
      nbclocal(8) = 0.1356d0
      nbclocal(9) = 14.034d0
      masslocal = 12.011d0 + 1.0079d0
      elelocal = ' C'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C_R1'
      cbnlocal = 'C_R'
      canlocal = 'C_R'
      ctnlocal = 'C_R'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return
c     ******************************************************************
c     * transform into nonbond                                         *
c     ******************************************************************
c     here we follow the suggestions in the DREIDING paper to use
c     (a) zeta = 12.0
c     (c) use the same R0 and D0 for 12-6 as exp-6

      sixroottwo = 2.0d0**(1.0d0/6.0d0)
      do itype = 1,ntypes
         ii = twh_getnbtype(itype,itype)
c        --- convert from R_0 units to sigma
         nbclocal(1) = twh_get_nbcoeff(7,ii)/sixroottwo
c        --- convert Di from kcal to kelvin
         nbclocal(2) = twh_get_nbcoeff(8,ii)*CON_KCALTOK
c        --- UFF includes full VDW and coulombic interactions for 1-4
         nbclocal(3) = nbclocal(1)
         nbclocal(4) = nbclocal(2)
c        --- set manually so we do not overwrite the other parameters
         do icoeff = 1,4
            call twh_nbcoeff(GLB_SET,icoeff,ii,nbclocal(icoeff))
         enddo
      enddo

c     ******************************************************************
c     * bonds                                                          *
c     ******************************************************************

c     --- initialize type
      type = 0
      do itype = 1,nallatom
         ii = twh_getnbtype(itype,itype)
         code = twh_get_cbname(itype)
         lettera = code(3:3)
         do jtype = itype,nallatom
            jj = twh_getnbtype(jtype,jtype)
            code = twh_get_cbname(jtype)
            letterb = code(3:3)
            type = type + 1
            call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
c           --- compute the equilibrium bond length
            vclocal(0) = twh_get_nbcoeff(5,ii) + twh_get_nbcoeff(5,jj)
     &           - 0.01d0
            if ( lettera .eq. '1' .and. letterb .eq. '1' ) then
c               --- triple bond
               bondorder = 3.0d0
            elseif ( lettera .eq. '2' .and. letterb .eq. '2' ) then
c              --- double bond
               bondorder = 2.0d0
            elseif ( lettera .eq. 'R' .and. letterb .eq. 'R' ) then
c               --- resonant
               bondorder = 1.50d0
            else
c               --- every thing else is a single bond
               bondorder = 1.0d0
            endif
            vclocal(1) = bondorder*CON_KCALTOK*(0.5d0)*(700.0d0)
c           --- set the rest of the information for this vibration
            call twh_set_vibcoeffs(type,vclocal)
            call twh_vibff(GLB_SET,type,ffname)
            nvslocal = 1
            vnlocal(1,1) = twh_get_cbname(itype)
            vnlocal(1,2) = twh_get_cbname(jtype)
            call twh_set_vibnames(type,nvslocal,vnlocal)
         enddo
      enddo

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * angles                                                         *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

      do itype = 1,nallatom
         ii = twh_getnbtype(itype,itype)
         theta = twh_get_nbcoeff(6,ii)
         type = type + 1
c        --- set equilibrium angle (degrees)
         bclocal(0) = theta
         if ( theta .eq. 180.0d0 ) then
c           --- special case where we use K[ 1 + Cos(theta)]
            call twh_anglestyle(GLB_SET,type,ANG_STYLE_DREIDING)
c           --- set force constant to default
            bclocal(1) = 100.0d0*CON_KCALTOK
         else
c           --- general case with harmonic cosine
            call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARM_COS)
c           --- set force constant
            bclocal(1) = 0.5d0*100.0d0*CON_KCALTOK
     &           /(dsin(theta*onepi/180.0d0)**2)
         endif
c        --- set the rest of the information
         call twh_set_bencoeff(type,bclocal) 
         call twh_bendff(GLB_SET,type,ffname)
         nbslocal = 1
         bnlocal(1,1) = 'wild'
         bnlocal(1,2) = twh_get_caname(itype)
         bnlocal(1,3) = 'wild'
         call twh_set_bendnames(type,nbslocal,bnlocal)
      enddo

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * dihedrals                                                      *
c     ******************************************************************

      type = 0

      do itype = 1,nallatom
         code = twh_get_cbname(itype)
         lettera = code(3:3)
         if ( lettera .ne. ' ' .and. code(1:1) .ne. 'H' ) then
            do jtype = itype,nallatom
               code = twh_get_cbname(jtype)
               letterb = code(3:3)
               if ( letterb .ne. ' ' .and. code(1:1) .ne. 'H' ) then
c                 --- increase type (perhaps temporarily)
                  type = type + 1
c                 --- initialize lfound
                  lfound = .true.
c                 --- check to see if atoms are group 16
                  lgroupi = .false.
                  lgroupj = .false.
                  if ( twh_get_elementname(itype) .eq. ' O' .or.
     &                 twh_get_elementname(itype) .eq. ' S' .or.
     &                 twh_get_elementname(itype) .eq. 'Se' .or.
     &                 twh_get_elementname(itype) .eq. 'Te' .or.
     &                 twh_get_elementname(itype) .eq. 'Po' ) then
                     lgroupi = .true.
                  endif
                  if ( twh_get_elementname(jtype) .eq. ' O' .or.
     &                 twh_get_elementname(jtype) .eq. ' S' .or.
     &                 twh_get_elementname(jtype) .eq. 'Se' .or.
     &                 twh_get_elementname(jtype) .eq. 'Te' .or.
     &                 twh_get_elementname(jtype) .eq. 'Po' ) then
                     lgroupj = .true.
                  endif
                  if ( lettera .eq. '3' .and. letterb .eq. '3' ) then
                     if ( lgroupi .and. lgroupj ) then
c                       (h) A dihedral single bond involving two sp3
c                           atoms of the oxygen column
                        tclocal(1) = 2.0d0*(0.5d0)*CON_KCALTOK
                        tclocal(2) = 2.0d0
                        tclocal(3) = 90.0d0*(onepi/180.0d0)
                     else
c                       (a) A dihedral single bond involving two sp3 
c                           atoms
                        tclocal(1) = 2.0d0*(0.5d0)*CON_KCALTOK
                        tclocal(2) = 3.0d0
                        tclocal(3) = 180.0d0*(onepi/180.0d0)
                     endif
                  elseif ( 
     &                    (lettera .eq. '3' .and. letterb .eq. '2') .or.
     &                    (lettera .eq. '3' .and. letterb .eq. 'R') .or.
     &                    (lettera .eq. '2' .and. letterb .eq. '3') .or.
     &                    (lettera .eq. 'R' .and. letterb .eq. '3') 
     &                    ) then

                     if ( (lettera .eq. '3' .and. lgroupi) .or.
     &                    (letterb .eq. '3' .and. lgroupj) ) then
c                      (i) dihedral bonds involving an sp3 atom of
c                          the oxygen column with an sp2 or resonant 
c                          atom of another column
                        tclocal(1) = 2.0d0*(0.5d0)*CON_KCALTOK
                        tclocal(2) = 2.0d0
                        tclocal(3) = 180.0d0*(onepi/180.0d0)
                     else
c                       (b) A dihedral single bond involving one sp2 
c                           center and one sp3 center
                        tclocal(1) = 1.0d0*(0.5d0)*CON_KCALTOK
                        tclocal(2) = 6.0d0
                        tclocal(3) = 0.0d0*(onepi/180.0d0)
                     endif
                  elseif ( lettera .eq. '2' .and. letterb .eq. '2') then
c                    (c) A dihedral double bond involving two sp2 atoms
                     tclocal(1) = 45.0d0*(0.5d0)*CON_KCALTOK
                     tclocal(2) = 2.0d0
                     tclocal(3) = 180.0d0*(onepi/180.0d0)
                  elseif ( lettera .eq. 'R' .and. letterb .eq. 'R') then
c                    (d) A dihedral resonance bond
                     tclocal(1) = 25.0d0*(0.5d0)*CON_KCALTOK
                     tclocal(2) = 2.0d0
                     tclocal(3) = 180.0d0*(onepi/180.0d0)
                  else
c                    --- doesn't match with the general cases here
                     lfound = .false.
c                    --- decrement type
                     type = type - 1
                  endif
                  if ( lfound ) then
c                    --- set all of the logic
                     call twh_torstyle(GLB_SET,type,15)
                     tclocal(-1) = 1.0d0
                     call twh_loftor(GLB_SET,type,.true.)
                     call twh_set_torcoeff(type,tclocal)
                     call twh_torff(GLB_SET,type,ffname)
                     ntslocal = 1
                     tnlocal(1,1) = 'wild'
                     tnlocal(1,2) = twh_get_ctname(itype)
                     tnlocal(1,3) = twh_get_ctname(jtype)
                     tnlocal(1,4) = 'wild'
                     call twh_set_tornames(type,ntslocal,tnlocal)
                  endif
               endif
            enddo
         endif
      enddo

c     (j) an exception to the above is made for a single bond involving
c         one sp2 atom and one sp3 atom when the sp2 atom is bonded 
c         to another sp3 atom
c         I .ne. X_2,X_R
c         J = X_2,X_R
c         K = X_3
c     --- initialize isame to dummy value
      isame = MAXSAMET + 1
      do jtype = 1,nallatom
         code = twh_get_ctname(jtype)
         letterj = code(3:3)
         if ( letterj .eq. '2' .or. letterj .eq. 'R' ) then
            do ktype = 1,nallatom
               code = twh_get_ctname(ktype)
               letterk = code(3:3)
               if ( letterk .eq. '3' ) then
                  do itype = 1,nallatom
                     code = twh_get_ctname(itype)
                     letteri = code(3:3)
                     if ( letteri .ne. '2' .and. letteri .ne. 'R') then
                        if ( isame .ge. MAXSAMET ) then
                           type = type + 1
                           isame = 0
                        endif
c                       --- found one of our special cases
c                       V0 = 2.0 kcal/mol, n = 3, phi = 180
                        tclocal(1) = 2.0d0*(0.5d0)*CON_KCALTOK
                        tclocal(2) = 3.0d0
                        tclocal(3) = 180.0d0*(onepi/180.0d0)
                        call twh_torstyle(GLB_SET,type,15)
                        tclocal(-1) = 1.0d0
                        call twh_loftor(GLB_SET,type,.true.)
                        call twh_set_torcoeff(type,tclocal)
                        call twh_torff(GLB_SET,type,ffname)
                        isame = isame + 1
                        tnlocal(isame,1) = twh_get_ctname(itype)
                        tnlocal(isame,2) = twh_get_ctname(jtype)
                        tnlocal(isame,3) = twh_get_ctname(ktype)
                        tnlocal(isame,4) = 'wild'
c                       --- exception to our normal means of setting
c                       --- torsion names, need to do this by hand
                        call twh_ntsame(GLB_SET,type,isame)
                        do index = 1,4
                           call twh_tornames(GLB_SET,type,isame,index
     &                          ,tnlocal(isame,index))
                        enddo
                     endif
                  enddo
               endif
            enddo
         endif
      enddo

c     (g) dihedrals involving one or two sp1 atoms, monovalent atoms
c         or metals
c     i.e. everything else gets the null torsion
      type = type + 1
      call twh_torstyle(GLB_SET,type,8)
      tclocal(-1) = 1.0d0
      call twh_loftor(GLB_SET,type,.true.)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'wild'
      tnlocal(1,3) = 'wild'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- unable to implement the logic for (e) a dihedral single bond
c     --- provide the parameters here in case somebody wants to hard 
c     --- wire
      type = type + 1
      call twh_torstyle(GLB_SET,type,15)
      tclocal(-1) = 1.0d0
      tclocal(1) = 5.0d0*(0.5d0)*CON_KCALTOK
      tclocal(2) = 2.0d0
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_loftor(GLB_SET,type,.true.)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
      tnlocal(1,1) = 'special'
      tnlocal(1,2) = 'special'
      tnlocal(1,3) = 'special'
      tnlocal(1,4) = 'special'
      call twh_set_tornames(type,ntslocal,tnlocal)

c     --- unable to implement the logic for (f) exocyclic dihedral
c     --- provide the parameters here in case somebody wants to hard 
c     --- wire
      type = type + 1
      call twh_torstyle(GLB_SET,type,15)
      tclocal(-1) = 1.0d0
      tclocal(1) = 10.0d0*(0.5d0)*CON_KCALTOK
      tclocal(2) = 2.0d0
      tclocal(3) = 180.0d0*(onepi/180.0d0)
      call twh_loftor(GLB_SET,type,.true.)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
      tnlocal(1,1) = 'special'
      tnlocal(1,2) = 'special'
      tnlocal(1,3) = 'special'
      tnlocal(1,4) = 'special'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * impropers                                                      *
c     ******************************************************************

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,7)
      iclocal(1) =  (40.0d0)*(CON_KCALTOK)/(3.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0

c     --- improper exists for all of the _2 and _R atomtypes
      do itype = 1,nallatom
         code = twh_get_ctname(itype)
         lettera = code(3:3)
         if ( lettera .eq. '2' .or. lettera .eq. 'R' ) then
c           --- add this to the list
            isame = isame + 1
c           --- safety check isame
            if ( isame .gt. MAXSAMEIMP ) then
               write(6,*) 'FFDREIDING: MAXSAMEIMP exceeded:',MAXSAMEIMP
               lfailure = .true.
               return
            endif
c           --- assign the names
            ninlocal(isame,1) = twh_get_ctname(itype)
            ninlocal(isame,2) = 'wild'
            ninlocal(isame,3) = 'wild'
            ninlocal(isame,4) = 'wild'
         endif
      enddo
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return
      
      return
      end
