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
      subroutine twh_ffuff(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes,nimprotypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the uff          *
c     * forcefield                                                     *
c     * Note: UFF is a bit unusual in that it transforms its basic     *
c     * parameters into use with a variety of functional forms.        *
c     *                                                                *
c     * originally written 12-01-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#define FUNCTION_GET_BONDPATT
#define FUNCTION_GET_CANAME
#define FUNCTION_GET_CBNAME
#define FUNCTION_GET_CLASSICAL_POTENTIAL
#define FUNCTION_GET_CTNAME
#define FUNCTION_GET_ELEMENTNAME
#define FUNCTION_GET_NBCOEFF
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes,nimprotypes
c     --- local variables
      logical lfound,lgsixj,lgsixk
      character*1 lettera,letterb,letterc,letterj,letterk
      character*(FFNAMELEN) code,ffname
      integer nvslocal,nbslocal,ntslocal,icoeff
      character*(FFNAMELEN) vnlocal
      dimension vnlocal(MAXSAMEV,MAXVNAME) 
      character*(FFNAMELEN) bnlocal 
      dimension bnlocal(MAXSAMEB,MAXBNAME)
      character*(FFNAMELEN) tnlocal
      dimension tnlocal(MAXSAMET,MAXTNAME)
      double precision vclocal
      dimension vclocal(MINVIBCOEFF:MAXVIBCOEFF)
      integer ii,jj,kk,cc,isame,index
      integer type,itype,jtype,ctype,ktype
      double precision sixroottwo,bondorder,lamda,ri,rj,rbo,ren
      double precision theta,vj,vk,phi
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

      ffname = 'UFF'
      call twh_classical_potential(GLB_SET,'UFF 12-6')
      call twh_classical_mixrule(GLB_SET,'Geometric')
c     --- initialize types
      type = 0

c     ******************************************************************
c     * UFF                                                            *
c     ******************************************************************
c     References:
c     [jacs1992]: A.K. Rappe; C.J. Casewit; K.S. Colwell; 
c     W.A. Goddard III; W.M. Skiff;
c     "UFF, a Full Periodic Table Force Field for Molecular Mechanics
c     and Molecular Dynamics Simulations", J. Am. Chem. Soc. 114,
c     10024-10035 (1992)

c     --- setup all of the basic UFF parameters for later conversion
c     --- into the actual forcefield parameters
      
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.354d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 2.886d0
      nbclocal(8) = 0.044d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 0.712d0
      nbclocal(11) = 4.5280d0
c     --- no default value provided for period 1, setting to 0.0
      nbclocal(12) = 0.0d0
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
      nbclocal(5) = 0.460d0
      nbclocal(6) = 83.5d0
      nbclocal(7) = 2.886d0
      nbclocal(8) = 0.044d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 0.712d0
      nbclocal(11) = 4.5280d0
c     --- no default value provided for period 1, setting to 0.0
      nbclocal(12) = 0.0d0
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
      nbclocal(5) = 0.849d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.362d0
      nbclocal(8) = 0.056d0
      nbclocal(9) = 15.24d0
      nbclocal(10) = 0.098d0
      nbclocal(11) = 9.66d0
c     --- no default value provided for period 1, setting to 0.0
      nbclocal(12) = 0.0d0
      masslocal = 4.003d0
      elelocal = 'He'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='He4+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.336d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 2.451d0
      nbclocal(8) = 0.025d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.026d0
      nbclocal(11) = 3.006d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 6.941d0
      elelocal = 'Li'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Li'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.074d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 2.745d0
      nbclocal(8) = 0.085d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.565d0
      nbclocal(11) = 4.877d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 9.012d0
      elelocal = 'Be'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Be3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.838d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.083d0
      nbclocal(8) = 0.180d0
      nbclocal(9) = 12.052d0
      nbclocal(10) = 1.755d0
      nbclocal(11) = 5.11d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
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
      nbclocal(5) = 0.828d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 4.083d0
      nbclocal(8) = 0.180d0
      nbclocal(9) = 12.052d0
      nbclocal(10) = 1.755d0
      nbclocal(11) = 5.11d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 10.811d0
      elelocal = ' B'
      bplocal = 'sp2'
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
      nbclocal(5) = 0.757d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.851d0
      nbclocal(8) = 0.105d0
      nbclocal(9) = 12.73d0
      nbclocal(10) = 1.912d0
      nbclocal(11) = 5.343d0
c     --- this value is directly specified in table III
      nbclocal(12) = 2.119d0
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
      nbclocal(5) = 0.729d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.851d0
      nbclocal(8) = 0.105d0
      nbclocal(9) = 12.73d0
      nbclocal(10) = 1.912d0
      nbclocal(11) = 5.343d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
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
      nbclocal(5) = 0.732d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.851d0
      nbclocal(8) = 0.105d0
      nbclocal(9) = 12.73d0
      nbclocal(10) = 1.912d0
      nbclocal(11) = 5.343d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
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
      nbclocal(5) = 0.706d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.851d0
      nbclocal(8) = 0.105d0
      nbclocal(9) = 12.73d0
      nbclocal(10) = 1.912d0
      nbclocal(11) = 5.343d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'null'
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
      nbclocal(5) = 0.700d0
      nbclocal(6) = 106.7d0
      nbclocal(7) = 3.660d0
      nbclocal(8) = 0.069d0
      nbclocal(9) = 13.407d0
      nbclocal(10) = 2.544d0
      nbclocal(11) = 6.899d0
c     --- specified directly in table III
      nbclocal(12) = 0.450d0
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
      nbclocal(5) = 0.699d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.660d0
      nbclocal(8) = 0.069d0
      nbclocal(9) = 13.407d0
      nbclocal(10) = 2.544d0
      nbclocal(11) = 6.899d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
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
      nbclocal(5) = 0.685d0
      nbclocal(6) = 111.2d0
      nbclocal(7) = 3.660d0
      nbclocal(8) = 0.069d0
      nbclocal(9) = 13.407d0
      nbclocal(10) = 2.544d0
      nbclocal(11) = 6.899d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
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
      nbclocal(5) = 0.656d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.660d0
      nbclocal(8) = 0.069d0
      nbclocal(9) = 13.407d0
      nbclocal(10) = 2.544d0
      nbclocal(11) = 6.899d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 14.007d0
      elelocal = ' N'
      bplocal = 'null'
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
      nbclocal(5) = 0.658d0
      nbclocal(6) = 104.51d0
      nbclocal(7) = 3.500d0
      nbclocal(8) = 0.060d0
      nbclocal(9) = 14.085d0
      nbclocal(10) = 2.300d0
      nbclocal(11) = 8.741d0
c     --- directly specified in table III
      nbclocal(12) = 0.018d0
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
      nbclocal(5) = 0.528d0
      nbclocal(6) = 146.0d0
      nbclocal(7) = 3.500d0
      nbclocal(8) = 0.060d0
      nbclocal(9) = 14.085d0
      nbclocal(10) = 2.300d0
      nbclocal(11) = 8.741d0
c     --- directly specified in table III
      nbclocal(12) = 0.018d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='O_3_z'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.680d0
      nbclocal(6) = 110.0d0
      nbclocal(7) = 3.500d0
      nbclocal(8) = 0.060d0
      nbclocal(9) = 14.085d0
      nbclocal(10) = 2.300d0
      nbclocal(11) = 8.741d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 's'
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
      nbclocal(5) = 0.634d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 3.500d0
      nbclocal(8) = 0.060d0
      nbclocal(9) = 14.085d0
      nbclocal(10) = 2.300d0
      nbclocal(11) = 8.741d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'null'
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
      nbclocal(5) = 0.639d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.500d0
      nbclocal(8) = 0.060d0
      nbclocal(9) = 14.085d0
      nbclocal(10) = 2.300d0
      nbclocal(11) = 8.741d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 15.999d0
      elelocal = ' O'
      bplocal = 'null'
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
      nbclocal(5) = 0.668d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.364d0
      nbclocal(8) = 0.050d0
      nbclocal(9) = 14.762d0
      nbclocal(10) = 1.735d0
      nbclocal(11) = 10.874d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
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
      nbclocal(5) = 0.920d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.243d0
      nbclocal(8) = 0.042d0
      nbclocal(9) = 15.440d0
      nbclocal(10) = 0.194d0
      nbclocal(11) = 11.04d0
c     --- default value for period 2 is 2.0
      nbclocal(12) = 2.0d0
      masslocal = 20.180d0
      elelocal = 'Ne'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ne4+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.539d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 2.983d0
      nbclocal(8) = 0.030d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.081d0
      nbclocal(11) = 2.843d0
c     --- default value for period 3 is 1.25
      nbclocal(12) = 1.25d0
      masslocal = 22.990d0
      elelocal = 'Na'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Na'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.421d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.021d0
      nbclocal(8) = 0.111d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.787d0
      nbclocal(11) = 3.951d0
c     --- default value for period 3 is 1.25
      nbclocal(12) = 1.25d0
      masslocal = 24.305d0
      elelocal = 'Mg'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Mg3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.244d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.499d0
      nbclocal(8) = 0.505d0
      nbclocal(9) = 11.278d0
      nbclocal(10) = 1.792d0
      nbclocal(11) = 4.06d0
c     --- default value for period 3 is 1.25
      nbclocal(12) = 1.25d0
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
      nbclocal(5) = 1.117d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.295d0
      nbclocal(8) = 0.402d0
      nbclocal(9) = 12.175d0
      nbclocal(10) = 2.323d0
      nbclocal(11) = 4.168d0
c     --- directly specified in table III
      nbclocal(12) = 1.225d0
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
      nbclocal(5) = 1.101d0
      nbclocal(6) = 93.8d0
      nbclocal(7) = 4.147d0
      nbclocal(8) = 0.305d0
      nbclocal(9) = 13.072d0
      nbclocal(10) = 2.863d0
      nbclocal(11) = 5.463d0
c     --- directly specified in table III
      nbclocal(12) = 2.400d0
      masslocal = 30.974d0
      elelocal = ' P'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P_3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.056d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.147d0
      nbclocal(8) = 0.305d0
      nbclocal(9) = 13.072d0
      nbclocal(10) = 2.863d0
      nbclocal(11) = 5.463d0
c     --- directly specified in table III
      nbclocal(12) = 2.400d0
      masslocal = 30.974d0
      elelocal = ' P'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P_3+5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.056d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.147d0
      nbclocal(8) = 0.305d0
      nbclocal(9) = 13.072d0
      nbclocal(10) = 2.863d0
      nbclocal(11) = 5.463d0
c     --- directly specified in table III
      nbclocal(12) = 2.400d0
      masslocal = 30.974d0
      elelocal = ' P'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='P_3+q'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.064d0
      nbclocal(6) = 92.1d0
      nbclocal(7) = 4.035d0
      nbclocal(8) = 0.274d0
      nbclocal(9) = 13.969d0
      nbclocal(10) = 2.703d0
      nbclocal(11) = 6.928d0
c     --- directly specified in table III
      nbclocal(12) = 0.484d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S_3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.049d0
      nbclocal(6) = 103.20d0
      nbclocal(7) = 4.035d0
      nbclocal(8) = 0.274d0
      nbclocal(9) = 13.969d0
      nbclocal(10) = 2.703d0
      nbclocal(11) = 6.928d0
c     --- directly specified in table III
      nbclocal(12) = 0.484d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S_3+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.027d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.035d0
      nbclocal(8) = 0.274d0
      nbclocal(9) = 13.969d0
      nbclocal(10) = 2.703d0
      nbclocal(11) = 6.928d0
c     --- directly specified in table III
      nbclocal(12) = 0.484d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S_3+6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.077d0
      nbclocal(6) = 92.2d0
      nbclocal(7) = 4.035d0
      nbclocal(8) = 0.274d0
      nbclocal(9) = 13.969d0
      nbclocal(10) = 2.703d0
      nbclocal(11) = 6.928d0
c     --- default value for period 3 is 1.25
      nbclocal(12) = 1.25d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S_R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 0.854d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 4.035d0
      nbclocal(8) = 0.274d0
      nbclocal(9) = 13.969d0
      nbclocal(10) = 2.703d0
      nbclocal(11) = 6.928d0
c     --- default for period 3 is 1.25
      nbclocal(12) = 1.25d0
      masslocal = 32.066d0
      elelocal = ' S'
      bplocal = 'sp2'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='S_2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.044d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.947d0
      nbclocal(8) = 0.227d0
      nbclocal(9) = 14.866d0
      nbclocal(10) = 2.348d0
      nbclocal(11) = 8.564d0
c     --- default value for period 3 is 1.25
      nbclocal(12) = 1.25d0
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
      nbclocal(5) = 1.032d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.868d0
      nbclocal(8) = 0.185d0
      nbclocal(9) = 15.763d0
      nbclocal(10) = 0.300d0
      nbclocal(11) = 9.465d0
c     --- default value for period 3 is 1.25
      nbclocal(12) = 1.25d0
      masslocal = 39.948d0
      elelocal = 'Ar'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ar4+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.953d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.812d0
      nbclocal(8) = 0.035d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.165d0
      nbclocal(11) = 2.421d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 39.098d0
      elelocal = ' K'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='K_'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.761d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.399d0
      nbclocal(8) = 0.238d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.141d0
      nbclocal(11) = 3.231d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 40.08d0
      elelocal = 'Ca'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ca6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.513d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.295d0
      nbclocal(8) = 0.019d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.592d0
      nbclocal(11) = 3.395d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 44.956d0
      elelocal = 'Sc'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sc3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.412d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.175d0
      nbclocal(8) = 0.017d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.659d0
      nbclocal(11) = 3.47d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 47.88d0
      elelocal = 'Ti'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ti3+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.412d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.175d0
      nbclocal(8) = 0.017d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.659d0
      nbclocal(11) = 3.47d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 47.88d0
      elelocal = 'Ti'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ti6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.402d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.144d0
      nbclocal(8) = 0.016d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.679d0
      nbclocal(11) = 3.65d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 50.942d0
      elelocal = ' V'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='V_3+5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.345d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.023d0
      nbclocal(8) = 0.015d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.463d0
      nbclocal(11) = 3.415d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 51.996d0
      elelocal = 'Cr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cr6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.382d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.961d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.43d0
      nbclocal(11) = 3.325d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 54.938d0
      elelocal = 'Mn'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Mn6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.270d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 2.912d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.43d0
      nbclocal(11) = 3.76d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 55.847d0
      elelocal = 'Fe'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fe3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.335d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.912d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.43d0
      nbclocal(11) = 3.76d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 55.847d0
      elelocal = 'Fe'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fe6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.241d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.872d0
      nbclocal(8) = 0.014d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.43d0
      nbclocal(11) = 4.105d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 58.933d0
      elelocal = 'Co'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Co6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.164d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.834d0
      nbclocal(8) = 0.015d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.43d0
      nbclocal(11) = 4.465d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 58.69d0
      elelocal = 'Ni'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ni4+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.302d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.495d0
      nbclocal(8) = 0.005d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.756d0
      nbclocal(11) = 3.729d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 63.546d0
      elelocal = 'Cu'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cu3+1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.193d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 2.763d0
      nbclocal(8) = 0.124d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.308d0
      nbclocal(11) = 5.106d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 65.39d0
      elelocal = 'Zn'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Zn3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.260d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.383d0
      nbclocal(8) = 0.415d0
      nbclocal(9) = 11.0d0
      nbclocal(10) = 1.821d0
      nbclocal(11) = 3.641d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 69.723d0
      elelocal = 'Ga'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ga3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.197d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.280d0
      nbclocal(8) = 0.379d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.789d0
      nbclocal(11) = 4.051d0
c     --- directly specified in table III
      nbclocal(12) = 0.701d0
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
      nbclocal(5) = 1.211d0
      nbclocal(6) = 92.1d0
      nbclocal(7) = 4.230d0
      nbclocal(8) = 0.309d0
      nbclocal(9) = 13.0d0
      nbclocal(10) = 2.864d0
      nbclocal(11) = 5.188d0
c     --- directly specified in table III
      nbclocal(12) = 1.5d0
      masslocal = 74.922d0
      elelocal = 'As'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='As3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.190d0
      nbclocal(6) = 90.6d0
      nbclocal(7) = 4.205d0
      nbclocal(8) = 0.291d0
      nbclocal(9) = 14.0d0
      nbclocal(10) = 2.764d0
      nbclocal(11) = 6.428d0
c     --- directly specified in table III
      nbclocal(12) = 0.335d0
      masslocal = 78.96d0
      elelocal = 'Se'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Se3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.192d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.189d0
      nbclocal(8) = 0.251d0
      nbclocal(9) = 15.0d0
      nbclocal(10) = 2.519d0
      nbclocal(11) = 7.790d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
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
      nbclocal(5) = 1.147d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 4.141d0
      nbclocal(8) = 0.220d0
      nbclocal(9) = 16.0d0
      nbclocal(10) = 0.452d0
      nbclocal(11) = 8.505d0
c     --- default value for period 4 is 0.7
      nbclocal(12) = 0.7d0
      masslocal = 83.80d0
      elelocal = 'Kr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Kr4+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 2.260d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.114d0
      nbclocal(8) = 0.04d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.592d0
      nbclocal(11) = 2.331d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 85.47d0
      elelocal = 'Rb'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Rb'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 2.052d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.641d0
      nbclocal(8) = 0.235d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.449d0
      nbclocal(11) = 3.024d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 87.62d0
      elelocal = 'Sr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sr6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.698d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.345d0
      nbclocal(8) = 0.072d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.257d0
      nbclocal(11) = 3.83d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 88.906d0
      elelocal = ' Y'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Y_3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.564d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.124d0
      nbclocal(8) = 0.069d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.667d0
      nbclocal(11) = 3.40d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 91.224d0
      elelocal = 'Zr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Zr3+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.473d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.165d0
      nbclocal(8) = 0.059d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.618d0
      nbclocal(11) = 3.55d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 92.906d0
      elelocal = 'Nb'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Nb3+5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.467d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.052d0
      nbclocal(8) = 0.056d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.40d0
      nbclocal(11) = 3.465d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 95.94d0
      elelocal = 'Mo'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Mo6+6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.484d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.052d0
      nbclocal(8) = 0.056d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.40d0
      nbclocal(11) = 3.465d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 95.94d0
      elelocal = 'Mo'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Mo3+6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.322d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.998d0
      nbclocal(8) = 0.048d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.40d0
      nbclocal(11) = 3.29d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 98.d0
      elelocal = 'Tc'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Tc6+5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.478d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.963d0
      nbclocal(8) = 0.056d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.40d0
      nbclocal(11) = 3.575d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 101.07d0
      elelocal = 'Ru'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ru6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.332d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.929d0
      nbclocal(8) = 0.053d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.508d0
      nbclocal(11) = 3.975d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 102.91d0
      elelocal = 'Rh'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Rh6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.338d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.899d0
      nbclocal(8) = 0.048d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.21d0
      nbclocal(11) = 4.32d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 106.42d0
      elelocal = 'Pd'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pd4+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.386d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 3.148d0
      nbclocal(8) = 0.036d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.956d0
      nbclocal(11) = 4.439d0
      nbclocal(11) = 4.436d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 107.87d0
      elelocal = 'Ag'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ag1+1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.403d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 2.848d0
      nbclocal(8) = 0.228d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.65d0
      nbclocal(11) = 5.034d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 112.41d0
      elelocal = 'Cd'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cd3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.459d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.463d0
      nbclocal(8) = 0.599d0
      nbclocal(9) = 11.0d0
      nbclocal(10) = 2.07d0
      nbclocal(11) = 3.506d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 114.82d0
      elelocal = 'In'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='In3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.398d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.392d0
      nbclocal(8) = 0.567d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.961d0
      nbclocal(11) = 3.987d0
c     --- directly specified in table III
      nbclocal(12) = 0.199d0
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
      nbclocal(5) = 1.407d0
      nbclocal(6) = 91.6d0
      nbclocal(7) = 4.420d0
      nbclocal(8) = 0.449d0
      nbclocal(9) = 13.0d0
      nbclocal(10) = 2.704d0
      nbclocal(11) = 4.899d0
c     --- directly specified in table III
      nbclocal(12) = 1.1d0
      masslocal = 121.75d0
      elelocal = 'Sb'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sb3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.386d0
      nbclocal(6) = 90.25d0
      nbclocal(7) = 4.470d0
      nbclocal(8) = 0.398d0
      nbclocal(9) = 14.0d0
      nbclocal(10) = 2.882d0
      nbclocal(11) = 5.816d0
c     --- directly specified in table III
      nbclocal(12) = 0.3d0
      masslocal = 127.60d0
      elelocal = 'Te'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Te3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.382d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.50d0
      nbclocal(8) = 0.339d0
      nbclocal(9) = 15.0d0
      nbclocal(10) = 2.65d0
      nbclocal(11) = 6.822d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
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
      nbclocal(5) = 1.267d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 4.404d0
      nbclocal(8) = 0.332d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 0.556d0
      nbclocal(11) = 7.595d0
c     --- default value for period 5 is 0.2
      nbclocal(12) = 0.2d0
      masslocal = 131.29d0
      elelocal = 'Xe'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Xe4+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 2.570d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.517d0
      nbclocal(8) = 0.045d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.573d0
      nbclocal(11) = 2.183d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 132.90d0
      elelocal = 'Cs'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cs'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 2.277d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.703d0
      nbclocal(8) = 0.364d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.727d0
      nbclocal(11) = 2.814d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 137.33d0
      elelocal = 'Ba'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ba6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.943d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.522d0
      nbclocal(8) = 0.017d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.8355d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 138.91d0
      elelocal = 'La'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='La3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.841d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.556d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.774d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 140.12d0
      elelocal = 'Ce'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ce6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.823d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.606d0
      nbclocal(8) = 0.010d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.858d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 140.91d0
      elelocal = 'Pr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pr6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.816d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.575d0
      nbclocal(8) = 0.010d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.8685d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 144.24d0
      elelocal = 'Nd'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Nd6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.801d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.547d0
      nbclocal(8) = 0.009d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.881d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 145.0d0
      elelocal = 'Pm'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pm6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.780d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.520d0
      nbclocal(8) = 0.008d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.9115d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 150.36d0
      elelocal = 'Sm'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Sm6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.771d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.493d0
      nbclocal(8) = 0.008d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 2.8785d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 151.96d0
      elelocal = 'Eu'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Eu6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.735d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.368d0
      nbclocal(8) = 0.009d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 3.1665d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 157.25d0
      elelocal = 'Gd'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Gd6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.732d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.451d0
      nbclocal(8) = 0.007d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 3.018d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 158.92d0
      elelocal = 'Tb'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Tb6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.710d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.428d0
      nbclocal(8) = 0.007d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 3.0555d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 162.50d0
      elelocal = 'Dy'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Dy6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.696d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.409d0
      nbclocal(8) = 0.007d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.416d0
      nbclocal(11) = 3.127d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 164.93d0
      elelocal = 'Ho'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ho6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.673d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.391d0
      nbclocal(8) = 0.007d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 3.1865d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 167.26d0
      elelocal = 'Er'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Er6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.660d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.374d0
      nbclocal(8) = 0.006d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.30d0
      nbclocal(11) = 3.2514d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 168.93d0
      elelocal = 'Tm'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Tm6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.637d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.355d0
      nbclocal(8) = 0.228d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.618d0
      nbclocal(11) = 3.2889d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 173.04d0
      elelocal = 'Yb'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Yb6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.671d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.640d0
      nbclocal(8) = 0.041d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.271d0
      nbclocal(11) = 2.9629d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 174.97d0
      elelocal = 'Lu'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Lu6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.611d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.141d0
      nbclocal(8) = 0.072d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.921d0
      nbclocal(11) = 3.70d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 178.49d0
      elelocal = 'Hf'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Hf3+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.511d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.170d0
      nbclocal(8) = 0.081d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 4.075d0
      nbclocal(11) = 5.10d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 180.95d0
      elelocal = 'Ta'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ta3+5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.392d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.069d0
      nbclocal(8) = 0.067d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.70d0
      nbclocal(11) = 4.63d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 183.85d0
      elelocal = ' W'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='W_6+6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.526d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.069d0
      nbclocal(8) = 0.067d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.70d0
      nbclocal(11) = 4.63d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 183.85d0
      elelocal = ' W'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='W_3+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.380d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 3.069d0
      nbclocal(8) = 0.067d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.70d0
      nbclocal(11) = 4.63d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 183.85d0
      elelocal = ' W'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='W_3+6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.372d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.954d0
      nbclocal(8) = 0.066d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.70d0
      nbclocal(11) = 3.96d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 186.21d0
      elelocal = 'Re'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Re6+5'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.314d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 2.954d0
      nbclocal(8) = 0.066d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.70d0
      nbclocal(11) = 3.96d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 186.21d0
      elelocal = 'Re'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Re3+7'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.372d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.120d0
      nbclocal(8) = 0.037d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.70d0
      nbclocal(11) = 5.14d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 190.2d0
      elelocal = 'Os'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Os6+6'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.371d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.840d0
      nbclocal(8) = 0.073d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.731d0
      nbclocal(11) = 5.00d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 192.22d0
      elelocal = 'Ir'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ir6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.364d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 2.754d0
      nbclocal(8) = 0.080d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.382d0
      nbclocal(11) = 4.79d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 195.08d0
      elelocal = 'Pt'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pt4+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.262d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.293d0
      nbclocal(8) = 0.039d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.625d0
      nbclocal(11) = 4.894d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 196.97d0
      elelocal = 'Au'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Au4+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.340d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 2.705d0
      nbclocal(8) = 0.385d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.75d0
      nbclocal(11) = 6.27d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 200.59d0
      elelocal = 'Hg'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Hg1+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.518d0
      nbclocal(6) = 120.0d0
      nbclocal(7) = 4.347d0
      nbclocal(8) = 0.680d0
      nbclocal(9) = 11.0d0
      nbclocal(10) = 2.068d0
      nbclocal(11) = 3.20d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 204.38d0
      elelocal = 'Tl'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Tl3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.459d0
      nbclocal(6) = 109.47d0
      nbclocal(7) = 4.297d0
      nbclocal(8) = 0.663d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.846d0
      nbclocal(11) = 3.90d0
c     --- directly specified in table III
      nbclocal(12) = 0.1d0
      masslocal = 207.2d0
      elelocal = 'Pb'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pb3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.512d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 4.370d0
      nbclocal(8) = 0.518d0
      nbclocal(9) = 13.0d0
      nbclocal(10) = 2.470d0
      nbclocal(11) = 4.69d0
c     --- directly specified in table III
      nbclocal(12) = 1.0d0
      masslocal = 208.98d0
      elelocal = 'Bi'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Bi3+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.50d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 4.709d0
      nbclocal(8) = 0.325d0
      nbclocal(9) = 14.0d0
      nbclocal(10) = 2.33d0
      nbclocal(11) = 4.21d0
c     --- directly specified in table III
      nbclocal(12) = 0.3d0
      masslocal = 209.0d0
      elelocal = 'Po'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Po3+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.545d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.750d0
      nbclocal(8) = 0.284d0
      nbclocal(9) = 15.0d0
      nbclocal(10) = 2.24d0
      nbclocal(11) = 4.75d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 210.0d0
      elelocal = 'At'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='At'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.420d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 4.765d0
      nbclocal(8) = 0.248d0
      nbclocal(9) = 16.0d0
      nbclocal(10) = 0.583d0
      nbclocal(11) = 5.37d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 222.0d0
      elelocal = 'Rn'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Rn4+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 2.880d0
      nbclocal(6) = 180.0d0
      nbclocal(7) = 4.90d0
      nbclocal(8) = 0.050d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 1.847d0
      nbclocal(11) = 2.00d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 223.0d0
      elelocal = 'Fr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fr'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 2.512d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.677d0
      nbclocal(8) = 0.404d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 2.92d0
      nbclocal(11) = 2.843d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 226.0d0
      elelocal = 'Ra'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ra6+2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.983d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.478d0
      nbclocal(8) = 0.033d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 2.835d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 227.0d0
      elelocal = 'Ac'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Ac6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.721d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.396d0
      nbclocal(8) = 0.026d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 4.202d0
      nbclocal(11) = 3.175d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 232.04d0
      elelocal = 'Th'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Th6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.711d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.424d0
      nbclocal(8) = 0.022d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 2.985d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 231.0d0
      elelocal = 'Pa'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pa6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.684d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.395d0
      nbclocal(8) = 0.022d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.341d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 238.03d0
      elelocal = ' U'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='U_6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.666d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.424d0
      nbclocal(8) = 0.019d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.549d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 237.0d0
      elelocal = 'Np'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Np6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.657d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.424d0
      nbclocal(8) = 0.016d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.243d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 244.0d0
      elelocal = 'Pu'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Pu6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.660d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.381d0
      nbclocal(8) = 0.014d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 2.9895d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 243.0d0
      elelocal = 'Am'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Am6+4'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.801d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.326d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 2.8315d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 247.0d0
      elelocal = 'Cm'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cm6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.761d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.339d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.1935d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 247.0d0
      elelocal = 'Bk'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Bk6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.750d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.313d0
      nbclocal(8) = 0.013d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.197d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 251.0d0
      elelocal = 'Cf'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Cf6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.724d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.299d0
      nbclocal(8) = 0.012d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.333d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 252.0d0
      elelocal = 'Es'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Es6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.712d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.286d0
      nbclocal(8) = 0.012d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.40d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 257.0d0
      elelocal = 'Fm'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Fm6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.689d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.274d0
      nbclocal(8) = 0.011d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.47d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 258.0d0
      elelocal = 'Md'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Md6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.679d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.248d0
      nbclocal(8) = 0.011d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.475d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 259.0d0
      elelocal = 'No'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='No6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(5) = 1.698d0
      nbclocal(6) = 90.0d0
      nbclocal(7) = 3.236d0
      nbclocal(8) = 0.011d0
      nbclocal(9) = 12.0d0
      nbclocal(10) = 3.90d0
      nbclocal(11) = 3.500d0
c     --- default value for period 6 is 0.1
      nbclocal(12) = 0.1d0
      masslocal = 260.0d0
      elelocal = 'Lr'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='Lr6+3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * transform into nonbond                                         *
c     ******************************************************************

      sixroottwo = 2.0d0**(1.0d0/6.0d0)
      do itype = 1,ntypes
         ii = twh_getnbtype(itype,itype)
         if ( twh_get_classical_potential() .eq. 'UFF 12-6' ) then
c           --- convert from Xi units to sigma
            nbclocal(1) = twh_get_nbcoeff(7,ii)/sixroottwo
c           --- convert Di from kcal to kelvin
            nbclocal(2) = twh_get_nbcoeff(8,ii)*CON_KCALTOK
c           --- UFF includes full VDW and coulombic interactions for 1-4
            nbclocal(3) = nbclocal(1)
            nbclocal(4) = nbclocal(2)
c           --- set manually so we do not overwrite the variables that
c           --- were already set above
            do icoeff = 1,4
               call twh_nbcoeff(GLB_SET,icoeff,ii,nbclocal(icoeff))
            enddo
         elseif ( twh_get_classical_potential() 
     &           .eq. 'UFF Exponential-6' ) then
c           --- assign the UFF C6 term and convert from kcal
            nbclocal(1) = CON_KCALTOK*twh_get_nbcoeff(8,ii)
     &           *(twh_get_nbcoeff(7,ii)**6)
     &           *( twh_get_nbcoeff(9,ii)
     &           / (twh_get_nbcoeff(9,ii)-6.0d0) )
c           --- assign the UFF A term (prefactor) convert from kcal
            nbclocal(2) = CON_KCALTOK*twh_get_nbcoeff(8,ii)
     &           *(6.0d0/(twh_get_nbcoeff(9,ii)-6.0d0))
     &           *dexp(twh_get_nbcoeff(9,ii))
c           --- assign the UFF B term
            nbclocal(3) = -twh_get_nbcoeff(9,ii)/twh_get_nbcoeff(7,ii)
c           --- set manually so we do not overwrite the variables that
c           --- were already set above
            do icoeff = 1,3
               call twh_nbcoeff(GLB_SET,icoeff,ii,nbclocal(icoeff))
            enddo
         else
            write(6,*) 'FFUFF: unknown classical_potential'
            write(6,*) twh_get_classical_potential()
            lfailure = .true.
            return
         endif
      enddo

c     ******************************************************************
c     * bonds                                                          *
c     ******************************************************************

c     --- initialize type
      type = 0
c     --- set some bond constants
      lamda = 0.1332d0
c     --- assign all of the bonds
      do itype = 1,ntypes
         ii = twh_getnbtype(itype,itype)
         do jtype = itype,ntypes
            jj = twh_getnbtype(jtype,jtype)
c           --- determine bond order
            code = twh_get_cbname(itype)
            lettera = code(3:3)
            code = twh_get_cbname(jtype)
            letterb = code(3:3)
            if ( twh_get_cbname(itype) .eq. 'C_R' .and. 
     &           twh_get_cbname(jtype) .eq. 'N_R' ) then
c              --- special hack for C-N amide bond
               bondorder = 1.41d0
            elseif ( lettera .eq.'R' .and. letterb .eq. 'R' ) then
c              --- aromatic bond
               bondorder = 1.5d0
            elseif ( lettera .eq. '1' .and. letterb .eq. '1' ) then
c              --- triple bond
               bondorder = 3.0d0
            elseif ( lettera .eq. '2' .and. letterb .eq. '2' ) then
c              --- double bond
               bondorder = 2.0d0
            else
c              --- any other combination is a single bond (for now)
               bondorder = 1.0d0
            endif
c           --- increment type
            type = type + 1
c           --- harmonic
            call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
c           --- compute the equilibrium bond length
c           --- store ri
            ri = twh_get_nbcoeff(5,ii)
c           --- store rj
            rj = twh_get_nbcoeff(5,jj)
c           --- calculate r_BO
            rbo = - lamda * ( ri + rj ) * dlog(bondorder)
c           --- r_EN
            ren = ri*rj
     &           *(dsqrt(twh_get_nbcoeff(11,ii))
     &           - dsqrt(twh_get_nbcoeff(11,jj)))**2
     &           / (twh_get_nbcoeff(11,ii)*ri
     &           + twh_get_nbcoeff(11,jj)*rj)
            vclocal(0) = ri + rj + rbo - ren
c           --- compute force constant
            vclocal(1) = 
     &            CON_KCALTOK*(0.5d0)*
     &           (664.12d0)*twh_get_nbcoeff(10,ii)
     &           *twh_get_nbcoeff(10,jj)/(vclocal(0)**3)
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

      do ctype = 1,ntypes
c        --- determine whether this atom can be at the center of an
c        --- angle via the 3rd character of the code
         code = twh_get_caname(ctype)
         letterc = code(3:3)
         cc = twh_getnbtype(ctype,ctype)
         if ( letterc .ne. ' ' ) then
c           --- viable angle center found, run through all combinations
c           --- increment angle counter
            type = type + 1
c           --- equilibrium bond angle depends only on center
            bclocal(0) = twh_get_nbcoeff(6,cc)
c           --- determine type of interaction based on center
            if ( letterc .eq. '1' ) then
               call twh_anglestyle(GLB_SET,type,ANG_STYLE_FOURIER)
               bclocal(1) = 1.0d0
               bclocal(2) = 1.0d0
            elseif ( letterc .eq. '2' ) then
               call twh_anglestyle(GLB_SET,type,ANG_STYLE_FOURIER)
               bclocal(1) = 3.0d0
               bclocal(2) = -1.0d0
            elseif ( letterc .eq. '3' 
     &              .and. bclocal(0) .eq. 90.0d0 ) then
               call twh_anglestyle(GLB_SET,type,ANG_STYLE_FOURIER)
               bclocal(1) = 2.0d0
               bclocal(2) = 1.0d0
            elseif ( letterc .eq. '4' ) then
               call twh_anglestyle(GLB_SET,type,ANG_STYLE_FOURIER)
               bclocal(1) = 4.0d0
               bclocal(2) = -1.0d0
            elseif ( letterc .eq. '6' ) then
               call twh_anglestyle(GLB_SET,type,ANG_STYLE_FOURIER)
               bclocal(1) = 4.0d0
               bclocal(2) = -1.0d0
            else
c              --- general nonlinear case
               call twh_anglestyle(GLB_SET,type,ANG_STYLE_TTFOURIER)
               theta = bclocal(0)*onepi/180.0d0
               bclocal(3) = 1.0d0/(4.0d0*(dsin(theta)**2))
               bclocal(2) = -4.0d0*bclocal(3)*dcos(theta)
               bclocal(1) = bclocal(3)
     &              *(2.0d0*(dcos(theta)**2)+1.0d0)
            endif
c           --- set lookup information about this angle
            call twh_set_bencoeff(type,bclocal) 
            call twh_bendff(GLB_SET,type,ffname)
            nbslocal = 1
            bnlocal(1,1) = 'wild'
            bnlocal(1,2) = twh_get_caname(ctype)
            bnlocal(1,3) = 'wild'
            call twh_set_bendnames(type,nbslocal,bnlocal)
         endif
      enddo

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * dihedrals                                                      *
c     ******************************************************************

      type = 0

      do jtype = 1,ntypes
         code = twh_get_ctname(jtype)
         letterj = code(3:3)
         jj = twh_getnbtype(jtype,jtype)
c        --- determine whether this is a group 6 sp3 atom
         lgsixj = .false.
         if ( twh_get_elementname(jtype) .eq. ' O' .or.
     &        twh_get_elementname(jtype) .eq. ' S' .or.
     &        twh_get_elementname(jtype) .eq. 'Se' .or.
     &        twh_get_elementname(jtype) .eq. 'Te' .or.
     &        twh_get_elementname(jtype) .eq. 'Po' ) then
            lgsixj = .true.
         endif
         do ktype = jtype,ntypes
            code = twh_get_ctname(ktype)
            letterk = code(3:3)
            kk = twh_getnbtype(ktype,ktype)
c           --- determine whether this is a group 6 sp3 atom
            lgsixk = .false.
            if ( twh_get_elementname(ktype) .eq. ' O' .or.
     &           twh_get_elementname(ktype) .eq. ' S' .or.
     &           twh_get_elementname(ktype) .eq. 'Se' .or.
     &           twh_get_elementname(ktype) .eq. 'Te' .or.
     &           twh_get_elementname(ktype) .eq. 'Po' ) then
               lgsixk = .true.
            endif
            lfound = .true.
            if ( lgsixj .and. twh_get_bondpatt(jtype) .eq. 'sp3' .and. 
     &           lgsixk .and. twh_get_bondpatt(ktype) .eq. 'sp3' ) then
c              --- very special case of
c              For a single bond involving a pair of group 6 sp3 atoms
               type = type + 1
               tclocal(3) = 2.0d0
               tclocal(2) = -1.0d0
               if ( twh_get_elementname(jtype) .eq. ' O' ) then
                  vj = 2.0d0
               else
                  vj = 6.8d0
               endif
               if ( twh_get_elementname(ktype) .eq. ' O' ) then
                  vk = 2.0d0
               else
                  vk = 6.8d0
               endif
               tclocal(1) = (0.5d0)*CON_KCALTOK*dsqrt(vj*vk)
            elseif ( lgsixj .and. twh_get_bondpatt(jtype) .eq. 'sp3'
     &              .and. .not. lgsixk 
     &              .and. twh_get_bondpatt(ktype) .eq. 'sp2' ) then
c              --- very special case of
c              For a single bond involving a sp3 atom of the oxygen 
c              column and an sp2 or resonant atom of another column
               type = type + 1
               tclocal(3) = 2.0d0
               tclocal(2) = -1.0d0
               bondorder = 1.0d0
               tclocal(1) = (0.5d0)*CON_KCALTOK*(5.0d0)
     &              *dsqrt(twh_get_nbcoeff(12,jj)
     &              *twh_get_nbcoeff(12,kk))
     &              *(1.0d0 + (4.18d0)*dlog(bondorder))
            elseif ( twh_get_bondpatt(jtype) .eq. 'sp3' .and. 
     &           twh_get_bondpatt(ktype) .eq. 'sp3' ) then
c              --- specific general case of
c              (a) j = an sp3 center and k = an sp3 center
               type = type + 1
               tclocal(3) = 3.0d0
               tclocal(2) = -1.0d0
               tclocal(1) = (0.5d0)*CON_KCALTOK
     &              *dsqrt(twh_get_nbcoeff(12,jj)
     &              *twh_get_nbcoeff(12,kk))
            elseif ( (twh_get_bondpatt(jtype) .eq. 'sp3' .and. 
     &              twh_get_bondpatt(ktype) .eq. 'sp2' ) .or. 
     &              ( twh_get_bondpatt(jtype) .eq. 'sp2' .and. 
     &              twh_get_bondpatt(ktype) .eq. 'sp3' ) ) then
c              --- specific general case of
c              (b) j = an sp2 center and k = an sp3 center
               type = type + 1
               tclocal(3) = 6.0d0
               tclocal(2) = 1.0d0
               tclocal(1) = (0.5d0)*CON_KCALTOK*(1.0d0)
            elseif ( twh_get_bondpatt(jtype) .eq. 'sp2' .and.
     &              twh_get_bondpatt(ktype) .eq. 'sp2' ) then
c              --- specific general case of
c              (c) j = an sp2 center and k = and sp 2 center of variable
c              bond order
               type = type + 1
               tclocal(3) = 2.0d0
               tclocal(2) = 1.0d0
c              --- guess at the bond order
               if ( (twh_get_ctname(jtype) .eq. 'C_R' 
     &              .and. twh_get_ctname(ktype) .eq. 'N_R') .or.
     &              (twh_get_ctname(jtype) .eq. 'N_R' 
     &              .and. twh_get_ctname(ktype) .eq. 'C_R') ) then
c                 --- special bondorder
                  bondorder = 1.41d0
               elseif ( letterj .eq. 'R' .and. letterk .eq. 'R' ) then
c                 --- resonant, use 1.5
                  bondorder = 1.5d0
               else
c                 --- boldly assume the bond order is 2, not always true
c                 --- but not sure what else to do
                  bondorder = 2.0d0
               endif
               tclocal(1) = (0.5d0)*CON_KCALTOK*(5.0d0)
     &              *dsqrt(twh_get_nbcoeff(12,jj)
     &              *twh_get_nbcoeff(12,kk))
     &              *(1.0d0 + (4.18d0)*dlog(bondorder))
            else
               lfound = .false.
            endif
            if ( lfound ) then
               call twh_torstyle(GLB_SET,type,14)
               tclocal(-1) = 1.0d0
               call twh_loftor(GLB_SET,type,.true.)
               call twh_set_torcoeff(type,tclocal) 
               call twh_torff(GLB_SET,type,ffname)
               ntslocal = 1
               tnlocal(1,1) = 'wild'
               tnlocal(1,2) = twh_get_ctname(jtype)
               tnlocal(1,3) = twh_get_ctname(ktype)
               tnlocal(1,4) = 'wild'
               call twh_set_tornames(type,ntslocal,tnlocal)
            endif
         enddo
      enddo

c     --- another special case for torsions like the following
c     --- sp2 - sp2 - sp3 - wild
c     --- initialize isame to dummy value
      isame = MAXSAMET + 1
      do itype = 1,ntypes
         if ( twh_get_bondpatt(itype) .eq. 'sp2' ) then
            do jtype = 1,ntypes
               if ( twh_get_bondpatt(jtype) .eq. 'sp2' ) then
                  do ktype = 1,ntypes
                     if ( twh_get_bondpatt(ktype) .eq. 'sp3' ) then
c                       --- decide whether to increment type or isame
                        if ( isame .ge. MAXSAMET ) then
                           type = type + 1
                           isame = 0
                        endif
c                       --- found one of our special cases
c                       V0 = 2.0 kcal/mol, n = 3, phi = 180
                        tclocal(3) = 3.0d0
                        tclocal(2) = -1.0d0
                        tclocal(1) = (0.5d0)*CON_KCALTOK*(2.0d0)
                        call twh_torstyle(GLB_SET,type,14)
                        tclocal(-1) = 1.0d0
                        call twh_loftor(GLB_SET,type,.true.)
                        call twh_set_torcoeff(type,tclocal)
                        call twh_torff(GLB_SET,type,ffname)
                        isame = isame + 1
                        tnlocal(isame,1) = twh_get_ctname(itype)
                        tnlocal(isame,2) = twh_get_ctname(jtype)
                        tnlocal(isame,3) = twh_get_ctname(ktype)
                        tnlocal(isame,4) = 'wild'
c                       --- set manually
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

c     --- now add in a null torsion for everything else (a dreaded 
c     --- quadruple wildcard torsion)
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

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

c     ******************************************************************
c     * impropers                                                      *
c     ******************************************************************

c     --- list the specific cases with non-zero impropers

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,6)
      iclocal(0) =  1.0d0
      iclocal(1) =  -1.0d0
      iclocal(2) =  0.0d0
      iclocal(3) =  (50.0d0)*(CON_KCALTOK)/(3.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- C_R - O_2 - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'C_R'
      ninlocal(isame,2) = 'O_2'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- C_R - wild - O_2 - wild
      isame = isame + 1
      ninlocal(isame,1) = 'C_R'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'O_2'
      ninlocal(isame,4) = 'wild'
c     --- C_R - wild - wild - O_2
      isame = isame + 1
      ninlocal(isame,1) = 'C_R'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'O_2'
c     --- C_2 - O_2 - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'C_2'
      ninlocal(isame,2) = 'O_2'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- C_2 - wild - O_2 - wild
      isame = isame + 1
      ninlocal(isame,1) = 'C_2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'O_2'
      ninlocal(isame,4) = 'wild'
c     --- C_2 - wild - wild - O_2
      isame = isame + 1
      ninlocal(isame,1) = 'C_2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'O_2'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,6)
      iclocal(0) =  1.0d0
      iclocal(1) =  -1.0d0
      iclocal(2) =  0.0d0
      iclocal(3) =  (6.0d0)*(CON_KCALTOK)/(3.0d0)
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- C_R - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'C_R'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- C_2 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'C_2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- N_2 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'N_2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- N_R - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'N_R'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- O_2 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'O_2'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- O_R - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'O_R'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,6)
c     --- position of minimum energy
      phi = 84.4339d0
c     --- convert to radians
      phi = phi*onepi/180.0d0
c     --- arbitrarily set c_2 to 1.0d0
      iclocal(2) = 1.0d0
c     --- c_1 = -4 C2 Cos(phi)
      iclocal(1) = -4.0d0*dcos(phi)
c     --- set c_0 so the minimum energy is zero
      iclocal(0) =  -(iclocal(1)*dcos(phi)
     &     + iclocal(2)*dcos(2.0d0*phi))
c     --- force constant is 22 kcal/mol (divided by 3) when phi=0
      iclocal(3) =  CON_KCALTOK*((22.0d0)/(3.0d0)) / 
     &     (iclocal(0)+iclocal(1)+iclocal(2))
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- P_3+3 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'P_3+3'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,6)
c     --- position of minimum energy
      phi = 86.9735d0
c     --- convert to radians
      phi = phi*onepi/180.0d0
c     --- arbitrarily set c_2 to 1.0d0
      iclocal(2) = 1.0d0
c     --- c_1 = -4 C2 Cos(phi)
      iclocal(1) = -4.0d0*dcos(phi)
c     --- set c_0 so the minimum energy is zero
      iclocal(0) =  -(iclocal(1)*dcos(phi)
     &     + iclocal(2)*dcos(2.0d0*phi))
c     --- force constant is 22 kcal/mol (divided by 3) when phi=0
      iclocal(3) =  CON_KCALTOK*((22.0d0)/(3.0d0)) / 
     &     (iclocal(0)+iclocal(1)+iclocal(2))
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- As3+3 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'As3+3'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,6)
c     --- position of minimum energy
      phi = 87.7047d0
c     --- convert to radians
      phi = phi*onepi/180.0d0
c     --- arbitrarily set c_2 to 1.0d0
      iclocal(2) = 1.0d0
c     --- c_1 = -4 C2 Cos(phi)
      iclocal(1) = -4.0d0*dcos(phi)
c     --- set c_0 so the minimum energy is zero
      iclocal(0) =  -(iclocal(1)*dcos(phi)
     &     + iclocal(2)*dcos(2.0d0*phi))
c     --- force constant is 22 kcal/mol (divided by 3) when phi=0
      iclocal(3) =  CON_KCALTOK*((22.0d0)/(3.0d0)) / 
     &     (iclocal(0)+iclocal(1)+iclocal(2))
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- Sb3+3 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'Sb3+3'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      type = type + 1
      call twh_impform(GLB_SET,type,4)
      call twh_impstyle(GLB_SET,type,6)
c     --- position of minimum energy
      phi = 90.0d0
c     --- convert to radians
      phi = phi*onepi/180.0d0
c     --- arbitrarily set c_2 to 1.0d0
      iclocal(2) = 1.0d0
c     --- c_1 = -4 C2 Cos(phi), manually set to zero for phi of Pi
      iclocal(1) = 0.0d0
c     --- set c_0 so the minimum energy is zero
      iclocal(0) =  -(iclocal(1)*dcos(phi)
     &     + iclocal(2)*dcos(2.0d0*phi))
c     --- force constant is 22 kcal/mol (divided by 3) when phi=0
      iclocal(3) =  CON_KCALTOK*((22.0d0)/(3.0d0)) / 
     &     (iclocal(0)+iclocal(1)+iclocal(2))
      call twh_set_impcoeff(type,iclocal)
      call twh_impff(GLB_SET,type,ffname)
      isame = 0
c     --- Bi3+3 - wild - wild - wild
      isame = isame + 1
      ninlocal(isame,1) = 'Bi3+3'
      ninlocal(isame,2) = 'wild'
      ninlocal(isame,3) = 'wild'
      ninlocal(isame,4) = 'wild'
c     --- set nimpsame
      nislocal = isame
      call twh_set_impnames(type,nislocal,ninlocal)

      nimprotypes = type
      call twh_ffcheck(lfailure,4,type)
      if ( lfailure ) return

      return
      end
