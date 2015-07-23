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
      subroutine twh_ffmm2(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the MM2          *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 12-20-2004 by M.G. Martin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer type,isame,j,ij
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
      double precision sumrstar,epsilon
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

      ffname = 'MM2'
      call twh_classical_potential(GLB_SET,'Exponential-6')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     ******************************************************************
c     * MM2                                                            *
c     ******************************************************************
c     References:
c     [jacs1977] 
c     N.L. Allinger; "Conformational Analysis. 130. MM2. A Hydrocarbon
c     Force Field Utilizing V1 and V2 Torsional Terms",
c     J. Am. Chem. Soc. 99, 8127-8134 (1977).

      type = 1
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     --- general carbon bonded to 2 or more other carbons
      nbnlocal = 'Ccc**(sp3)'
      cbnlocal = 'C'
      canlocal = nbnlocal
      ctnlocal = 'C'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of C(gen) with C(gen) (1-1) [jacs1977]
      j = 1
      ij = twh_getnbtype(type,j)
      sumrstar = 3.80d0
      epsilon = 0.044d0*CON_KCALTOK
      nbclocal(1) = -(2.25d0)*epsilon*(sumrstar**6)
      nbclocal(2) = (2.90d5)*epsilon
      nbclocal(3) = -(12.50d0)/sumrstar
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of C(gen) with C(methyl) (1-2) [jacs1977]
      j = 2
      ij = twh_getnbtype(type,j)
      sumrstar = 3.80d0
      epsilon = 0.044d0*CON_KCALTOK
      nbclocal(1) = -(2.25d0)*epsilon*(sumrstar**6)
      nbclocal(2) = (2.90d5)*epsilon
      nbclocal(3) = -(12.50d0)/sumrstar
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of C(gen) with H (1-3) [jacs1977]
      j = 3
      ij = twh_getnbtype(type,j)
      sumrstar = 3.34d0
      epsilon = 0.046d0*CON_KCALTOK
      nbclocal(1) = -(2.25d0)*epsilon*(sumrstar**6)
      nbclocal(2) = (2.90d5)*epsilon
      nbclocal(3) = -(12.50d0)/sumrstar
      call twh_set_nbcoeff(ij,nbclocal) 

      type = 2
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
c     --- methyl carbon bonded to 1 C and 3 H
      nbnlocal = 'Cchhh(sp3)'
      cbnlocal = 'C'
      canlocal = nbnlocal
      ctnlocal = 'C'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of C(methyl) with C(methyl) (2-2) [jacs1977]
      j = 2
      ij = twh_getnbtype(type,j)
      sumrstar = 3.80d0
      epsilon = 0.044d0*CON_KCALTOK
      nbclocal(1) = -(2.25d0)*epsilon*(sumrstar**6)
      nbclocal(2) = (2.90d5)*epsilon
      nbclocal(3) = -(12.50d0)/sumrstar
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of C(methyl) with H (3-3) [jacs1977]
      j = 3
      ij = twh_getnbtype(type,j)
      sumrstar = 3.34d0
      epsilon = 0.046d0*CON_KCALTOK
      nbclocal(1) = -(2.25d0)*epsilon*(sumrstar**6)
      nbclocal(2) = (2.90d5)*epsilon
      nbclocal(3) = -(12.50d0)/sumrstar
      call twh_set_nbcoeff(ij,nbclocal) 

      type = 3
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H'
      cbnlocal = 'H'
      canlocal = nbnlocal
      ctnlocal = 'H'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of H with H (3-3) [jacs1977]
      j = 3
      ij = twh_getnbtype(type,j)
      sumrstar = 3.00d0
      epsilon = 0.047d0*CON_KCALTOK
      nbclocal(1) = -(2.25d0)*epsilon*(sumrstar**6)
      nbclocal(2) = (2.90d5)*epsilon
      nbclocal(3) = -(12.50d0)/sumrstar
      call twh_set_nbcoeff(ij,nbclocal) 
      
      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * MM2 bond lengths                                               *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_MMTWO)
      vclocal(0) = 1.523d0
      vclocal(1) = (4.4d0)*(71.94d0)*CON_KCALTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C - C [jacs1977]
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_MMTWO)
      vclocal(0) = 1.113d0
      vclocal(1) = (4.6d0)*(71.94d0)*CON_KCALTOK
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C - H [jacs1997]
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'H'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * MM2 angles - Sixth power with multiple types and Angle-Bond    *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SIXTH_AUTO)
c     --- angle terms
      bclocal(0) = 109.5d0
      bclocal(1) = (0.45d0)*(0.021914d0)*CON_KCALTOK
      bclocal(2) = 7.0d-8
c     --- angle-bond terms
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(3) = (0.12d0)*(2.51124d0)*CON_KCALTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-C-C is not a methyl group (obviously) and also not
c     cyclobutane [jacs1997]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Ccc**(sp3)'
      bnlocal(isame,3) = 'Cwild'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SIXTH_AUTO)
c     --- angle terms
      bclocal(0) = 109.4d0
      bclocal(1) = (0.36d0)*(0.021914d0)*CON_KCALTOK
      bclocal(2) = 7.0d-8
c     --- angle-bond terms
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(3) = (0.09d0)*(2.51124d0)*CON_KCALTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-C-H where center C is not a methyl group [jacs1997]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Ccc**(sp3)'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SIXTH_AUTO)
c     --- angle terms
      bclocal(0) = 110.0d0
      bclocal(1) = (0.36d0)*(0.021914d0)*CON_KCALTOK
      bclocal(2) = 7.0d-8
c     --- angle-bond terms
      call twh_lbonang(GLB_SET,type,.true.)
      bclocal(3) = (0.09d0)*(2.51124d0)*CON_KCALTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-C-H where center C is a methyl group [jacs1997]
      isame = isame + 1
      bnlocal(isame,1) = 'Cwild'
      bnlocal(isame,2) = 'Cchhh(sp3)'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SIXTH_AUTO)
c     --- angle terms
      bclocal(0) = 109.4d0
      bclocal(1) = (0.32d0)*(0.021914d0)*CON_KCALTOK
      bclocal(2) = 7.0d-8
c     --- no angle-bond terms
      call twh_lbonang(GLB_SET,type,.false.)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H-C-H where center C is not a methyl group [jacs1997]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'Ccc**(sp3)'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_SIXTH_AUTO)
c     --- angle terms
      bclocal(0) = 109.0d0
      bclocal(1) = (0.32d0)*(0.021914d0)*CON_KCALTOK
      bclocal(2) = 7.0d-8
c     --- no angle-bond terms
      call twh_lbonang(GLB_SET,type,.false.)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H-C-H where center C is a methyl group [jacs1997]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'Cchhh(sp3)'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * MM2 torsions cosine series                                     *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      tclocal(1) = (0.20d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) = (0.27d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) = (0.093d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C-C-C-C not cyclobutane [jacs1977]
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'C'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      tclocal(1) = (0.00d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) = (0.00d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) = (0.267d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     C-C-C-H [jacs1977]
      tnlocal(1,1) = 'C'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      tclocal(1) = (0.00d0)*CON_KCALTOK*(0.5d0)
      tclocal(2) = (0.00d0)*CON_KCALTOK*(0.5d0)
      tclocal(3) = (0.237d0)*CON_KCALTOK*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     H-C-C-H [jacs1977]
      tnlocal(1,1) = 'H'
      tnlocal(1,2) = 'C'
      tnlocal(1,3) = 'C'
      tnlocal(1,4) = 'H'
      call twh_set_tornames(type,ntslocal,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
