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
      subroutine twh_ffepm(lfailure,ntypes,nbondtypes,nangletypes
     &     ,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the EPM family   *
c     * of CO2 forcefields                                             *
c     *                                                                *
c     * originally written 02-08-2005 by M.G. Martin                   *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,nbitypes
c     --- local variables
      integer ii,type,isame
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

      ffname = 'EPM'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Geometric')
      type = 0

c     ******************************************************************
c     * EPM                                                            *
c     ******************************************************************
c     References:
c     [jpc1995]
c     J.G. Harris; K.H. Yung;
c     "Carbon Dioxide's Liquid-Vapor Coexistence Curve and Critical 
c     Properties As Predicted by a Simple Molecular Model",
c     J. Phys. Chem. 99, (1995) 12021-12024.

c     *** EPM-Rigid ***
      
c     C_EPM-R [jpc1995]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.785d0
      nbclocal(2) = 28.999d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C_EPM-R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     O_EPM-R [jpc1995]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.064d0
      nbclocal(2) = 82.997d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O_EPM-R'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** EPM-Flexible ***
      
c     C_EPM-F [jpc1995]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.785d0
      nbclocal(2) = 28.999d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C_EPM-F'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     O_EPM-F [jpc1995]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.064d0
      nbclocal(2) = 82.997d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O_EPM-F'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     *** EPM2 ***
      
c     C_EPM2 [jpc1995]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 2.757d0
      nbclocal(2) = 28.129d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 12.011d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' C'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C_EPM2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
      
c     O_EPM2 [jpc1995]
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = 3.033d0
      nbclocal(2) = 80.507d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      masslocal = 15.999d0
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      elelocal =  ' O'
      bplocal = 'null'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'O_EPM2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = 'none'
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed bond lengths                                             *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.161d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C_EPM-R - O_EPM-R [jpc1995]
      isame = isame + 1
      vnlocal(isame,1) = 'C_EPM-R'
      vnlocal(isame,2) = 'O_EPM-R'
c     C_EPM-F - O_EPM-F [jpc1995]
      isame = isame + 1
      vnlocal(isame,1) = 'C_EPM-F'
      vnlocal(isame,2) = 'O_EPM-F'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 1.149d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     C_EPM2 - O_EPM2 [jpc1995]
      isame = isame + 1
      vnlocal(isame,1) = 'C_EPM2'
      vnlocal(isame,2) = 'O_EPM2'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed angles                                                   *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_RIGID)
      bclocal(0) = 180.0d0
      bclocal(1) = 1.0d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O - C - O in CO2 EPM-R [jpc1995]
      isame = isame + 1
      bnlocal(isame,1) = 'O_EPM-R'
      bnlocal(isame,2) = 'C_EPM-R'
      bnlocal(isame,3) = 'O_EPM-R'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ******************************************************************
c     * harmonic angles                                                *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (0.5d0)*(1275.0d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O - C - O in CO2 EPM-F [jpc1995]
      isame = isame + 1
      bnlocal(isame,1) = 'O_EPM-F'
      bnlocal(isame,2) = 'C_EPM-F'
      bnlocal(isame,3) = 'O_EPM-F'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 180.0d0
      bclocal(1) = (0.5d0)*(1236.0d0)*CON_KJTOK
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     O - C - O in CO2 EPM2 [jpc1995]
      isame = isame + 1
      bnlocal(isame,1) = 'O_EPM2'
      bnlocal(isame,2) = 'C_EPM2'
      bnlocal(isame,3) = 'O_EPM2'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increments                                                *
c     * inferred from the charge examples in the paper                 *
c     ******************************************************************

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = (0.6645d0)*(0.5d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - O in EPM-R [jpc1995]
      qbinlocal(1) = 'C_EPM-R'
      qbinlocal(2) = 'O_EPM-R'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = (0.6645d0)*(0.5d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - O in EPM-F [jpc1995]
      qbinlocal(1) = 'C_EPM-F'
      qbinlocal(2) = 'O_EPM-F'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = (0.6512d0)*(0.5d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     C - O in EPM2 [jpc1995]
      qbinlocal(1) = 'C_EPM2'
      qbinlocal(2) = 'O_EPM2'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end
