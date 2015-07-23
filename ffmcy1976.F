#include "preproc.h"
c     ******************************************************************
c     * MCCCS - Towhee: A Monte Carlo molecular simulation program     *
c     * Copyright (C) 2007-2014 Marcus G. Martin                       *
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
      subroutine twh_ffmcy1976(lfailure,ntypes,nbondtypes,nangletypes
     &     ,nbitypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the MCY          *
c     * forcefield for water                                           *
c     *                                                                *
c     * originally written in 2006 by Jianhui Li                       *
c     * last modified 07-15-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,nbitypes
c     --- local variables
      integer ii,itype,jtype,type
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

c     --- set variables for this forcefield
      ffname = 'MCY1976'
      call twh_classical_potential(GLB_SET,'Double Exponential')
      call twh_classical_mixrule(GLB_SET,'Explicit')

c     --- reference
c     [jcp1976] O. Matsuoka; E. Clementi; M. Yoshimine 
c     'CI study of the water dimer potential surface',
c     J. Chem. Phys. 64, 1351-1361 (1976)

c     -- Hydrogen for the fit denoted CI in table VIII
      itype = 1
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,itype,masslocal) 
      elelocal =  ' H'
      call twh_elementname(GLB_SET,itype,elelocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,itype,bplocal) 
      call twh_nbff(GLB_SET,itype,ffname)
      nbnlocal = 'H:CI'
      cbnlocal = 'H:CI'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(itype,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- H-H (1-1)
      jtype = 1
      ii = twh_getnbtype(itype,jtype)
c     --- denoted as a2 in the original reference
      nbclocal(1) = (0.662712d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
c     --- b2
      nbclocal(2) = 1.29982d0/(CON_BOHRRADIUS)
c     --- no attractive portion
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- H-O (1-2)
      jtype = 2
      ii = twh_getnbtype(itype,jtype)
c     --- a3
      nbclocal(1) = (2.684452d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
c     --- b3
      nbclocal(2) = 1.439787d0/(CON_BOHRRADIUS)
c     --- attractive portion is a4,b4 in original referece
      nbclocal(3) = (0.675342d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
      nbclocal(4) = 1.141494d0/(CON_BOHRRADIUS)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- H-M (1-3)
      jtype = 3
      ii = twh_getnbtype(itype,jtype)
c     --- no vdw interaction with the midpoint, only coulomb
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      
c     -- Oxygen for the fit denoted CI in Table VIII
      itype = 2
      masslocal = 15.999d0
      call twh_mass(GLB_SET,itype,masslocal) 
      elelocal = ' O'
      call twh_elementname(GLB_SET,itype,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,itype,bplocal) 
      call twh_nbff(GLB_SET,itype,ffname)
      nbnlocal = 'O:CI'
      cbnlocal = 'O:CI'
      canlocal = 'O'
      ctnlocal = 'O'
      call twh_set_names(itype,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- O-O (2-2)
      jtype = 2
      ii = twh_getnbtype(itype,jtype)
c     --- denoted as a1 in the original reference
      nbclocal(1) = (1864.271482d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
c     --- b1
      nbclocal(2) = 2.753110d0/(CON_BOHRRADIUS)
c     --- no attractive portion
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- O-M (2-3)
      jtype = 3
      ii = twh_getnbtype(itype,jtype)
c     --- no vdw interaction with the midpoint, only coulomb
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
 
c     --- charge midpoint for the CI model in Table VIII
      itype = 3
      masslocal = 0.0d0
      call twh_mass(GLB_SET,itype,masslocal) 
      elelocal =  '  '
      call twh_elementname(GLB_SET,itype,elelocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,itype,bplocal) 
      call twh_nbff(GLB_SET,itype,ffname)
      nbnlocal = 'M:CI'
      cbnlocal = 'M:CI'
      canlocal = 'M'
      ctnlocal = 'M'
      call twh_set_names(itype,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- M-M (3-3)
      jtype = 3
      ii = twh_getnbtype(itype,jtype)
c     --- no vdw interaction with the midpoint, only coulomb
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

c     -- Hydrogen for the fit denoted inter in table VIII
      itype = 4
      masslocal = 1.0079d0
      call twh_mass(GLB_SET,itype,masslocal) 
      elelocal =  ' H'
      call twh_elementname(GLB_SET,itype,elelocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,itype,bplocal) 
      call twh_nbff(GLB_SET,itype,ffname)
      nbnlocal = 'H:inter'
      cbnlocal = 'H:inter'
      canlocal = 'H'
      ctnlocal = 'H'
      call twh_set_names(itype,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- H-H (4-4)
      jtype = 4
      ii = twh_getnbtype(itype,jtype)
c     --- denoted as a2 in the original reference
      nbclocal(1) = (1.061887d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
c     --- b2
      nbclocal(2) = 1.460975d0/(CON_BOHRRADIUS)
c     --- no attractive portion
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- H-O (4-5)
      jtype = 5
      ii = twh_getnbtype(itype,jtype)
c     --- a3
      nbclocal(1) = (2.319395d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
c     --- b3
      nbclocal(2) = 1.567367d0/(CON_BOHRRADIUS)
c     --- attractive portion is a4,b4 in original referece
      nbclocal(3) = (0.436006d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
      nbclocal(4) = 1.181792d0/(CON_BOHRRADIUS)
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- H-M (4-6)
      jtype = 6
      ii = twh_getnbtype(itype,jtype)
c     --- no vdw interaction with the midpoint, only coulomb
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
      
c     -- Oxygen for the fit denoted inter in Table VIII
      itype = 5
      masslocal = 15.999d0
      call twh_mass(GLB_SET,itype,masslocal) 
      elelocal = ' O'
      call twh_elementname(GLB_SET,itype,elelocal) 
      bplocal = 'sp3'
      call twh_bondpatt(GLB_SET,itype,bplocal) 
      call twh_nbff(GLB_SET,itype,ffname)
      nbnlocal = 'O:inter'
      cbnlocal = 'O:inter'
      canlocal = 'O'
      ctnlocal = 'O'
      call twh_set_names(itype,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- O-O (5-5)
      jtype = 5
      ii = twh_getnbtype(itype,jtype)
c     --- denoted as a1 in the original reference
      nbclocal(1) = (1734.196000d0)*(CON_HARTREETOKJ)*(CON_KJTOK)
c     --- b1
      nbclocal(2) = 2.726696d0/(CON_BOHRRADIUS)
c     --- no attractive portion
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
c     --- O-M (5-6)
      jtype = 6
      ii = twh_getnbtype(itype,jtype)
c     --- no vdw interaction with the midpoint, only coulomb
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 
 
c     --- charge midpoint for the inter model in Table VIII
      itype = 6
      masslocal = 0.0d0
      call twh_mass(GLB_SET,itype,masslocal) 
      elelocal =  '  '
      call twh_elementname(GLB_SET,itype,elelocal) 
      bplocal = 's'
      call twh_bondpatt(GLB_SET,itype,bplocal) 
      call twh_nbff(GLB_SET,itype,ffname)
      nbnlocal = 'M:inter'
      cbnlocal = 'M:inter'
      canlocal = 'M'
      ctnlocal = 'M'
      call twh_set_names(itype,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- M-M (6-6)
      jtype = 6
      ii = twh_getnbtype(itype,jtype)
c     --- no vdw interaction with the midpoint, only coulomb
      nbclocal(1) = 0.0d0
      nbclocal(2) = 0.0d0
      nbclocal(3) = 0.0d0
      nbclocal(4) = 0.0d0
      call twh_set_nbcoeff(ii,nbclocal) 

      ntypes = itype
      call twh_ffcheck(lfailure,0,itype)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed bond lengths                                             *
c     ******************************************************************

      type = 0
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = 0.9572d0
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     --- CI version O-H [jcp1976]
      ii = ii + 1
      vnlocal(ii,1) = 'O:CI'
      vnlocal(ii,2) = 'H:CI'
c     --- inter version O-H [jcp1976]
      ii = ii + 1
      vnlocal(ii,1) = 'O:inter'
      vnlocal(ii,2) = 'H:inter'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.487741d0)*(CON_BOHRRADIUS)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     OW-M distance between the water oxygen and the midpoint charge
c     site located in the direction of the bisector of the HW
c     [jcp1976]
      ii = ii + 1
      vnlocal(ii,1) = 'O:CI'
      vnlocal(ii,2) = 'M:CI'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_FIXED)
      vclocal(0) = (0.505783d0)*(CON_BOHRRADIUS)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      ii = 0
c     OW-M distance between the water oxygen and the midpoint charge
c     site located in the direction of the bisector of the HW
c     [jcp1976]
      ii = ii + 1
      vnlocal(ii,1) = 'O:inter'
      vnlocal(ii,2) = 'M:inter'
c     --- set nvsame
      nvslocal = ii
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * fixed angles                                                   *
c     ******************************************************************

      type = 0
      type = type + 1
      call twh_anglestyle(GLB_SET,type,0)
      bclocal(0) = 104.52d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     H-O-H [jcp1976]
      ii = ii + 1
      bnlocal(ii,1) = 'H'
      bnlocal(ii,2) = 'O'
      bnlocal(ii,3) = 'H'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,0)
      bclocal(0) = 52.26d0
      bclocal(1) = 1d-5
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      ii = 0
c     M-O-H (along bisector) [jcp1976]
      ii = ii + 1
      bnlocal(ii,1) = 'M'
      bnlocal(ii,2) = 'O'
      bnlocal(ii,3) = 'H'
c     --- set nbsame
      nbslocal = ii
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond increment                                                 *
c     ******************************************************************

      type = 0

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = dsqrt(0.565117d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - O for the CI model
      qbinlocal(1) = 'H:CI'
      qbinlocal(2) = 'O:CI'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -2.0d0*dsqrt(0.565117d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - O for the CI model
      qbinlocal(1) = 'M:CI'
      qbinlocal(2) = 'O:CI'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = dsqrt(0.514783d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - O for the CI model
      qbinlocal(1) = 'H:inter'
      qbinlocal(2) = 'O:inter'
      call twh_set_qbinames(type,qbinlocal)

      type = type + 1
      call twh_qbiff(GLB_SET,type,ffname)
      qbivlocal = -2.0d0*dsqrt(0.514783d0)
      call twh_qbivalue(GLB_SET,type,qbivlocal)
c     H - O for the CI model
      qbinlocal(1) = 'M:inter'
      qbinlocal(2) = 'O:inter'
      call twh_set_qbinames(type,qbinlocal)

      nbitypes = type
      call twh_ffcheck(lfailure,6,type)
      if ( lfailure ) return

      return
      end


