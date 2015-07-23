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
      subroutine twh_ffjaramillo(lfailure,ntypes,nbondtypes,nangletypes
     &     ,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the jaramillo    *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 7-31-2006 by O. Yazaydin                    *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer type,j,ij,isame
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

      ffname = 'Jaramillo'
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Explicit')
      type = 0

c     --- references for jaramillo
c     [jpcb2001]
c     E. Jaramillo; C.P. Grey; S.M. Auerbach;
c     "Molecular Dynamics Studies of Hydrofluorocarbons in Faujasite-type 
c      Zeolites: Modeling Guest-Induced Cation Migration in Dry Zeolites",
c     J. Phys. Chem. B 105, 12319-12329 (2001).

c     Carbon of HFC 
      type = 1
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'C'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of C with C (1-1) [jpcb2001]
      j = 1
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 4.35d0
      nbclocal(2) = 1.691d0*(CON_EVTOK)/1000
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of C with H (1-2) [jpcb2001]
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.458d0
      nbclocal(2) = 1.669d0*(CON_EVTOK)/1000
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of C with F (1-3) [jpcb2001]
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.878d0
      nbclocal(2) = 2.246d0*(CON_EVTOK)/1000
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal)  

c     H of HFC [jpcb2001]
      type = 2
      masslocal = 1.0079d0
      elelocal = ' H'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'H'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of H with H (2-2) [jpcb2001]
      j = 2
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 2.75d0
      nbclocal(2) = 1.648d0*(CON_EVTOK)/1000
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 
c     --- cross term of H with F (2-3) [jpcb2001]
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.083d0
      nbclocal(2) = 2.217d0*(CON_EVTOK)/1000
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal)  

c     Methyl (CH3) united-atom [jpcb2004]
      type = 3
      masslocal = 18.998d0
      elelocal = ' F'
      bplocal = 's'
      call twh_qbasevalue(GLB_SET,type,0.0d0)
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal = 'F'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 
c     --- cross term of F with F (3-3) [jpcb2001]
      j = 3
      ij = twh_getnbtype(type,j)
      nbclocal(1) = 3.457d0
      nbclocal(2) = 2.982d0*(CON_EVTOK)/1000
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      call twh_set_nbcoeff(ij,nbclocal) 

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * bond lengths - simple harmonic                                 *
c     ******************************************************************

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.526d0
      vclocal(1) = (0.5d0)*(27.988d0)*(CON_EVTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- any bond between sp3 carbons
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'C'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.105d0
      vclocal(1) = (0.5d0)*(29.540d0)*(CON_EVTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- bond between carbon and hydrogen
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'H'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)

      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.363d0
      vclocal(1) = (0.5d0)*(43.016d0)*(CON_EVTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      isame = 0
c     --- bond between carbon and flourine
      isame = isame + 1
      vnlocal(isame,1) = 'C'
      vnlocal(isame,2) = 'F'
c     --- set nvsame
      nvslocal = isame
      call twh_set_vibnames(type,nvslocal,vnlocal)


      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * angles - simple harmonic                                       *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.8d0
      bclocal(1) = (0.5d0)*(8.586d0)*(CON_EVTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     F-C-C [jpcb2001]
      isame = isame + 1
      bnlocal(isame,1) = 'F'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 107.1d0
      bclocal(1) = (0.5d0)*(5.378d0)*(CON_EVTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     F-C-H [jpcb2001]
      isame = isame + 1
      bnlocal(isame,1) = 'F'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 106.4d0
      bclocal(1) = (0.5d0)*(3.426d0)*(CON_EVTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H-C-H [jpcb2001]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'H'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 110.0d0
      bclocal(1) = (0.5d0)*(3.816d0)*(CON_EVTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H-C-C [jpcb2001]
      isame = isame + 1
      bnlocal(isame,1) = 'H'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'C'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 109.5d0
      bclocal(1) = (0.5d0)*(4.042d0)*(CON_EVTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     H-C-C [jpcb2001]
      isame = isame + 1
      bnlocal(isame,1) = 'F'
      bnlocal(isame,2) = 'C'
      bnlocal(isame,3) = 'F'
c     --- set nbsame
      nbslocal = isame
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ******************************************************************
c     * more angles                                                    * 
c     * Jaramillo paper does not list an angle term for C-C-C          *
c     * I took this one from OPLS2001 CTf-CTf-CTf                      *
c     ******************************************************************

      type = type + 1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_HARMONIC)
      bclocal(0) = 112.7d0
      bclocal(1) = (58.35d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal)
      call twh_bendff(GLB_SET,type,ffname)
      isame = 0
c     C-C-C [jpcb2001]
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
c     * torsions - OPLS style Cosine Series                            *
c     ******************************************************************

      type = type + 1
      call twh_torstyle(GLB_SET,type,2)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal (-1) = (1.0d0)
      tclocal(1) =  (0.00d0)*(CON_EVTOK)*(0.5d0)
      tclocal(2) =  (0.00d0)*(CON_EVTOK)*(0.5d0)
      tclocal(3) =  (0.0617d0)*(CON_EVTOK)*(0.5d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      isame=0
c     type 1, any 1-4 that has the C-C as 2,3
c     wild - C - C - wild [jpcb2001]
      isame = isame + 1
      tnlocal(isame,1) = 'wild'
      tnlocal(isame,2) = 'C'
      tnlocal(isame,3) = 'C'
      tnlocal(isame,4) = 'wild'
c     --- set ntsame
      call twh_set_tornames(type,isame,tnlocal)

      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return

      return
      end
