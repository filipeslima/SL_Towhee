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
      subroutine twh_ffcharmm22fe(lfailure,ntypes,nbondtypes
     &     ,nangletypes,ndihedtypes)
c     ******************************************************************
c     * sets up all of the force field parameters for the charmm22fe   *
c     * forcefield                                                     *
c     *                                                                *
c     * originally written 08-01-2006 by O. Yazaydin                   *
c     * last modified 07-21-2014 by M.G. Martin                        *
c     ******************************************************************
      implicit none
#define FUNCTION_GETNBTYPE
#include "functions.h"
c     --- variables passed to/from the subroutine
      logical lfailure
      integer ntypes,nbondtypes,nangletypes,ndihedtypes
c     --- local variables
      integer ii,type
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

      ffname = 'Charmm22'
      type = 0
      call twh_classical_potential(GLB_SET,'Lennard-Jones')
      call twh_classical_mixrule(GLB_SET,'Lorentz-Berthelot')
      
c     **************************************************************
c     * Charmm22fe force field                                     *
c     **************************************************************
c     references:
c     [chen2001] I.J. Chen; D. Yin; A.D. MacKerrell Jr.; 
c     ---------- 'Combined Ab initiao/Empirical Approach for 
c     ---------- Optimization of Lennard-Jones Parameters for
c     ---------- Polar-Neutral Compounds'
c     ---------- J. Comp. Chem. 23 (2), 199-213 (2001) 
c     ---------- toppar_all22_prot_fluoro_alkanes_jan04.str file
c     ---------- which was obtained by accessing to MacKerell website
c     ---------- on 8-1-2006

c     ******************************************************************
c     * Charmm22fe nonbonded potentials                                *
c     ******************************************************************

c     *** important note for non bonden potentials*** 
c     NOTE: In the paper cited as reference for this forcefield
c     Rmin is directly reported as opposed to the case Rmin/2 in the  
c     main Charmm reference. Therefore multiplication by 2 is removed.
      
c     (CF1) Carbon in Monofluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.80d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.060d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CF1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal) 

c     (CF2) Carbon in Difluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.10d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.042d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CF2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

c     (CF3) Carbon in Trifluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (4.60d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.020d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 12.011d0
      elelocal = ' C'
      bplocal = 'sp3'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='CF3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

c     (HF1) Hydrogen in Monofluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.64d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.028d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HF1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

c     (HF2) Hydrogen in Difluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (2.60d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.030d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 1.00800d0
      elelocal = ' H'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='HF2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

c     (F1) Fluorine in Monofluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.26d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.135d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 18.99800d0
      elelocal = ' F'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F1'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

c     (F2) Fluorine in  Difluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.26d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.105d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 18.99800d0
      elelocal = ' F'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F2'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

c     (F3) Fluorine in Trifluoroethane
      type = type + 1
      ii = twh_getnbtype(type,type)
      nbclocal(1) = (3.20d0)/(2.0d0**(1.0/6.0))
      nbclocal(2) = (0.097d0)*(CON_KCALTOK)
      nbclocal(3) = nbclocal(1)
      nbclocal(4) = nbclocal(2)
      masslocal = 18.99800d0
      elelocal = ' F'
      bplocal = 's'
      call twh_set_nbcoeff(ii,nbclocal) 
      call twh_mass(GLB_SET,type,masslocal) 
      call twh_elementname(GLB_SET,type,elelocal) 
      call twh_bondpatt(GLB_SET,type,bplocal) 
      call twh_nbff(GLB_SET,type,ffname)
      nbnlocal ='F3'
      cbnlocal = nbnlocal
      canlocal = nbnlocal
      ctnlocal = nbnlocal
      call twh_set_names(type,nbnlocal,cbnlocal,canlocal,ctnlocal)

      ntypes = type
      call twh_ffcheck(lfailure,0,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm22fe bond potentials                                     *
c     ******************************************************************

c     --- F1 - CF1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3740d0
      vclocal(1) = (420.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'F1'
      vnlocal(1,2) = 'CF1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- F2 - CF2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.3530d0
      vclocal(1) = (349.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'F2'
      vnlocal(1,2) = 'CF2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- F3 - CF3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.340d0
      vclocal(1) = (265.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'F3'
      vnlocal(1,2) = 'CF3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CF1 - CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.520d0
      vclocal(1) = (170.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CF1'
      vnlocal(1,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CF2 - CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.520d0
      vclocal(1) = (198.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CF2'
      vnlocal(1,2) = 'CT3'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CF3 - CT3
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.520d0
      vclocal(1) = (250.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 3
      vnlocal(1,1) = 'CF3'
      vnlocal(1,2) = 'CT3'
c     --- CF3 - CT2
      vnlocal(2,1) = 'CF3'
      vnlocal(2,2) = 'CT2'
c     this one is added to extend the force field
c     to create HFC-227 (1,1,1,2,3,3,3)heptafluoropropane
c     it is not in the original charmm supplementary file
c     --- CF3 - CF1
      vnlocal(3,1) = 'CF3'
      vnlocal(3,2) = 'CF1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CF1 - HF1
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0828d0
      vclocal(1) = (342.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CF1'
      vnlocal(1,2) = 'HF1'
      call twh_set_vibnames(type,nvslocal,vnlocal)

c     --- CF2 - HF2
      type = type + 1
      call twh_vibstyle(GLB_SET,type,BOND_STYLE_HARMONIC)
      vclocal(0) = 1.0828d0
      vclocal(1) = (346.0d0)*(CON_KCALTOK)
      call twh_set_vibcoeffs(type,vclocal)
      call twh_vibff(GLB_SET,type,ffname)
      nvslocal = 1
      vnlocal(1,1) = 'CF2'
      vnlocal(1,2) = 'HF2'
      call twh_set_vibnames(type,nvslocal,vnlocal)

      nbondtypes = type
      call twh_ffcheck(lfailure,1,type)
      if ( lfailure ) return

c     ******************************************************************
c     * charmm22fe angles with urey-bradley                            *
c     ******************************************************************

      type = type + 1
c     --- F2 - CF2 - F2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.0d0
      bclocal(1) = (150.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.170d0
      bclocal(3) = (10.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'F2'
      bnlocal(1,2) = 'CF2'
      bnlocal(1,3) = 'F2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- F3 - CF3 - F3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 107.0d0
      bclocal(1) = (118.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.155d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'F3'
      bnlocal(1,2) = 'CF3'
      bnlocal(1,3) = 'F3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- F1 - CF1 - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.0d0
      bclocal(1) = (44.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.369d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'F1'
      bnlocal(1,2) = 'CF1'
      bnlocal(1,3) = 'CT3'
c     this one is added to extend the force field
c     to create HFC-227 (1,1,1,2,3,3,3)heptafluoropropane
c     it is not in the original charmm supplementary file
c     --- F1 - CF1 - CF3
      bnlocal(2,1) = 'F1'
      bnlocal(2,2) = 'CF1'
      bnlocal(2,3) = 'CF3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- F2 - CF2 - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.0d0
      bclocal(1) = (50.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.357d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'F2'
      bnlocal(1,2) = 'CF2'
      bnlocal(1,3) = 'CT3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- F3 - CF3 - CT3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.0d0
      bclocal(1) = (42.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.357d0
      bclocal(3) = (30.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 3
      bnlocal(1,1) = 'F3'
      bnlocal(1,2) = 'CF3'
      bnlocal(1,3) = 'CT3'
c     --- F3 - CF3 - CT2
      bnlocal(2,1) = 'F3'
      bnlocal(2,2) = 'CF3'
      bnlocal(2,3) = 'CT2'
c     this one is added to extend the force field
c     to create HFC-227 (1,1,1,2,3,3,3)heptafluoropropane
c     it is not in the original charmm supplementary file
c     --- F3 - CF3 - CF1
      bnlocal(3,1) = 'F3'
      bnlocal(3,2) = 'CF3'
      bnlocal(3,3) = 'CF1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CT3 - CF1 - HF1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.0d0
      bclocal(1) = (31.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.168d0
      bclocal(3) = (3.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 2
      bnlocal(1,1) = 'CT3'
      bnlocal(1,2) = 'CF1'
      bnlocal(1,3) = 'HF1'
c     this one is added to extend the force field
c     to create HFC-227 (1,1,1,2,3,3,3)heptafluoropropane
c     it is not in the original charmm supplementary file
      bnlocal(2,1) = 'CF3'
      bnlocal(2,2) = 'CF1'
      bnlocal(2,3) = 'HF1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- CT3 - CF2 - HF2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 112.0d0
      bclocal(1) = (32.00d0)*(CON_KCALTOK)
      bclocal(2) = 2.168d0
      bclocal(3) = (3.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CT3'
      bnlocal(1,2) = 'CF2'
      bnlocal(1,3) = 'HF2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- F1 - CF1 - HF1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.89d0
      bclocal(1) = (57.5d0)*(CON_KCALTOK)
      bclocal(2) = 1.997d0
      bclocal(3) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'F1'
      bnlocal(1,2) = 'CF1'
      bnlocal(1,3) = 'HF1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- F2 - CF2 - HF2
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.89d0
      bclocal(1) = (41.9d0)*(CON_KCALTOK)
      bclocal(2) = 1.980d0
      bclocal(3) = (5.00d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'F2'
      bnlocal(1,2) = 'CF2'
      bnlocal(1,3) = 'HF2'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      type = type + 1
c     --- HF1 - CF1 - HF1
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 108.40d0
      bclocal(1) = (35.5d0)*(CON_KCALTOK)
      bclocal(2) = 1.746d0
      bclocal(3) = (10.40d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'HF1'
      bnlocal(1,2) = 'CF1'
      bnlocal(1,3) = 'HF1'
      call twh_set_bendnames(type,nbslocal,bnlocal)

c     ******************************************************************
c     * more angle potentials                                          *
c     * to have the CF3-CF1-CF3 working I adapted                      *
c     * CT3-CT2-CT3 from Charmm22                                      *
c     ******************************************************************

      type = type + 1
c     --- CF3 - CF1 - CF3
      call twh_anglestyle(GLB_SET,type,ANG_STYLE_CHARMM_UB)
      bclocal(0) = 114.0d0
      bclocal(1) = (53.35d0)*(CON_KCALTOK)
      bclocal(2) = 2.561d0
      bclocal(3) = (8.0d0)*(CON_KCALTOK)
      call twh_set_bencoeff(type,bclocal) 
      call twh_bendff(GLB_SET,type,ffname)
      nbslocal = 1
      bnlocal(1,1) = 'CF3'
      bnlocal(1,2) = 'CF1'
      bnlocal(1,3) = 'CF3'
      call twh_set_bendnames(type,nbslocal,bnlocal)

      nangletypes = type
      call twh_ffcheck(lfailure,2,type)
      if ( lfailure ) return

c     ******************************************************************
c     * Charmm22fe standard dihedral                                   *
c     ******************************************************************

      call twh_constant_pi(GLB_GET,onepi)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1850d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 2
c     --- wild - CT3 - CF1 - wild
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'CF1'
      tnlocal(1,4) = 'wild'
c     this one is added to extend the force field
c     to create HFC-227 (1,1,1,2,3,3,3)heptafluoropropane
c     it is not in the original charmm supplementary file
c     --- wild - CT3 - CF1 - wild
      tnlocal(2,1) = 'wild'
      tnlocal(2,2) = 'CF3'
      tnlocal(2,3) = 'CF1'
      tnlocal(2,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1780d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- wild - CT3 - CF2 - wild
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'CF2'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)

      type = type + 1
      call twh_torstyle(GLB_SET,type,3)
      call twh_loftor(GLB_SET,type,.true.)
      tclocal(-1) = 1.0d0
      call twh_ntorloop(GLB_SET,type,1)
      tclocal(1) = (0.1580d0)*(CON_KCALTOK)
      tclocal(2) = 3.0d0
      tclocal(3) = (0.0d0)*(onepi/180.0d0)
      call twh_set_torcoeff(type,tclocal)
      call twh_torff(GLB_SET,type,ffname)
      ntslocal = 1
c     --- wild - CT3 - CF3 - wild
      tnlocal(1,1) = 'wild'
      tnlocal(1,2) = 'CT3'
      tnlocal(1,3) = 'CF3'
      tnlocal(1,4) = 'wild'
      call twh_set_tornames(type,ntslocal,tnlocal)
            
      ndihedtypes = type
      call twh_ffcheck(lfailure,3,type)
      if ( lfailure ) return
      
      return
      end
